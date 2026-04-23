//! Core implementation of `foli_add_tags` — shared by both the standalone
//! binary (`rust/src/bin/foli_add_tags.rs`) and the Python extension module
//! (`rust/src/lib.rs`). See the module docs there for context.

use std::collections::{BTreeSet, HashSet};
use std::fs::OpenOptions;
use std::io::Write;

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use rust_htslib::bam::{
    self,
    record::{Aux, Record},
    Format, Header, Read, Reader, Writer,
};

#[derive(Parser, Debug)]
#[command(
    name = "foli_add_tags",
    about = "Add primer/UMI/gene-with-primer tags to a collated BAM from foli map.",
    long_about = None,
)]
pub struct Args {
    /// Input BAM path, or '-' for stdin.
    #[arg(short = 'i', long = "input", default_value = "-")]
    pub input: String,

    /// Output BAM path, or '-' for uncompressed BAM to stdout.
    #[arg(short = 'o', long = "output", default_value = "-")]
    pub output: String,

    /// Name of the cell barcode tag.
    #[arg(long = "cell_tag_name", default_value = "CB")]
    pub cell_tag_name: String,

    /// Value for the cell barcode tag. If absent, no CB tag is written.
    #[arg(long = "cell_tag")]
    pub cell_tag: Option<String>,

    /// Append-mode log file for non-fatal SAM-compliance warnings.
    /// If omitted, warnings are silently dropped.
    #[arg(long = "log")]
    pub log: Option<String>,
}

const UMI_LEN: usize = 6;

// 2-byte SAM aux tag constants.
const TAG_US: &[u8; 2] = b"US";
const TAG_UC: &[u8; 2] = b"UC";
const TAG_PR: &[u8; 2] = b"PR";
const TAG_XF: &[u8; 2] = b"XF";
const TAG_XT: &[u8; 2] = b"XT";
const TAG_XN: &[u8; 2] = b"XN";
const TAG_HI: &[u8; 2] = b"HI";

/// Parse `argv` (argv[0] is the program name) and run the tag-adding pipeline.
/// Returns 0 on success, non-zero on usage or runtime failure.
pub fn run_cli(argv: Vec<String>) -> i32 {
    let args = match Args::try_parse_from(argv) {
        Ok(a) => a,
        Err(e) => {
            // clap already formats help/version output nicely; let it do its thing.
            e.print().ok();
            return match e.kind() {
                clap::error::ErrorKind::DisplayHelp | clap::error::ErrorKind::DisplayVersion => 0,
                _ => 2,
            };
        }
    };
    match run(args) {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("{:#}", e);
            1
        }
    }
}

/// Programmatic entry point; prefer this from other Rust code.
pub fn run(args: Args) -> Result<()> {
    let cell_tag_name_bytes = args.cell_tag_name.as_bytes();
    if cell_tag_name_bytes.len() != 2 {
        return Err(anyhow!(
            "cell_tag_name must be exactly 2 bytes, got {:?}",
            args.cell_tag_name
        ));
    }

    let mut reader = if args.input == "-" {
        Reader::from_stdin().context("open stdin for BAM read")?
    } else {
        Reader::from_path(&args.input).with_context(|| format!("open {}", args.input))?
    };

    let header = Header::from_template(reader.header());

    let (fmt, uncompressed) = if args.output == "-" {
        (Format::Bam, true) // match Python's "wb0" (stdout = uncompressed BAM)
    } else if args.output.ends_with(".bam") {
        (Format::Bam, false)
    } else if args.output.ends_with(".sam") {
        (Format::Sam, false)
    } else {
        return Err(anyhow!(
            "Output file must end with .bam or .sam (got {})",
            args.output
        ));
    };

    let mut writer = if args.output == "-" {
        Writer::from_stdout(&header, fmt).context("open stdout BAM writer")?
    } else {
        Writer::from_path(&args.output, &header, fmt)
            .with_context(|| format!("open writer for {}", args.output))?
    };
    if uncompressed {
        writer
            .set_compression_level(bam::CompressionLevel::Uncompressed)
            .context("set stdout compression to uncompressed")?;
    }

    let mut log_fh = match args.log.as_deref() {
        Some(path) => Some(
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .with_context(|| format!("open log file {}", path))?,
        ),
        None => None,
    };

    let mut current_qname: Vec<u8> = Vec::new();
    let mut current_primers: Vec<u8> = Vec::new();
    let mut primary_batch: Vec<Record> = Vec::new();
    let mut xt_values: Vec<Vec<u8>> = Vec::new();

    // Per-QNAME counters, incremented on the primary R1 only so each read pair
    // contributes exactly once. Emitted at program end via the SUMMARY line in
    // --log, consumed by folitools.summary.summary_stats.
    let mut total_r1_count: u64 = 0;
    let mut good_umi_count: u64 = 0;
    let mut not_na_adapter_count: u64 = 0;

    let mut record = Record::new();

    loop {
        match reader.read(&mut record) {
            None => break,
            Some(Ok(())) => {}
            Some(Err(e)) => return Err(anyhow::Error::from(e).context("reading input BAM")),
        }

        let qname = record.qname().to_vec();

        if !current_qname.is_empty() && current_qname.as_slice() != qname.as_slice() {
            finalize_group(
                &mut primary_batch,
                &xt_values,
                &current_primers,
                &current_qname,
                &mut writer,
                log_fh.as_mut(),
            )?;
            xt_values.clear();
            current_qname.clear();
            current_primers.clear();
        }

        let mut parts = qname.splitn(4, |&b| b == b'_');
        let id = parts
            .next()
            .ok_or_else(|| anyhow!("empty qname"))?
            .to_vec();
        let umi1 = parts
            .next()
            .ok_or_else(|| bad_qname(&qname, "missing UMI1"))?
            .to_vec();
        let umi2 = parts
            .next()
            .ok_or_else(|| bad_qname(&qname, "missing UMI2"))?
            .to_vec();
        let primers = parts
            .next()
            .ok_or_else(|| bad_qname(&qname, "missing primers"))?
            .to_vec();

        let mut primer_split = primers.splitn(2, |&b| b == b'+');
        let primer_fwd = primer_split.next().unwrap().to_vec();
        let primer_rev = primer_split
            .next()
            .ok_or_else(|| bad_qname(&qname, "primers missing '+'"))?
            .to_vec();

        if current_qname.is_empty() {
            current_qname.extend_from_slice(&qname);
            current_primers.clone_from(&primers);
        }

        record.set_qname(&id);

        let flag = record.flags();
        let is_read1 = (flag & 0x40) != 0;
        let is_primary = (flag & 0x900) == 0;

        // Count per-QNAME stats on the primary R1 so each read pair contributes
        // once. good_umi is UMI-only; not_na_adapter additionally requires both
        // primers assigned (== criteria_ok below), making it a subset of
        // good_umi so summary_stats' monotonic-non-increasing assert holds.
        if is_read1 && is_primary {
            total_r1_count += 1;
            let umi_ok = umi1.len() == UMI_LEN
                && umi2.len() == UMI_LEN
                && !umi1.contains(&b'N')
                && !umi2.contains(&b'N');
            let adapter_ok = primer_fwd.as_slice() != b"no_adapter"
                && primer_rev.as_slice() != b"no_adapter";
            if umi_ok {
                good_umi_count += 1;
            }
            if umi_ok && adapter_ok {
                not_na_adapter_count += 1;
            }
        }

        // Custom tags (CB/US/PR/UC/XN/XT default) are written to R1 only:
        // umi_tools group/count in --paired mode only inspects R1, so mirroring
        // them onto R2 was redundant work.
        if is_read1 {
            if let Some(cv) = args.cell_tag.as_deref() {
                let _ = record.remove_aux(cell_tag_name_bytes);
                record.push_aux(cell_tag_name_bytes, Aux::String(cv))?;
            }

            let criteria_ok = umi1.len() == UMI_LEN
                && umi2.len() == UMI_LEN
                && !umi1.contains(&b'N')
                && !umi2.contains(&b'N')
                && primer_fwd.as_slice() != b"no_adapter"
                && primer_rev.as_slice() != b"no_adapter";

            let us_str = std::str::from_utf8(&umi1).context("UMI not UTF-8")?;
            let _ = record.remove_aux(TAG_US);
            record.push_aux(TAG_US, Aux::String(us_str))?;

            let pr_str = std::str::from_utf8(&primers).context("primers not UTF-8")?;
            let _ = record.remove_aux(TAG_PR);
            record.push_aux(TAG_PR, Aux::String(pr_str))?;

            if criteria_ok {
                let mut uc = Vec::with_capacity(umi1.len() + umi2.len());
                uc.extend_from_slice(&umi1);
                uc.extend_from_slice(&umi2);
                let uc_str = std::str::from_utf8(&uc).context("UC not UTF-8")?;
                let _ = record.remove_aux(TAG_UC);
                record.push_aux(TAG_UC, Aux::String(uc_str))?;
            }

            if record.aux(TAG_XT).is_err() {
                record.push_aux(TAG_XN, Aux::I32(-1))?;
                record.push_aux(TAG_XT, Aux::String("Unassigned"))?;
            }
        }

        // XT is aggregated across every alignment of this QNAME (both mates,
        // plus any secondary/supplementary) so XF reflects the full
        // feature-assignment set for the read pair. R2's XT comes straight
        // from featureCounts; fall back to "Unassigned" if absent.
        let xt_val: Vec<u8> = match record.aux(TAG_XT) {
            Ok(Aux::String(s)) => s.as_bytes().to_vec(),
            _ => b"Unassigned".to_vec(),
        };
        xt_values.push(xt_val);

        if (flag & 0x900) != 0 {
            writer.write(&record).context("write non-primary record")?;
        } else {
            primary_batch.push(record.clone());
        }
    }

    if !current_qname.is_empty() {
        finalize_group(
            &mut primary_batch,
            &xt_values,
            &current_primers,
            &current_qname,
            &mut writer,
            log_fh.as_mut(),
        )?;
    }

    // Single-line summary consumed by folitools.summary.summary_stats. Written
    // only when --log is set so we don't silently create state without user
    // intent. cell_tag=- when --cell_tag is unset so the field position stays
    // stable for parsers.
    if let Some(fh) = log_fh.as_mut() {
        let cell_tag_field = args.cell_tag.as_deref().unwrap_or("-");
        writeln!(
            fh,
            "SUMMARY cell_tag={} total_r1={} good_umi={} not_na_adapter={}",
            cell_tag_field, total_r1_count, good_umi_count, not_na_adapter_count
        )
        .context("write SUMMARY line to log")?;
    }

    Ok(())
}

fn bad_qname(qname: &[u8], what: &str) -> anyhow::Error {
    anyhow!(
        "qname {} ({})",
        what,
        String::from_utf8_lossy(qname).into_owned()
    )
}

fn finalize_group(
    primaries: &mut Vec<Record>,
    xt_values: &[Vec<u8>],
    primers: &[u8],
    qname: &[u8],
    writer: &mut Writer,
    log_fh: Option<&mut std::fs::File>,
) -> Result<()> {
    let mut r1_idxs: Vec<usize> = Vec::new();
    let mut r2_idxs: Vec<usize> = Vec::new();
    for (i, r) in primaries.iter().enumerate() {
        if (r.flags() & 0x40) != 0 {
            r1_idxs.push(i);
        } else {
            r2_idxs.push(i);
        }
    }

    if r1_idxs.len() != 1 || r2_idxs.len() != 1 {
        if let Some(fh) = log_fh {
            let dump: Vec<String> = primaries
                .iter()
                .map(|r| format!("flag=0x{:x} tid={} pos={}", r.flags(), r.tid(), r.pos()))
                .collect();
            let _ = writeln!(
                fh,
                "WARNING non-compliant primary count for {}: R1={} R2={} [{}]",
                String::from_utf8_lossy(qname),
                r1_idxs.len(),
                r2_idxs.len(),
                dump.join(", "),
            );
        }
        let mut keep: HashSet<usize> = HashSet::new();
        if let Some(&i) = r1_idxs.first() {
            keep.insert(i);
        }
        if let Some(&i) = r2_idxs.first() {
            keep.insert(i);
        }
        // The extras we're about to downgrade are the "ghost" unmapped
        // mate of a SECONDARY alignment that STAR forgot to mark
        // secondary. See STAR issue #2190
        // (https://github.com/alexdobin/STAR/issues/2190), reported
        // 2024-08 against 2.7.10b and 2.7.11b, still open and unresolved
        // as of the latest release (2.7.11b, 2024-01-26). Same bug also
        // leaves `HI:i:<n>` at 0 on these ghosts — which, after we flip
        // 0x100, looks like a multimap secondary carrying the primary
        // pair's hit index. We tried fixing at the STAR layer with
        // `--chimSegmentMin 12` but that controls chimeric-junction
        // emission, not the mate-pairing code path; verified on a
        // targeted reproduction with the exact CSF2RA+TIRAP read pairs
        // from a prior production warning (e.g.
        // LH00328:76:22G7CJLT3:4:2187:9174:22278) — 7/7 still produce
        // the anomaly even with the flag set. Only downstream patching
        // neutralizes it, which is what this block does.
        //
        // We can't reconstruct the correct HI without buffering every
        // secondary record of the QNAME (the downgraded record's PNEXT
        // points at a secondary that was already flushed to the output
        // stream). Leaving a stale `HI:i:0` on a non-primary is worse
        // than having no HI: it collides with the primary pair's HI
        // under naive readers. Strip it — RNEXT/PNEXT still uniquely
        // identifies the mate.
        for (i, r) in primaries.iter_mut().enumerate() {
            if !keep.contains(&i) {
                let f = r.flags();
                r.set_flags(f | 0x100);
                let _ = r.remove_aux(TAG_HI);
            }
        }
    }

    let mut gene_set: BTreeSet<&[u8]> = BTreeSet::new();
    for tag in xt_values {
        if tag.as_slice() != b"Unassigned" {
            for part in tag.split(|&b| b == b',') {
                if !part.is_empty() {
                    gene_set.insert(part);
                }
            }
        }
    }
    let mut xf_value: Vec<u8> = Vec::new();
    if gene_set.is_empty() {
        xf_value.extend_from_slice(b"Unassigned,");
        xf_value.extend_from_slice(primers);
    } else {
        let mut first = true;
        for g in &gene_set {
            if !first {
                xf_value.push(b',');
            }
            xf_value.extend_from_slice(g);
            first = false;
        }
        xf_value.push(b',');
        xf_value.extend_from_slice(primers);
    }
    let xf_str = std::str::from_utf8(&xf_value).context("XF contains non-UTF-8 bytes")?;

    // XF is stamped on primary R1 only; R2 is passed through untouched
    // (umi_tools group/count in --paired mode only inspects R1).
    for r in primaries.iter_mut() {
        let f = r.flags();
        let is_primary = (f & 0x900) == 0;
        let is_read1 = (f & 0x40) != 0;
        if is_primary && is_read1 {
            let _ = r.remove_aux(TAG_XF);
            r.push_aux(TAG_XF, Aux::String(xf_str))?;
        }
        writer.write(r).context("write primary record")?;
    }

    primaries.clear();
    Ok(())
}
