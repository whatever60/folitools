from pathlib import Path
import subprocess
import sys

from cyclopts import App

SCRIPT_DIR = Path(__file__).parent / "scripts"
app = App(help="Foli Tools CLI")

def run(script_name: str, args: tuple[str, ...]) -> None:
    script_path = SCRIPT_DIR / script_name
    if not script_path.exists():
        sys.exit(f"ERROR: script not found: {script_path}")
    subprocess.run([str(script_path), *args], check=True)

@app.command(help="Run fastp preprocessing")
def fastp(*args: str):
    """Arguments passed to foli_01_fastp.sh"""
    run("foli_01_fastp.sh", args)

@app.command(help="Run cutadapt demultiplexing")
def cutadapt(*args: str):
    """Arguments passed to foli_02_cutadapt.sh"""
    run("foli_02_cutadapt.sh", args)

@app.command(help="Run genome mapping")
def map(*args: str):
    """Arguments passed to foli_03_map.sh"""
    run("foli_03_map.sh", args)

@app.command(help="Run read counting")
def count(*args: str):
    """Arguments passed to foli_04_count.sh"""
    run("foli_04_count.sh", args)

if __name__ == "__main__":
    app()
