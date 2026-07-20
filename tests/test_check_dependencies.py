from folitools import dependencies


def test_check_reports_versions_and_paths(monkeypatch, capsys):
    """Report resolved paths and parsed versions for compatible programs."""
    monkeypatch.setattr(
        dependencies,
        "PROGRAMS",
        (
            {
                "name": "tool-a",
                "command": "tool-a",
                "version_args": ("--version",),
                "checks": ((("--help",), ("--needed",)),),
                "required": True,
            },
        ),
    )
    monkeypatch.setattr(
        dependencies.shutil, "which", lambda command: f"/tools/{command}"
    )
    monkeypatch.setattr(
        dependencies,
        "_command_output",
        lambda path, args: "tool-a version 3.4.5\n--needed",
    )

    assert dependencies.check_dependencies()
    output = capsys.readouterr().out
    assert "OK tool-a 3.4.5 /tools/tool-a" in output


def test_check_fails_for_missing_required_program(monkeypatch, capsys):
    """Return failure when a required command cannot be resolved."""
    monkeypatch.setattr(
        dependencies,
        "PROGRAMS",
        (
            {
                "name": "required-tool",
                "command": "required-tool",
                "version_args": ("--version",),
                "checks": (),
                "required": True,
            },
        ),
    )
    monkeypatch.setattr(dependencies.shutil, "which", lambda command: None)

    assert not dependencies.check_dependencies()
    assert "MISSING required-tool" in capsys.readouterr().out


def test_check_allows_missing_optional_program(monkeypatch, capsys):
    """Report but tolerate an unavailable optional command."""
    monkeypatch.setattr(
        dependencies,
        "PROGRAMS",
        (
            {
                "name": "optional-tool",
                "command": "optional-tool",
                "version_args": ("--version",),
                "checks": (),
                "required": False,
            },
        ),
    )
    monkeypatch.setattr(dependencies.shutil, "which", lambda command: None)

    assert dependencies.check_dependencies()
    assert "OPTIONAL-MISSING optional-tool" in capsys.readouterr().out


def test_check_fails_for_missing_capability(monkeypatch, capsys):
    """Return failure when a resolved command lacks an option Folitools uses."""
    monkeypatch.setattr(
        dependencies,
        "PROGRAMS",
        (
            {
                "name": "old-tool",
                "command": "old-tool",
                "version_args": ("--version",),
                "checks": ((("--help",), ("--needed",)),),
                "required": True,
            },
        ),
    )
    monkeypatch.setattr(
        dependencies.shutil, "which", lambda command: f"/tools/{command}"
    )
    monkeypatch.setattr(
        dependencies,
        "_command_output",
        lambda path, args: "old-tool 1.0.0",
    )

    assert not dependencies.check_dependencies()
    assert "INCOMPATIBLE old-tool 1.0.0" in capsys.readouterr().out
