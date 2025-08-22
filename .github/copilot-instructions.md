---
applyTo: "**"
---

# Workspace rules

- Python: I use uv to manage the current project. Use the interpreter at `.venv` for everything.
- Test: My test framework is pytest, as you can see in the in the [test] group  of pyproject.toml. Therefore, if you think some test is worth adding to the codebase, add to the `tests` folder, but make sure you keep the test structure clean and intuitive.
- Style: do not throw test files/readme/scripts in the project root without cleaning up.
- Don’t propose `sudo` or global installs. Use `uv add`, `uv pip`.
