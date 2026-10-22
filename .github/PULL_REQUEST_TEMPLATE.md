## Summary
<!-- What changed and why. Name the workflows, scripts or images affected. -->


## Checks
Tick what applies to this change; see [CONTRIBUTING.md](../CONTRIBUTING.md) for the commands.
- [ ] `womtool validate` passes on every changed `.wdl` file, and `python .github/scripts/check_wdl_style.py` passes
- [ ] New or changed workflows have a `meta` block and a `parameter_meta` entry for every input and output; `docs/workflows.md` is not edited by hand
- [ ] A new directly-run workflow has a `.dockstore.yml` entry, and the entry lists `main` and no leftover feature branch
- [ ] `ruff check` and `ruff format --check` pass on `scripts`, `.github/scripts` and `tests`, and `pytest tests/unit` passes
- [ ] The WDL smoke tests pass, if a workflow they cover or `Dockerfile.utils` changed
- [ ] `shellcheck` and `hadolint` pass, if a shell script or Dockerfile changed
- [ ] The affected images build with `dockerfiles/build_docker.sh --no-push <image>`, if a Dockerfile or `scripts/` changed
- [ ] `python .github/scripts/check_markdown_style.py` passes, and `docs/` is updated for anything user-visible
- [ ] `python .github/scripts/check_dockstore_sync.py` and `python .github/scripts/check_agents_sync.py` pass, if `.dockstore.yml` or `AGENTS.md` changed
