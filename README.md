# foldflute

Folded low-D whistle design notes and tooling.

## MCP Setup

This project uses `build123d-mcp` as an AI/MCP server via `uv tool run`; it is
not vendored or tracked as a submodule.

Project-scoped MCP config lives in `.mcp.json` and runs:

```sh
uv tool run --python 3.12 \
  --with openwind \
  --with scipy \
  build123d-mcp@latest \
  --allow-imports openwind,scipy,numpy
```

Openwind and SciPy are included for the acoustic helper workflow described in
`SPEC.md`.

## Layout

- `spec/` — InstrumentSpec JSONC documents (seed spec from SPEC.md §5)
- `foldflute/prelude.py` — acoustic engine, hole optimizer, bend corrections,
  fold-path solver; written to run inside the build123d-mcp `execute()`
  sandbox (load the spec host-side, inject as a dict)
- `scripts/` — parametric CAD builds for the solved designs + raw session log
- `reports/` — REPORT.md (first end-to-end run), solved specs, renders
- `artifacts/` — gated STL/STEP exports (M1 straight, M2 folded)

See `reports/REPORT.md` for the current design state: both M1 (straight) and
M2 (folded foot U-bend) are solved and export print-ready parts; the low A
hole layout is formally infeasible and was relaxed per the spec policy.
