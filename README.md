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
