# Repository Guidelines

## Project Structure & Module Organization

This repository is currently a clean slate (only Git metadata). When adding code, use this layout to keep things predictable:

```
repo-root/
  src/       # application/library code
  tests/     # automated tests mirroring src/
  scripts/   # dev/ops utilities and one-off tools
  configs/   # app/CI configuration (yaml/json)
  assets/    # static files (images, styles)
  data/      # local data (gitignored)
  docs/      # additional docs and ADRs
```

Keep modules small and cohesive. Public entry points live under `src/`.

## Build, Test, and Development Commands

Define common tasks in a `Makefile` (preferred) or your language’s toolchain. Examples:

```
make setup   # install dependencies
make build   # compile/package the project
make test    # run tests with coverage
make fmt lint# format and lint the codebase
```

Language-specific examples (use what applies):
- Python: `uv sync` or `pip install -r requirements.txt`, then `pytest -q`.
- Node.js: `npm ci` and `npm test` (or `pnpm i` / `pnpm test`).

## Coding Style & Naming Conventions

- Indentation: 4 spaces (Python); 2 spaces (JS/TS). 
- Names: files/modules `snake_case` (Python) or `kebab-case` (scripts); classes `PascalCase`; constants `UPPER_SNAKE_CASE`.
- Formatting/Linting: adopt auto-formatters and enforce via pre-commit hooks. Python: Black + Ruff. JS/TS: Prettier + ESLint.

## Testing Guidelines

- Framework: prefer `pytest` (Python) or `Vitest/Jest` (TS/JS).
- Location: place tests in `tests/`, mirroring `src/` structure.
- Naming: `test_*.py` (pytest) or `*.spec.ts`/`*.test.ts` (TS/JS).
- Coverage: target ≥80%; run tests locally before opening a PR.

## Commit & Pull Request Guidelines

- Commits: follow Conventional Commits (`feat:`, `fix:`, `docs:`, `refactor:`, `chore:`, `test:`). Scope as needed (e.g., `feat(api): ...`).
- PRs: include a clear description, link issues (`Fixes #123`), screenshots/logs for UX/CLI changes, and note any breaking changes.
- Checklist: tests pass, docs updated, secrets excluded, small focused diff.

## Security & Configuration Tips

- Never commit secrets; use `.env` and provide `.env.example`.
- Git-ignore large or generated artifacts (`data/`, `dist/`, `*.env`, `.venv/`, `node_modules/`).
- Pin dependencies and update via dedicated PRs.

