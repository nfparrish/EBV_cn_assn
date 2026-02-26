# Technology Stack

**Analysis Date:** 2026-02-26

## Languages

**Primary:**
- JavaScript (CommonJS) — installer and workflow utility logic.
- Markdown — command, workflow, reference, and template assets.

**Secondary:**
- JSON — package metadata and generated planning/config files.

## Runtime

**Environment:**
- Node.js >=16.7.0 (from `get-shit-done/package.json` engines).
- CLI-first runtime; no browser runtime required for core features.

**Package Manager:**
- npm (lockfile present: `get-shit-done/package-lock.json`).

## Frameworks

**Core:**
- Vanilla Node.js scripts and CLI routing.
- Workflow-as-Markdown architecture (commands call workflow/reference/template docs).

**Testing:**
- Built-in Node test runner (`node --test tests/*.test.cjs`).
- Assertions via `node:assert`.

**Build/Dev:**
- `esbuild` (dev dependency) used by hook build script (`build:hooks`).

## Key Dependencies

**Critical:**
- Node built-ins (`fs`, `path`, `os`, `readline`, `crypto`, `child_process`) power install and tooling logic.
- `esbuild` ^0.24.0 for hook bundling.

**Infrastructure:**
- npm package distribution via `bin/install.js` (`get-shit-done-cc`).
- Optional web search integration via `gsd-tools` websearch command (Brave API if configured).

## Configuration

**Environment:**
- Install behavior is controlled by CLI flags (`--claude`, `--opencode`, `--gemini`, `--codex`, `--global`, `--local`).
- Runtime-specific environment variables supported by installer (e.g., `CODEX_HOME`, `CLAUDE_CONFIG_DIR`, `GEMINI_CONFIG_DIR`).
- Project execution state lives in `.planning/`.

**Build:**
- `get-shit-done/package.json` defines scripts and runtime requirements.
- No ESLint/Prettier config files detected in repository root.

## Platform Requirements

**Development:**
- macOS/Linux/Windows compatible with Node.js and npm installed.
- Git recommended for commit-capable GSD workflows.

**Production / Distribution:**
- Delivered as npm package (`get-shit-done-cc`).
- Installs into runtime config directories (`~/.codex`, `~/.claude`, etc.) or local project directories.

---
*Stack analysis: 2026-02-26*
*Update after major dependency or runtime changes*
