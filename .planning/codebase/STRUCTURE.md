# Codebase Structure

**Analysis Date:** 2026-02-26

## Directory Layout

```
GSD/
├── .codex/                    # Installed local Codex runtime assets
│   ├── get-shit-done/         # Runtime workflows/templates/references
│   ├── skills/                # gsd-* skill directories
│   └── agents/                # Agent definitions
├── .planning/                 # Project planning state and maps
│   ├── codebase/              # Generated codebase map docs
│   ├── PROJECT.md
│   ├── REQUIREMENTS.md
│   ├── ROADMAP.md
│   ├── STATE.md
│   └── config.json
└── get-shit-done/             # Upstream source checkout
    ├── bin/                   # Installer entrypoint
    ├── commands/gsd/          # User-facing command docs
    ├── get-shit-done/         # Core workflows/templates/tooling
    │   ├── bin/               # gsd-tools CLI + lib modules
    │   ├── workflows/
    │   ├── templates/
    │   └── references/
    ├── tests/                 # Node test files
    ├── docs/                  # User/developer docs
    └── package.json
```

## Directory Purposes

**`.codex/`**
- Purpose: Active runtime installation used by Codex in this workspace.
- Contains: executable skill assets, workflow docs, templates, and agent specs.
- Key files: `.codex/get-shit-done/workflows/*.md`, `.codex/skills/*/SKILL.md`.

**`.planning/`**
- Purpose: This project’s durable state for planning/execution.
- Contains: project context, requirements, roadmap, status, and codebase map docs.
- Key files: `.planning/PROJECT.md`, `.planning/ROADMAP.md`, `.planning/STATE.md`.

**`get-shit-done/`**
- Purpose: Source repository for the GSD package itself.
- Contains: installer code, command/workflow assets, docs, tests.
- Key files: `bin/install.js`, `get-shit-done/bin/gsd-tools.cjs`, `tests/*.test.cjs`.

## Key File Locations

**Entry Points:**
- `get-shit-done/bin/install.js` — npm/npx installation entrypoint.
- `get-shit-done/get-shit-done/bin/gsd-tools.cjs` — workflow utility command router.

**Configuration:**
- `get-shit-done/package.json` — package metadata, engines, scripts.
- `.planning/config.json` — active project workflow configuration.

**Core Logic:**
- `get-shit-done/get-shit-done/bin/lib/*.cjs` — state/roadmap/phase/verify logic.
- `get-shit-done/get-shit-done/workflows/*.md` — orchestrated process definitions.

**Testing:**
- `get-shit-done/tests/*.test.cjs` — CLI and workflow utility tests.
- `get-shit-done/tests/helpers.cjs` — shared test helpers.

**Documentation:**
- `get-shit-done/README.md` — installation and usage guide.
- `get-shit-done/docs/` — extended user documentation.

## Naming Conventions

**Files:**
- kebab-case for most Markdown workflow/command files (`new-project.md`, `plan-phase.md`).
- `.test.cjs` suffix for tests.
- Uppercase for major docs (`README.md`, `CHANGELOG.md`, `LICENSE`).

**Directories:**
- kebab-case folder naming (`get-shit-done`).
- grouped-by-role under `workflows/`, `templates/`, `references/`, `skills/`.

## Where to Add New Code

**New workflow behavior:**
- Primary: `get-shit-done/get-shit-done/workflows/`.
- Supporting references/templates: matching files in `references/` or `templates/`.

**New utility command:**
- Router changes: `get-shit-done/get-shit-done/bin/gsd-tools.cjs`.
- Implementation: new/updated module in `get-shit-done/get-shit-done/bin/lib/`.
- Tests: `get-shit-done/tests/`.

**New runtime command surface:**
- Command docs: `get-shit-done/commands/gsd/`.
- Codex skill mapping: `.codex/skills/*/SKILL.md` after install.

---
*Structure analysis: 2026-02-26*
*Update when directory layout changes*
