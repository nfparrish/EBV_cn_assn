# Architecture

**Analysis Date:** 2026-02-26

## Pattern Overview

**Overall:** Monolithic CLI installer + workflow/document-driven orchestration system.

**Key Characteristics:**
- Installer script materializes runtime-specific assets into `.codex`, `.claude`, `.gemini`, or `.opencode`.
- Operational behavior is encoded in Markdown skills/workflows/templates rather than only in JS modules.
- `gsd-tools.cjs` acts as a central utility router for state, roadmap, phase, verification, and scaffolding operations.

## Layers

**Installation Layer:**
- Purpose: Runtime selection, path resolution, install/uninstall operations.
- Contains: `get-shit-done/bin/install.js`.
- Depends on: Node FS/path APIs and package metadata.
- Used by: `npx get-shit-done-cc` invocation.

**Workflow Definition Layer:**
- Purpose: Define execution protocols and command behavior.
- Contains: `get-shit-done/commands/gsd/*.md`, `get-shit-done/get-shit-done/workflows/*.md`.
- Depends on: References and templates.
- Used by: Runtime skill systems and command invocations.

**Tooling Layer:**
- Purpose: Deterministic file/state operations and project bookkeeping.
- Contains: `get-shit-done/get-shit-done/bin/gsd-tools.cjs` + `bin/lib/*.cjs`.
- Depends on: `.planning` file structure.
- Used by: Workflows during plan/execute/progress/milestone operations.

## Data Flow

**Install + Runtime Activation Flow:**
1. User runs installer via `npx get-shit-done-cc@latest`.
2. `install.js` parses flags and resolves runtime/location.
3. Assets (skills, workflows, agents, templates) are copied into target config dir.
4. Manifest/version files are written.
5. Runtime command (`$gsd-*`, `/gsd:*`) becomes available in target environment.

**Planning Workflow Flow:**
1. User triggers a GSD skill (e.g., `$gsd-new-project`, `$gsd-plan-phase`).
2. Skill references workflow Markdown.
3. Workflow invokes `gsd-tools` commands for initialization, parsing, updates, and commits.
4. `.planning` docs are generated/updated.
5. Subsequent workflows consume these docs as state memory.

**State Management:**
- File-based and document-centric (`.planning/*`).
- No persistent service/database process.
- Each command invocation is independent and stateless beyond project files.

## Key Abstractions

**Skill/Workflow Documents:**
- Purpose: Declarative orchestration and behavioral policy.
- Examples: `new-project.md`, `plan-phase.md`, `map-codebase.md`.
- Pattern: Markdown contract with stepwise procedure and expected artifacts.

**Tool Commands (`gsd-tools`):**
- Purpose: Atomic, scriptable operations on planning state.
- Examples: `init`, `state`, `roadmap`, `phase`, `commit`, `verify`.
- Pattern: Router + domain modules in `bin/lib`.

## Entry Points

**Installer Entry:**
- Location: `get-shit-done/bin/install.js`.
- Triggers: npm/npx execution.
- Responsibilities: install/uninstall, runtime targeting, config file updates.

**Workflow Utility Entry:**
- Location: `get-shit-done/get-shit-done/bin/gsd-tools.cjs`.
- Triggers: workflow shell commands.
- Responsibilities: parse args, route subcommands, mutate/read `.planning` state.

## Error Handling

**Strategy:**
- CLI validation with explicit usage errors and non-zero exits.
- Fail-fast for invalid flags, missing paths, and malformed operations.

**Patterns:**
- Command-level argument checks before mutation.
- Test-backed behavior for missing files and malformed docs in `tests/*.test.cjs`.

## Cross-Cutting Concerns

**Logging/UX:**
- Human-readable terminal output with colorized status in installer.

**Validation:**
- Multiple `verify` and `init` commands enforce project consistency and required files.

**Git Integration:**
- `gsd-tools commit` centralizes planning-doc commits when git is available.

---
*Architecture analysis: 2026-02-26*
*Update when major workflow/tooling patterns change*
