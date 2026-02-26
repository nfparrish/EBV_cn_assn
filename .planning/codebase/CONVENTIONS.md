# Coding Conventions

**Analysis Date:** 2026-02-26

## Naming Patterns

**Files:**
- kebab-case for workflow/command Markdown files (`plan-phase.md`, `map-codebase.md`).
- `.cjs` for Node CommonJS modules.
- `.test.cjs` for tests in `tests/`.

**Functions/Variables:**
- camelCase for JS function and variable names (`runGsdTools`, `createTempProject`, `getGlobalDir`).
- UPPER_SNAKE_CASE mostly for constants in local scopes is limited; stylistic emphasis is clear descriptive names.

## Code Style

**Formatting (observed):**
- 2-space indentation in JavaScript.
- Semicolons used consistently.
- Single quotes preferred in JS strings.
- Comments are practical and often sectioned for readability in large scripts.

**Linting:**
- No explicit ESLint/Prettier config files detected at repository root.
- Style consistency appears maintained manually and via review/testing.

## Import / Module Organization

**Pattern:**
1. Built-in modules first (`fs`, `path`, `os`, etc.).
2. Local module imports after built-ins.
3. Helper utilities grouped in `bin/lib/*.cjs` by domain (state, phase, roadmap, verify).

**Module Boundaries:**
- `gsd-tools.cjs` acts as command router.
- Domain logic stays in dedicated lib modules rather than monolithic inline handlers.

## Error Handling

**Patterns:**
- Fail fast for invalid args and missing resources.
- CLI usage errors produce clear messages and non-zero exits.
- Commands that read project files return explicit error states when files are missing.

**Testing Relation:**
- Many error cases are explicitly tested (e.g., missing roadmap/state, invalid cwd).

## Logging and Output

**Installer UX:**
- Uses colorized output and success markers for install progress.
- Includes structured summaries (runtime target, created assets).

**Tool Output:**
- JSON output for machine-readable workflow integration.
- Human-readable output for direct CLI usage.

## Comments and Docs

**Commenting style:**
- Prefers intent-focused comments around routing blocks and complex argument parsing.
- Major files include comprehensive usage headers documenting supported commands.

**Documentation as code:**
- Behavior contracts are encoded in workflow and template Markdown documents.
- Changes to behavior often require paired updates in command docs + workflows + templates.

## Function and Module Design

**Design bias:**
- Small domain modules for mutable project state operations.
- Utility-first approach for shared file parsing and mutation tasks.
- Deterministic commands intended for repeated orchestration by AI agents.

**Public surface:**
- Primary external API is CLI commands and workflow file conventions, not imported JS library APIs.

---
*Convention analysis: 2026-02-26*
*Update when style/tooling norms change*
