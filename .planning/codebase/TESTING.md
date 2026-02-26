# Testing Patterns

**Analysis Date:** 2026-02-26

## Test Framework

**Runner:**
- Node built-in test runner (`node --test`).
- Invoked via `npm test` in `get-shit-done/package.json`.

**Assertion Library:**
- `node:assert` standard assertions.

**Run Commands:**
```bash
cd get-shit-done
npm test
node --test tests/*.test.cjs
```

## Test File Organization

**Location:**
- Centralized under `get-shit-done/tests/`.

**Naming:**
- Domain-oriented test files (`init.test.cjs`, `phase.test.cjs`, `roadmap.test.cjs`, `state.test.cjs`, `verify.test.cjs`, `commands.test.cjs`, `milestone.test.cjs`).

**Helpers:**
- Shared utility helpers in `tests/helpers.cjs`.
- Uses temp directories for isolated filesystem-based test scenarios.

## Test Structure

**Suite Organization:**
- `describe` blocks per command domain.
- `test` cases for both happy-path and failure-path behavior.
- Uses `beforeEach`/`afterEach` to create and clean temporary project fixtures.

**Patterns:**
- Command execution tests shell out to `gsd-tools.cjs` via helper wrappers.
- JSON outputs are parsed and validated with explicit field assertions.
- Error behavior is asserted for missing files, malformed docs, and invalid arguments.

## Mocking and Isolation

**Approach:**
- Minimal mocking; heavy use of real file I/O in temporary directories.
- Process-level command invocation is tested to verify end-to-end CLI behavior.

**What is isolated:**
- Project filesystem state per test via temp roots.
- Test artifacts removed after each case to prevent cross-test contamination.

## Coverage Focus

**Strongly covered areas:**
- `init` command outputs and field extraction.
- phase/roadmap parsing and numbering semantics.
- state/frontmatter mutations.
- command routing and utility operations.

**Likely weaker areas:**
- End-to-end runtime installs across all OS/runtime combinations.
- Interactive installer prompt branches (some are difficult to unit-test thoroughly).

## Test Types

**Unit-style CLI tests:**
- Validate individual command outputs and edge cases with synthetic files.

**Integration-style filesystem tests:**
- Validate behavior across multiple planning docs (`STATE`, `ROADMAP`, phase files).

---
*Testing analysis: 2026-02-26*
*Update when test framework or strategy changes*
