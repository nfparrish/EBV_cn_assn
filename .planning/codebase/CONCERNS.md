# Codebase Concerns

**Analysis Date:** 2026-02-26

## Tech Debt

**Large installer script (`bin/install.js`):**
- Issue: very large single-file implementation (~2000+ lines).
- Why: installer supports multiple runtimes, modes, and migration logic in one place.
- Impact: higher change risk and harder local reasoning/testing for modifications.
- Fix approach: incrementally extract runtime-specific and settings-mutation modules.

**Behavior split across docs + code:**
- Issue: operational behavior is distributed between workflow Markdown and JS utilities.
- Why: system is intentionally workflow-driven.
- Impact: easy to introduce drift if command docs and utility behavior diverge.
- Fix approach: enforce doc/code sync checks and regression tests for key workflow contracts.

## Known Bugs / Sharp Edges

**State coupling to file conventions:**
- Symptoms: command failures when expected `.planning` files/sections are malformed.
- Trigger: manual edits that break expected heading/frontmatter structures.
- Workaround: run `gsd-tools` init/verify commands and repair docs.
- Root cause: parser logic assumes convention-compliant document structure.

## Security Considerations

**Sensitive data in planning docs:**
- Risk: users may accidentally place secrets in markdown state/research docs.
- Current mitigation: map workflow includes explicit pre-commit secret scan guidance.
- Recommendations: add automated secret scanning in CI and stricter pre-commit hooks.

**Shell command orchestration:**
- Risk: workflows can execute shell commands in user environment.
- Current mitigation: runtime permission systems and user confirmation policies.
- Recommendations: use least-privilege runtime settings and avoid storing credentials in plaintext files.

## Performance Bottlenecks

**Very large Markdown parsing surface:**
- Problem: many workflows/templates and large files increase parse/search overhead for agents.
- Cause: feature growth over time with accumulated docs.
- Improvement path: archive obsolete assets, optimize init/lookup helpers, add focused indexing.

## Fragile Areas

**Phase/roadmap numbering semantics:**
- Why fragile: decimal phase insertion/removal and renumbering rules are nuanced.
- Common failures: mismatches between disk phase dirs and roadmap sections.
- Safe modification: rely on existing phase commands (`phase add/insert/remove`) and corresponding tests.
- Test coverage: good command-level coverage in `tests/phase.test.cjs`, but keep extending for edge cases.

## Scaling Limits

**Repository-scale context load:**
- Current capacity: practical for moderate-size planning repos.
- Limit: very large codebases + very large planning history can degrade workflow speed and context quality.
- Symptoms at limit: slower command execution and more brittle parsing/planning behavior.
- Scaling path: periodic archive/cleanup of old phase artifacts and concise state summaries.

## Dependencies at Risk

**Node version drift:**
- Risk: local environments below required Node baseline break installer/tools.
- Impact: installation or command failures.
- Migration plan: maintain explicit version checks and clear upgrade guidance in README.

## Test Coverage Gaps

**Interactive installer branches:**
- What's not fully tested: some interactive prompt flows and multi-runtime edge permutations.
- Risk: regressions in less-common runtime/flag combinations.
- Priority: Medium.
- Difficulty to test: requires robust prompt simulation and environment matrix setup.

---
*Concerns audit: 2026-02-26*
*Update as issues are fixed or newly discovered*
