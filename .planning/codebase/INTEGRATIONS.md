# External Integrations

**Analysis Date:** 2026-02-26

## APIs & External Services

**Package Distribution:**
- npm Registry — package publishing and install (`get-shit-done-cc`).
  - Integration method: npm/npx CLI.
  - Auth: Maintainer npm credentials (not stored in repo).

**Source Hosting / Collaboration:**
- GitHub — repository hosting, issues, PR workflow.
  - Integration method: git remote + optional GitHub workflows.
  - Auth: local git credentials/token/SSH.

**Optional Web Search:**
- Brave Search API (optional) used by `gsd-tools websearch` when configured.
  - Integration method: CLI utility command.
  - Auth: API key expected via environment setup outside this repository.

## Data Storage

**Databases:**
- None required for core GSD tooling.

**File Storage:**
- Local filesystem is primary storage medium for all runtime assets and project state.
- Key persisted paths: runtime config dirs (`~/.codex`, `~/.claude`, etc.) and project `.planning/`.

**Caching:**
- No dedicated external cache service observed.

## Authentication & Identity

**Runtime Identity:**
- Delegated to host runtime environments (Codex/Claude/Gemini/OpenCode) and local machine auth context.

**Git Identity:**
- Commit operations rely on local git config (`user.name`, `user.email`).

## Monitoring & Observability

**Error Tracking:**
- No third-party error tracker configured in repository code.

**Logs:**
- CLI stdout/stderr output is the primary operational telemetry.

## CI/CD & Deployment

**Distribution:**
- npm publish workflow implied by package metadata and `prepublishOnly` build script.

**CI:**
- Repository contains `.github/` directory; test/lint/release specifics depend on workflow files there.

## Environment Configuration

**Development:**
- Requires Node.js and npm.
- Optional runtime-specific env vars for install targeting (`CODEX_HOME`, `CLAUDE_CONFIG_DIR`, etc.).

**Production/Use:**
- Relies on user filesystem permissions and runtime config directories.
- No external service required for baseline local planning workflows.

## Webhooks & Callbacks

**Incoming:**
- None in current codebase scope.

**Outgoing:**
- None in current codebase scope (excluding optional websearch API call pattern).

---
*Integration audit: 2026-02-26*
*Update when new external services are added*
