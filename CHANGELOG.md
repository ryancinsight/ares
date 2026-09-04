# Changelog

## [Unreleased]

### Added

- Repository scaffold at the stack lint and gate floor: pinned toolchain,
  committed nextest budgets, `deny.toml` carrying the atlas ADR 0055 substrate
  prohibitions, pedantic clippy with `unwrap_used`, `indexing_slicing`, and
  `arithmetic_side_effects` denied, and `#![forbid(unsafe_code)]` with
  `#![deny(missing_docs)]`.

  No physics yet. The scaffold passes the full gate before any is added, so a
  later failure is attributable to the physics rather than to the floor.
