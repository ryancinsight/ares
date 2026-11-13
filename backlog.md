# Ares Backlog

Open items in priority order. The meta item
`ATLAS-ARES-PROMOTION-2026-09-03` (atlas `backlog.md#ares-promotion`) owns
creation, registration and first publication of the repository and is not
copied here. Items carry no `milestone:` until `ARES-RELEASE-LADDER` lands a
ladder. Evidence below was read at `basis:` unless stated otherwise.

<a id="ARES-RELEASE-LADDER"></a>
## ARES-RELEASE-LADDER — Draft the release scope and milestone ladder [patch] — todo
- priority: verification; needs: none; basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `docs/adr/` (release-scope ADR; ADR 0002 is taken by the traction coordinate binding), `docs/adr/README.md`, `backlog.md`
- Outcome: release-scope ADR and milestone ladder per the stack roadmap rule (sprint: roadmap). The ladder opens this board; every item then carries `milestone:`.
- Acceptance: the ADR records each capability area of `ares-solid`, `ares-operator` and `ares-coupling` as delivered, in scope or a non-goal; each milestone exit check is a runnable check that fails without its feature; an independent judge verdict maps every capability area and open feature item to a check or non-goal.
- Inputs: `README.md` scope paragraph, atlas ADRs 0055, 0057 and 0059, `docs/book`. Phase 0 non-goals (plasticity, contact, finite deformation, dynamics, fracture) come from the meta item.
- Next step: draft the ADR from those inputs with a recommended scope, name any scope no source determines as a veto point, and file the ladder under its own `judge-required` PR.

<a id="ARES-README-CRATE-TABLE"></a>
## ARES-README-CRATE-TABLE — List every workspace crate in the README [patch] — todo
- priority: verification; needs: none; basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `README.md`
- Defect (read at basis): the `## Crates` table lists `ares-solid` and `ares-operator`; the workspace is `members = ["crates/*"]` and `crates/ares-coupling` exists, which `docs/book/README.md` and `docs/book/stack_position.md` already list as the third crate.
- Outcome: the README table names all three crates with registry names and the ADR 0001 reference.
- Acceptance: every `name` under `crates/*/Cargo.toml` appears in the README table (`cargo metadata --no-deps` names compared with the table rows); `git grep -n 'ares-coupling' README.md` returns the row.
- Next step: add the row, then run the docs and link checks of the committed gate.

<a id="ARES-PACKAGE-CLOSURE"></a>
## ARES-PACKAGE-CLOSURE — Qualify packaged builds against registry providers [patch] — blocked
- priority: verification; needs: none (external providers); basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `Cargo.toml` provider requirements, `Cargo.lock`, package qualification in CI
- Blocker (crates.io API, 2026-11-13): `proteus-mat` 0.1.1, its latest, requires `eunomia ^0.10.0` while `Cargo.toml` requires eunomia 0.12.0; `leto-ops` `max_stable_version` is 0.46.4 while `Cargo.toml` requires `>=0.47.0, <0.49.0`; `ares-solid`, `ares-operator` and `ares-coupling` are not on the registry. A normalized package build cannot resolve `leto-ops` or the unpublished Ares names from the registry; the eunomia split additionally breaks compilation (PR #73's E0277).
- Re-open trigger: a `proteus-mat` release on eunomia 0.12 and a `leto-ops` 0.47 release are published by their owners. First publication of `ares-solid` is the A9 Ask-User step of atlas `backlog/ares-promotion.md`; `ares-operator` and `ares-coupling` are outside A9 while `harmonia` is `publish = false`.
- Outcome and acceptance: `cargo package --locked` verification passes for each publishable crate against registry providers, with no provider back-pin and no reduced test gate.
- Next step: when the trigger fires, run `cargo package --locked` for each publishable crate against the registry and record the result here; the E0277 PR #73 reports is not re-run.
