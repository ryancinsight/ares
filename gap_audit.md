# Ares Gap Audit

## Verified downstream provider migration

A new external Cargo project resolved the consuming package without a pre-existing lock or provider overrides. Its registry graph includes `horae-time 0.4.0`, `proteus-mat 0.3.0`. The committed lock passes 120 executed tests covering structural verification, analytical solutions and field exchange (0 tests excluded or skipped in that run). Ares pins the migrated Git-only Harmonia source. Its Cargo dependency follows the package name on that source. Harmonia’s registry identity and graduation remain a separate decision.

## Incremental release baseline

Reviewed [ares-solid 0.3.0](crates/ares/Cargo.toml) at [45346575397b](https://github.com/ryancinsight/ares/tree/45346575397b55d6c6fb71df134498363511daa3). The registry has no entry for this package and version. First publication requires a verified package graph and a clean consumer.

The next [patch] checkpoint uses [cylinder](crates/ares-operator/tests/verification/cylinder.rs). [Target version and acceptance criteria](backlog.md#ares-stable-qualification).

## Symmetric projection [patch]

At `a75eae35`, averaging two maximum finite entries returned infinity.
A 3×3 projection now computes three distinct off-diagonal means and
copies its three diagonal entries. The previous full traversal computed
nine means. The oracle compares native precision results with the
standard midpoint operation, including ties adjacent to zero. All 121
locked workspace tests pass, as do strict all-target Clippy and the
no-default-features Solid build. The Solid compatibility check passed
196 rules and skipped 57. A standalone consumer could not find
`ares-solid = "=0.3.1"` in the registry.

Run `cargo bench --locked -p ares-solid --bench symmetric_projection`
for 15 samples of 1,048,576 binary64 3D projections. Median times on
this host were 1.439 ns per projection and 3.439 ns for the full matrix
midpoint reference. Applying the same input workload to the original
unchecked implementation and candidate in separate executables measured
1.225 and 1.456 ns respectively. The finite-range check costs about
0.23 ns per projection against the original arithmetic on these inputs.

## Athena Leto 0.5 provider identity [patch]

Updating Athena alone to `2be05424` fails the `ares-operator` trait
implementation because the Athena adapter takes crates.io Leto views and
Ares still loads Git Leto. Select the same published Leto and Leto Ops 0.52
source in both workspaces. Ares and Harmonia continue to share Athena Core
0.4 from the same Git workspace. All 120 locked tests pass on Rust 1.95,
including second-order refinement and the cantilever benchmark.

Risks nobody is fixing in this repository. A planned fix belongs in
[backlog.md](backlog.md).

<a id="ARES-COUPLING-PUBLISH-CLOSURE"></a>
## ARES-COUPLING-PUBLISH-CLOSURE — publishable crate depends on an unpublishable one
- Risk: `crates/ares-coupling/Cargo.toml` has `publish.workspace = true` (workspace `publish = true`) and depends on `harmonia`, whose manifest at harmonia `origin/main` has `publish = false`; the publishable closure is not dependency-closed (engineering_gates: publish pipelines).
- Evidence: both manifests read at those revisions; ADR 0001 revision note of 2026-09-04 records the same defect. Atlas `backlog/sibling-named-crates.md` lists `ares-coupling` as `publish=false`, which contradicts the manifest at ares `ab7710580497f09222f3bff6b7518df2e69259e8`.
- Re-open trigger: `harmonia` becomes publishable, or the owner sets `publish = false` on `ares-coupling`. Owning item: atlas `backlog/sibling-named-crates.md`. The crates.io name `harmonia` (0.1.0) belongs to github.com/sogh/harmonia, so publishing it also needs another registry name.
