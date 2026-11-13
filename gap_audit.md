# Ares Gap Audit

Risks nobody is fixing in this repository. A planned fix belongs in
[backlog.md](backlog.md).

<a id="ARES-COUPLING-PUBLISH-CLOSURE"></a>
## ARES-COUPLING-PUBLISH-CLOSURE — publishable crate depends on an unpublishable one
- Risk: `crates/ares-coupling/Cargo.toml` has `publish.workspace = true` (workspace `publish = true`) and depends on `harmonia`, whose manifest at harmonia `origin/main` has `publish = false`; the publishable closure is not dependency-closed (engineering_gates: publish pipelines).
- Evidence: both manifests read at those revisions; ADR 0001 revision note of 2026-09-04 records the same defect. Atlas `backlog/sibling-named-crates.md` lists `ares-coupling` as `publish=false`, which contradicts the manifest at ares `ab7710580497f09222f3bff6b7518df2e69259e8`.
- Re-open trigger: `harmonia` becomes publishable, or the owner sets `publish = false` on `ares-coupling`. Owning item: atlas `backlog/sibling-named-crates.md`. The crates.io name `harmonia` (0.1.0) belongs to github.com/sogh/harmonia, so publishing it also needs another registry name.
