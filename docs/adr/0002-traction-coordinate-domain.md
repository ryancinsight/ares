# 0002 — A traction boundary keeps the nodes it validated

- **Status:** Accepted
- **Date:** 2026-11-13
- **Driver:** board item `ARES-TRACTION-COORDS`, delivered by PR #78
- **Class:** `[major]` — public signature change, workspace 0.2.0 to 0.3.0
  under Cargo's leftmost-nonzero rule.

## Context

`TractionBoundary::try_new(facets, nodes)` checks every facet against
`nodes`: each index is in range, each facet measure is finite and positive.
Until 0.2.0 the boundary kept only `facets`, and
`add_consistent_loads(&self, nodes, loads)` took a second node slice. Nothing
tied that slice to the one validation read, so the validation did not cover
the geometry assembly measured:

- A replacement slice shorter than a facet index passed the shape check,
  which used the replacement's `nodes.len() * D`, then panicked on
  `nodes[*node]`.
- A replacement whose facet nodes coincide gave facet measure 0, returned
  `Ok`, and added no load: a prescribed traction silently vanished, which is
  the failure `DegenerateFacet` exists to reject.

## Decision

The boundary borrows the validated nodes for its own lifetime:

```rust
pub struct TractionBoundary<'bc, T, const D: usize> {
    facets: &'bc [TractionFacet<T, D>],
    nodes: &'bc [[T; D]],
}

pub fn try_new(facets: &'bc [TractionFacet<T, D>], nodes: &'bc [[T; D]])
    -> Result<Self, InvalidBoundary>;
pub fn add_consistent_loads(&self, loads: &mut [T]) -> Result<(), MisshapedLoad>;
```

Assembly measures facets on `self.nodes`, and the shared borrow forbids
mutating them while the boundary lives. A caller with moved nodes validates a
new boundary. Every in-repo caller migrates in the same change; no method
keeps the two-slice form.

## Alternatives rejected

**Re-validate inside `add_consistent_loads`.** Keep the second slice and rerun
the range and measure checks on it. Rejected: it duplicates `try_new`, adds a
second error type to assembly, and leaves the API admitting a geometry other
than the one the caller validated.

**Store a fingerprint of the nodes and compare.** Rejected: a length or hash
check detects some substitutions at run time; the borrow makes every one a
compile error at no run-time cost.

**Copy the nodes into the boundary.** Rejected: it allocates in an
allocation-free `no_std` core, and a copy can drift from the mesh the caller
then mutates, which is the defect in another form.

## Consequences

- Breaking: `add_consistent_loads(&nodes, &mut loads)` becomes
  `add_consistent_loads(&mut loads)`, and the nodes passed to `try_new` must
  outlive the boundary. The workspace moves to 0.3.0.
- `ares-coupling`'s partition already validated and assembled on its own
  mesh's nodes, so its migration is the argument removal.
- Rustdoc `compile_fail` examples on `TractionBoundary` show that a foreign
  slice and node mutation during the borrow do not compile. Generic oracles
  at `f32` and `f64` assert the replacement cases by value.

## Evidence that would overturn this

A caller that must load one validated facet set onto several node arrays,
for example a moving mesh, where re-validation per step is measured to
dominate. That would justify a separate type, not a second slice here.
