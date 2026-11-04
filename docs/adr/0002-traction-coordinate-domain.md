# 0002 — Bind traction validation to its coordinate domain

- **Status:** Accepted
- **Class:** `[major]`, pre-1.0 release 0.3.0
- Date: 2026-11-04

## Context

`TractionBoundary::try_new` validates facet connectivity and measure against a
node slice, but retains only the facets. `add_consistent_loads` accepts another
node slice and indexes it under assumptions established for the first slice.
The boundary is documented as validated on one mesh; the type does not retain
that invariant.

At `b66cda11d160653223fc4e828269de81d5fd4476`, real public calls demonstrate:

- Construct a boundary for a unit edge with node indices `[0, 1]`, then supply
  one replacement node and a load vector sized for that replacement. The shape
  check passes and coordinate gathering panics at index 1.
- Supply two coincident replacement nodes instead. Construction previously
  accepted a unit-length edge, but assembly silently adds zero traction load.

Both witnesses execute in binary32 and binary64 with unchanged production
source. These are validation-domain failures, separate from the facet-measure
algorithm's numerical behavior.

## Decision

Store the validated node slice alongside the facet slice in
`TractionBoundary<'bc, T, D>`. Both are immutable borrows lasting for the
boundary's lifetime. Change assembly to `add_consistent_loads(&mut loads)`;
it reads the retained coordinates and checks the load shape against that
domain. Update every in-repository caller, including the coupling partition,
operator oracles and book example, in the same change.

The method cannot receive replacement geometry. Coordinates cannot be mutated
or destroyed while a live boundary still uses them. A moving or different mesh
constructs a new boundary and validates its geometry. There is no pointer-
identity check, repeated validation in the load loop, compatibility forwarding
method or copied coordinate buffer.

A facet has dimension `D - 1`, so `D >= 1` is a structural precondition. Reject
zero-dimensional instantiation at compile time rather than reaching an empty
chunk width or unsigned subtraction. Preserve the generic positive-dimension
formula; do not restrict it to named spatial dimensions.

## Invariant and consequences

Construction proves every facet index is below the retained slice length and
every accepted facet measure is finite and positive. Immutable borrowing keeps
that exact input valid. Assembly checks the caller's load length before writing.
Consequently each gathered coordinate and each nodal load slot is in bounds for
the whole operation. The arithmetic and additive load formula stay unchanged.

The boundary stores one additional slice reference: 16 bytes on a 64-bit host.
It still allocates nothing and copies no node data. Benchmark ordinary typed
assembly before/after; do not trade validation or test limits for speed.

This intentionally changes the public method signature and constructor borrow
contract. Deliver at 0.3.0 with a migration guide: remove the redundant nodes
argument, retain immutable coordinates for the boundary lifetime, and rebuild
the boundary after geometry changes. Update workspace internal requirements;
preserve external dependencies and every external lock record.

## Rejected alternatives

- Checking only replacement lengths cannot establish the validated geometry.
- Preflight revalidation on each call can be sound under a different contract,
  but needs geometry errors and repeated work; the existing contract describes
  a boundary validated on one mesh. Retaining that immutable domain encodes it
  directly without two coordinate authorities.
- Ignoring a retained legacy argument hides the API change behind a shim.
- Copying coordinates adds ownership and allocation to a borrowed hot path.
- Returning a load-shape error for equal-length but foreign geometry gives the
  wrong diagnosis and does not encode the actual invariant.

## Acceptance

- Retain original native panic/value witnesses and exact source identity.
- Add compile-fail regressions for foreign-coordinate assembly, mutation and
  coordinate lifetime escape, plus the invalid zero-dimensional case.
- Verify exact 2-D/3-D resultant loads, additive composition, and untouched
  output on short/long load rejection for native and transparent real scalars.
- Preserve and migrate all existing analytical oracles; run full workspace,
  minimal, release, MSRV, diagnostics, package/API and documentation gates.
- Obtain independent exact-head review and use unchanged normal hooks and CI.
- Inspect current registry, consumers, workflows and deployment settings before
  any conditional merge. Existing dependency, publication and owner holds
  remain separate; this design does not authorize protected work changes.

## Design review

The independent architecture review endorses the immutable coordinate borrow
as the smallest complete repair of the one-mesh contract. Its conditions are
continued coordinate use in the mutation/lifetime rejection tests, positive
`D` enforcement in both public entry points, a `D=1` control, preservation of
the existing statics oracles, and separate exact-revision source and release
qualification. Acceptance of this design does not declare release readiness.
