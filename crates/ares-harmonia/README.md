# ares-harmonia

The [Harmonia](https://github.com/ryancinsight/harmonia) coupling partition for
the [Ares](https://github.com/ryancinsight/ares) solid momentum balance.

Atlas ADR 0059 splits Phase 0 fluid-structure interaction so CFDrs produces the
interface traction from its own flow state and a structural solver consumes it.
Neither depends on the other; this is the structural half.

## What it owns

Nothing but the adapter. Ares balances, Athena solves, Harmonia drives the
partitions. This crate carries no dependency on CFDrs, kwavers, or any other
balance domain, and cannot: atlas ADR 0055 R7 forbids a direct
balance-to-balance edge and routes coupling through Harmonia.

## The exchange ordering contract

- **input** — traction, facet index major, component minor: `input[f * D + c]`.
- **output** — displacement, node index major, component minor.

The orderings differ because the quantities live on different entities. A
traction is a stress resolved on a surface, and the fluid side computes one per
face; a displacement is a property of a point. Forcing either onto the other's
entity would mean interpolating, which ADR 0050 places outside Harmonia.

## Verification

The headline is **interface work conservation**: the work the traction does on
the interface equals twice the strain energy the structure stores. It is a
conservation identity, so it holds at any resolution, and the two sides are
computed by routes that share no arithmetic — a facet integral against a
stiffness assembly.

It is exact rather than bounded, and exact *because* the nodal load is the
consistent one. Measured: replacing the consistent load with a lumped one of
the same resultant force breaks this test and no other.

## Licence

MIT OR Apache-2.0.
