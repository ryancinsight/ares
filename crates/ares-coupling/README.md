# ares-coupling

The [Harmonia](https://github.com/ryancinsight/harmonia) coupling partition for
the [Ares](https://github.com/ryancinsight/ares) solid momentum balance.

Atlas ADR 0059 splits Phase 0 fluid-structure interaction so `CFDrs` produces
the interface traction from its own flow state and a structural solver consumes
it. Neither depends on the other; this is the structural half, joined to the
fluid half by a Harmonia driver.

## Installation

```toml
[dependencies]
ares-coupling = "0.3.0"
```

## Usage

The example couples a unit square of two triangles, held along `x = 0`, with
the `x = 1` edge as the interface. The interface lists its nodes as `[2, 1]`,
so the exported displacement follows that list and not mesh order. It also uses
`ares` (`ares-solid`), `ares-operator`, `athena-core`, `harmonia`, `aequitas`
and `proteus-mat`.

```rust
use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{DirichletConditions, PrescribedDisplacement, SimplexMesh};
use ares_coupling::{StructuralInterface, StructuralPartition};
use ares_operator::ConstrainedStiffness;
use athena_core::ConvergencePolicy;
use harmonia::Partition as _;
use proteus::IsotropicModuli;

let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let cells = [[0, 1, 2], [0, 2, 3]];
let clamped = [
    PrescribedDisplacement::new(0, 0, 0.0),
    PrescribedDisplacement::new(0, 1, 0.0),
    PrescribedDisplacement::new(3, 0, 0.0),
    PrescribedDisplacement::new(3, 1, 0.0),
];
let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
let conditions = DirichletConditions::try_new(&clamped, mesh.node_count()).expect("valid");
let moduli = IsotropicModuli::from_young_poisson(
    Pressure::from_base(200e9),
    Dimensionless::from_base(0.3),
)
.expect("inside the positive-definite domain");
let operator = ConstrainedStiffness::new(mesh, moduli, conditions);

let interface = StructuralInterface::try_new(&[2, 1], &[[1, 2]], &mesh).expect("conforming");
let policy = ConvergencePolicy::<f64>::new(1e-18, 1e-12, 200).expect("valid policy");
let mut partition =
    StructuralPartition::try_new(mesh, operator, interface, policy).expect("workspace");

// One facet, so the traction exchange is `[t_x, t_y]`: 1 MPa in +x.
let mut state = vec![0.0; partition.state_dimension()];
partition.solve_for_traction(&mut state, &[1.0e6, 0.0]).expect("the solve converges");

// The displacement exchange is node major, in interface order `[2, 1]`.
let mut exported = vec![0.0; partition.output_dimension()];
partition.export(&state, &mut exported).expect("well-shaped");
assert_eq!(exported[..2], state[4..6]);
assert_eq!(exported[2..], state[2..4]);
assert!(exported[0] > 0.0 && exported[2] > 0.0);
```

## What it owns

Nothing but the adapter. Ares balances, Athena solves, Harmonia drives the
partitions. This crate carries no dependency on `CFDrs`, `kwavers`, or any other
balance domain, and cannot: atlas ADR 0055 R7 forbids a direct
balance-to-balance edge and routes coupling through Harmonia.

## Phase 0 is one-way

Traction flows fluid to solid. Displacement is exported but does not move the
fluid mesh, because `CFDrs` has no mechanism to move it (atlas ADR 0059). The
partition is quasi-static: `advance` solves equilibrium for the traction it was
handed and does not integrate in time, so the substep is unused.

## The exchange ordering contract

- **input** — traction, facet index major, component minor: `input[f * D + c]`.
- **output** — displacement, node index major, component minor.

The orderings differ because the quantities live on different entities. A
traction is a stress resolved on a surface, and the fluid side computes one per
face; a displacement is a property of a point. Forcing either onto the other's
entity would mean interpolating, which ADR 0050 places outside Harmonia.
ADR 0059 states the contract as "interface node index major, component minor",
which describes the displacement side and not the traction side.

## Verification

The headline is **interface work conservation**: the work the traction does on
the interface equals twice the strain energy the structure stores. The two sides
are computed by routes that share no arithmetic — a facet integral against a
stiffness assembly.

With traction constant per facet and displacement linear over it, the facet
integral equals `sum u_a . f_a` over the consistent nodal load, so the identity
holds at any resolution up to the solve's residual: the test
`interface_work_equals_the_strain_energy` asserts a relative difference below
`1e-9`. A lumped load with the same resultant force would not satisfy it
(derivation in the `StructuralInterface` documentation).

## Licence

MIT OR Apache-2.0.
