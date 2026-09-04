# 6. Assembly without a matrix

Chapter 5 gave the stiffness action on one cell. Assembly walks the mesh,
gathers each cell's displacements, applies that action, and scatters the
resulting forces back to the nodes.

## The mesh view

```rust
# extern crate ares;
use ares::SimplexMesh;

// A unit square, split into two triangles.
let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let cells = [[0, 1, 2], [0, 2, 3]];
let mesh = SimplexMesh::try_new(&nodes, &cells).expect("a conforming mesh");
assert_eq!(mesh.node_count(), 4);
assert_eq!(mesh.degrees_of_freedom(), 8);
```

`SimplexMesh` borrows; it does not own. Gaia owns mesh generation, geometry,
and proximity queries, and this is the shape assembly reads rather than a
representation competing with Gaia's. Any producer that can lend node
coordinates and connectivity satisfies it.

## What construction rejects, and why it checks that way

`try_new` refuses every cell assembly could not integrate: a node index outside
the mesh, a cell with no measure, a cell wound the wrong way, coordinates that
are not finite.

It establishes this by running the *same* `shape_gradients` call assembly will
run, rather than inferring success from the measure. The two share a
determinant but not a pivot sequence, so a measure test would leave the
inference one step short of a proof. Running the identical computation on the
identical data makes the invariant an observed outcome.

That is what makes `internal_forces` fail only on a misshaped field, which in
turn is what lets the Athena operator of [chapter 8](solving.md) report the
errors its backend defines rather than needing one of its own. The validation
is upstream because a constraint downstream demanded it.

## Applying the operator

```rust
# extern crate ares;
# extern crate aequitas;
# extern crate proteus;
use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::SimplexMesh;
use proteus::IsotropicModuli;

let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let cells = [[0, 1, 2], [0, 2, 3]];
let mesh = SimplexMesh::try_new(&nodes, &cells).expect("a conforming mesh");
let moduli = IsotropicModuli::from_young_poisson(
    Pressure::from_base(200.0e9),
    Dimensionless::from_base(0.3),
).expect("admissible");

// A rigid translation: every node moves by the same amount.
let displacement = [0.5_f64, -0.25].repeat(4);
let mut forces = [0.0_f64; 8];
mesh.internal_forces(&moduli, &displacement, &mut forces)
    .expect("well-shaped fields");

// Exactly zero, not approximately: the differencing of chapter 5.
assert!(forces.iter().all(|f| *f == 0.0));
```

## Why fields are flat and coordinates are not

Node coordinates arrive as `&[[T; D]]`; displacement and force fields as flat
`&[T]`. The asymmetry is deliberate.

Geometry is built once and never leaves the crate. A field crosses the solver
boundary on *every* Krylov iteration, and Athena's vector views are flat — so a
nodal field type would impose a copy per iteration. Flat storage makes that
boundary zero-copy, and `as_chunks` recovers the nodal view inside the loop for
nothing.

Degree of freedom `node * D + component`. That ordering is the contract every
exchange in this crate keeps.

## No global matrix

The textbook next step is to assemble `K` into a sparse matrix. Ares does not,
because nothing downstream reads one.

Athena's solvers take a `LinearOperator` and ask it for `K u`. Building the
sparse structure would mean computing a sparsity pattern, converting
coordinate to compressed form, and storing the result — all so it can be
multiplied by a vector once per iteration. Assembly does that multiplication
directly from the mesh instead.

The loop is allocation-free: each cell gathers into stack buffers sized by the
compile-time node count, so it touches the heap zero times regardless of how
large the mesh is.

`internal_forces` also *writes* rather than accumulates. A caller reusing a
buffer across Krylov iterations would otherwise sum every previous iteration
into the current one, which converges to nonsense rather than failing.

## Recovering stress

Linear simplices carry constant strain, so one tensor per cell is the complete
answer rather than a sample of one:

```rust
# extern crate ares;
# extern crate aequitas;
# extern crate proteus;
# use aequitas::systems::si::quantities::{Dimensionless, Pressure};
# use ares::SimplexMesh;
# use proteus::IsotropicModuli;
# let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
# let cells = [[0, 1, 2], [0, 2, 3]];
# let mesh = SimplexMesh::try_new(&nodes, &cells).expect("mesh");
# let moduli = IsotropicModuli::from_young_poisson(
#     Pressure::from_base(200.0e9), Dimensionless::from_base(0.3)).expect("ok");
// A uniform stretch: u_x = 1e-3 * x.
let displacement: Vec<f64> = nodes.iter().flat_map(|p| [1.0e-3 * p[0], 0.0]).collect();
for strain in mesh.cell_strains(&displacement).expect("well-shaped") {
    let stress = ares::isotropic_hooke(&moduli, &strain);
    assert!(*stress.mean_stress().as_base() > 0.0);
}
```

`cell_strains` stops at the kinematic quantity deliberately. Stress follows by
passing each strain through `isotropic_hooke`, and doing that step inside would
give the crate two constitutive closures to keep agreeing.
