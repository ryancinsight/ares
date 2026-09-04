# 7. Boundary conditions

Chapter 3 established that the balance law alone does not determine a
displacement. This chapter is how the two kinds of condition are imposed, and
why both are typed rather than expressed as index arithmetic on an assembled
system.

Both failure modes here are silent. A Dirichlet condition applied by striking
rows out of a matrix, and a Neumann load applied as a force where a traction
was meant, each produce a system that still solves.

## Neumann: prescribed traction

A traction is a **stress** — force per unit area — not a force. The two differ
by the area of the facet they act on, and on a small test mesh the difference
looks like a modelling choice rather than an error.

```rust
# extern crate ares;
use ares::{TractionBoundary, TractionFacet};

// A unit square; pull its right-hand edge outward at 1 MPa.
let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let facets = [TractionFacet::new([1, 2], [1.0e6_f64, 0.0])];
let boundary = TractionBoundary::try_new(&facets, &nodes).expect("a valid facet");

let mut loads = [0.0_f64; 8];
boundary.add_consistent_loads(&nodes, &mut loads).expect("well-shaped");

// The resultant is traction times edge length, split between the two nodes.
let total: f64 = (0..4).map(|n| loads[n * 2]).sum();
assert!((total - 1.0e6).abs() < 1.0);
```

A facet of a `D`-simplex is a `(D - 1)`-simplex with exactly `D` nodes — an edge
in 2-D, a triangle in 3-D — so both arrays are `D` long and the shape needs no
further parameter.

### The consistent load

From chapter 4, the load is `f_a = integral of N_a t dS`, not "the traction
assigned to a node". For linear shape functions on a `(D-1)`-simplex,
`integral of N_a dS = measure / D`, so a uniform traction distributes equally:

```text
f_a = t * measure / D
```

Equal distribution is a *property* of linear elements under uniform traction,
derived from that integral — not a simplification. A quadratic element's
consistent load is famously unequal, and negative at the corners. Writing it as
the integral rather than as a split is what makes the derivation portable to
elements where the answer differs.

It is also what makes interface work exactly conserved in
[chapter 9](coupling.md).

### The facet measure

The `(D-1)`-measure comes from the Gram determinant
`sqrt(det(E^T E)) / (D-1)!`, where `E` holds the facet's edge vectors. That is
the general form of the cross-product magnitude giving a triangle's area and
the difference giving an edge's length, so one expression covers both
dimensions rather than a match on `D`.

## Dirichlet: prescribed displacement

```rust
# extern crate ares;
use ares::{DirichletConditions, PrescribedDisplacement};

// Hold node 0 in both directions; stretch node 1 along x.
let prescribed = [
    PrescribedDisplacement::new(0, 0, 0.0_f64),
    PrescribedDisplacement::new(0, 1, 0.0),
    PrescribedDisplacement::new(1, 0, 2.5e-4),
];
let conditions = DirichletConditions::<f64, 2>::try_new(&prescribed, 4)
    .expect("valid against a four-node mesh");
assert_eq!(conditions.len(), 3);
```

Conditions must be **strictly increasing** in degree of freedom. That makes
duplicate detection a single pass, and it removes a silent dependence on input
order: two conditions on one degree of freedom would otherwise resolve by
whichever the caller happened to list last. Out-of-order input is rejected
rather than sorted, because sorting would need to allocate and the caller
already knows the order it built them in.

### How the constraint is imposed

The textbook approach strikes the constrained rows and columns out of the
matrix, producing a smaller system with a different numbering. Then every
vector crossing the boundary needs mapping in both directions, and two mappings
that must agree is one more chance to disagree.

Ares instead leaves **identity rows** in place. With `P` the projection that
zeroes constrained entries, the operator is

```text
A = P K P + (I - P)
```

The constrained part of the output is the constrained part of the input. `A` is
still symmetric, and still positive definite once the conditions remove the
rigid-body modes, so conjugate gradients applies unchanged. The solution vector
is indexed identically before and after the solve.

Symmetry needs the projection on **both** sides. Applying it only to the output
would leave the coupling column `K_cf` intact, giving a non-symmetric operator
on which conjugate gradients has no convergence guarantee — and which still
returns plausible numbers.

### The term that is easy to drop

The right-hand side is

```text
b = P (f_ext - K u_g) + (I - P) g
```

where `u_g` holds the prescribed values. The `K u_g` term is there because a
non-zero prescribed displacement does work on the free degrees of freedom
through the coupling block: moving one node strains the cells that touch it,
and that strain is a load on their other nodes.

Omitting it silently solves a *different* problem — the one where every
prescribed displacement is zero. That problem has a perfectly convergent
solution, so nothing reports the substitution.

And with every prescribed value zero the term vanishes, so the omission is
invisible in the fixed-at-zero case, which is the common one. That is why the
crate's tests prescribe a non-zero value: it is the only configuration in which
the defect can be detected at all.
