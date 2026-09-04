# 5. The linear simplex element

The weak form needs basis functions. Ares uses the simplest ones that work:
linear functions on triangles in 2-D and tetrahedra in 3-D.

## What a simplex is, and why it is the natural cell

A **simplex** in `D` dimensions is the shape with the fewest vertices that
encloses any volume at all: a triangle in 2-D, a tetrahedron in 3-D. `D + 1`
nodes.

That count is exactly what a linear function needs. A linear function in `D`
dimensions has `D + 1` coefficients — one constant plus one slope per axis — so
`D + 1` nodal values determine it uniquely. Simplices and linear interpolation
fit each other exactly, with nothing left over and nothing underdetermined.

The count appears in the type:

```rust
# extern crate ares;
use ares::Simplex;

// A triangle. `D = 2` and `N = 3` are both inferred from the argument.
let triangle = Simplex::new(&[[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]]);
assert_eq!(triangle.signed_measure(), 0.5);
```

`N` is `D + 1` and stable Rust cannot spell that as an array length, so it is a
second const parameter. That is not a workaround: it makes the node count a
*compile-time* fact, so construction cannot fail and the displacement and force
buffers cannot be misshaped. A mismatched pair does not compile.

## Shape functions

The basis function `N_a` is one at node `a`, zero at the others, and linear
between. On a simplex these are the **barycentric coordinates**.

Two properties do the work.

They **sum to one** at every point, so their gradients sum to zero. This is
what makes a constant displacement reproduce exactly — it is the algebraic form
of "translating the body strains nothing".

Their **gradients are constant** over the element, since the functions are
linear. Every quantity built from them is therefore constant too. Strain is
constant per cell, stress is constant per cell, and hence the name for this
family: **constant strain triangles**.

## No quadrature loop, and why that is not a shortcut

A general finite element integrates over each cell with a numerical quadrature
rule — evaluate at a few sample points, weight, sum.

Here every integrand is constant, so `integral = value × measure` is **exact**.
Not accurate enough: exact, with no error term to bound.

Ares therefore has no quadrature loop, and adding one would be adding a rule
whose error is identically zero. Higher-order elements will need one; these do
not.

## The measure carries a sign

```rust
# extern crate ares;
use ares::Simplex;

let forward = Simplex::new(&[[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]]);
let reversed = Simplex::new(&[[0.0_f64, 0.0], [0.0, 1.0], [1.0, 0.0]]);
assert_eq!(forward.signed_measure(), 0.5);
assert_eq!(reversed.signed_measure(), -0.5);
```

`signed_measure` returns a negative value for a cell whose nodes are wound the
wrong way, rather than taking an absolute value. A negative measure is not a
small cell, it is an **inverted** one, and inversion negates that cell's
stiffness contribution. A mesh with mixed winding assembles an operator that
stores negative energy under deformation — indefinite, and conjugate gradients
has no reason to converge on it.

[Chapter 6](assembly.md) rejects such a mesh at construction for that reason.
Silently taking the absolute value would hide an inverted mesh and produce a
plausible wrong answer.

## The stiffness action

Rather than build an element matrix, Ares computes its **action** on a
displacement:

```text
eps   = sum_a sym((u_a - u_0) (x) grad N_a)
sigma = C : eps
f_a   = measure * sigma . grad N_a
```

Three lines, each the responsibility of a different chapter: kinematics, then
the Proteus closure, then balance. It reuses the strain measure and the
constitutive law unchanged, so a defect in either surfaces here rather than
being duplicated in a separate path — and the separate path, the `B` matrix, is
exactly where the Voigt shear factor of chapter 1 is normally introduced.

## The subtraction in the first line

`(u_a - u_0)` rather than `u_a` is the one piece of this that looks like a typo
and is not.

The two are mathematically identical, because the shape gradients sum to zero,
so subtracting a constant from every nodal displacement changes nothing.
Numerically they differ, and the difference was found by a test failing.

The gradients cancel exactly only when summed in the order they were built.
Downstream they are re-accumulated in a different order, so their sum is zero
to *rounding* rather than identically — measured at `1.4e-17` for an ordinary
triangle. Under the plain form a rigid translation therefore leaves a residual
gradient proportional to the translation, which becomes a spurious stress that
grows with how far the body has moved and not with how fine the mesh is. No
refinement study would ever reveal it.

Differencing against node 0 removes it at the source. A uniform translation
makes every `u_a - u_0` exactly zero — the same float minus itself — so the
gradient is exactly zero whatever the shape gradients rounded to. Translation
invariance stops being a cancellation that happens to work out and becomes a
property of the formulation.

Rotation is a different matter and is **not** exact: its relative displacements
do not vanish, so the reconstructed gradient is antisymmetric only to rounding.
The crate's tests bound it rather than asserting equality, and an earlier
version that asserted exactness was passing only because the geometry it
happened to pick — the unit triangle, whose edge matrix is the identity —
carries no rounding at all.

## The cost of choosing the simplest element

Linear simplices are stiff in bending. Their strain is constant per cell, and
bending needs strain that varies linearly through the depth, which one element
cannot represent. So a bending problem converges from below: the computed
deflection is too small and rises toward the true value under refinement.

That is a real limitation with a real consequence — the cantilever oracle in
[chapter 10](oracles.md) asserts the *direction* of the error rather than its
absence, because a result that matched beam theory on a coarse mesh would be
evidence of a defect rather than of accuracy.
