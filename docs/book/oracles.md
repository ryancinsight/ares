# 10. Oracles, and what each one catches

Ares was built with **no reference implementation to compare against**. There
was no trusted solid solver in this stack to difference results against, so
every correctness claim rests on analytical oracles: closed-form solutions and
conservation identities.

That is a real risk, and the mitigation is breadth rather than depth. Each
oracle below is blind to something another catches. This chapter says what each
one is for and, as importantly, what it cannot see.

## The patch test

**The headline.** Impose a constant-strain displacement field on an arbitrary
distorted patch of elements. The discretisation must reproduce it: every cell
recovers the same constant strain, and every *interior* node carries zero
internal force.

It earns its status because it fails on almost any assembly defect — a wrong
gather or scatter index, a mapping error, a quadrature weight, a sign, a
transposed gradient — while a convergence study can pass with several of those
present.

The mechanism is exact cancellation. For an interior node `p`, the shape
function `N_p` vanishes on the patch boundary, so

```text
f_p = sum over cells of |Om_e| sigma . grad N_p = sigma . integral(grad N_p) = 0
```

by the divergence theorem, and the sum is over terms that are individually
large. Anything that perturbs one term without perturbing its partners leaves a
residual far above rounding.

The patch is deliberately **distorted** — the interior node sits off-centre, and
the 3-D corners are perturbed off the unit cube. A symmetric patch lets sign
errors cancel in mirror pairs, which is exactly the failure the test exists to
catch.

Measured, by mutation: reversing the scatter order fails eight tests; a Voigt
engineering-shear factor of two on the off-diagonal strain fails six; dropping
the cell measure from the nodal force fails six.

### "Exact to machine precision" is a derived bound

The cancellation is exact in exact arithmetic. In floating point the nodal
values `u = a + G x` already carry rounding, so the assertion is a bound derived
from the problem rather than an equality:

```text
|f_p| <= n |Om| (lambda + 2 mu) max|u| max|grad N|^2 eps
```

Every factor is measured from the patch itself — the measures, the shape
gradients, the displacement magnitude — rather than fitted.

### What the patch test cannot see

It cannot distinguish a formulation that is translation-invariant by
construction from one that merely rounds small. Both leave a residual of the
same order, because `u = a + G x` is rounded at construction with absolute
error `|a| eps` and no later subtraction recovers what the stored value never
held.

That property has its own test — a rigid translation of the whole mesh, where
every nodal value is the *same* float, the differences are identically zero,
and the assertion is exact equality. That test, and only that test, fails when
the reference-node differencing of chapter 5 is removed.

## Hand-computed element stiffness

Three columns of the unit triangle's stiffness and one of the unit
tetrahedron's, derived by hand through the textbook Voigt route
`K = A B^T D B`.

The route is the point. The implementation goes through the full stress tensor
and never forms a `B` matrix, so the two paths share no code and their
agreement is evidence rather than restatement. The engineering shear
`gamma_xy = 2 eps_xy` appears in the reference and not in the implementation, so
a version that confused the two would disagree here while remaining internally
consistent on both sides.

## The manufactured solution

Every other oracle depends on a special geometry — a slender beam, an
axisymmetric annulus — and so tests the solver only where that geometry's
assumptions hold. A manufactured solution has no such dependence: pick any
smooth field, differentiate it through the governing equation to get the body
force that produces it, and check the solver recovers it.

On the unit square with `f = sin(pi x) sin(pi y)`, the field `u = (A f, B f)`
vanishes on the whole boundary, so the Dirichlet data is homogeneous and no
boundary term contaminates the interior error.

### The convergence study, and a claim it corrected

Refinement should recover second-order convergence in the `L2` displacement
norm. The measured rates are **1.652, 1.884, 1.967**.

The first is below two, and the original assertion — a floor of 1.8 at every
refinement step — failed on it. The rates are monotonically approaching two
from below, which is the signature of a second-order method observed outside
the asymptotic regime, not of a first-order one: at `h = 0.35` a single sine
hump spans four elements.

The assertion now checks the *shape* of the approach — rates rising, finest
above 1.9 — which uses every data point and is harder to satisfy by accident. A
genuinely first-order method gives rates flat near one and fails both halves at
any starting mesh.

### A second claim the same study corrected

The body-force documentation asserted that a *lumped* load would cap refinement
at first order. Mutation testing falsified it: the vertex rule is also exact for
a linear integrand, so it is second order too, and substituting it leaves the
rate study unchanged.

Measured on the same problem, the lumped load is in fact the **more** accurate
of the two — `3.8e-3` relative against `1.2e-2` — because its error partly
cancels the discretisation's rather than adding to it.

The consistent load is kept for being the Galerkin integral rather than an
approximation of it, and because work conservation across a coupling interface
is a property of the consistent form. Not for an order it does not buy.

## Lamé's thick-walled cylinder

A cylinder under internal pressure, against the closed form

```text
u_r(r) = (1 + nu) a^2 p / (E (b^2 - a^2)) [ (1 - 2 nu) r + b^2 / r ]
```

It is the only oracle here with a **curved boundary**. Every other fixture is a
rectangle whose edges the mesh represents exactly, so none of them can detect a
defect that appears only when element edges approximate a surface — a traction
resolved along the wrong normal, or a facet measure right for an axis-aligned
edge and wrong for an oblique one.

The headline assertion is a convergence *rate* rather than a tolerance. The
straight-edged mesh approximates the circle, so geometry carries an error of its
own on top of the discretisation's; a fixed tolerance on one mesh would be a
fitted number standing in for two effects at once. Both are second order, so
their sum must be.

The direction is asserted separately: the wall must move outward, and the inner
wall further than the outer. An inward traction would still converge, still be
smooth, and still refine cleanly — differing from the closed form only by a
sign, which a relative error norm reports as roughly two rather than as the
reversal it is.

## Cantilever tip deflection

`delta = P L^3 / (3 E I)` — Euler-Bernoulli beam theory, an **independent
structural model** with its own assumptions. Every other oracle here is
elasticity checking elasticity.

Approached from the known direction, per chapter 5: linear triangles are stiff
in bending, so the computed deflection sits *below* the beam value and rises
toward it under refinement. A result that matched on a coarse mesh, or
overshot, would be evidence of a defect — most likely an element that has lost
stiffness somewhere.

The assertions are therefore structural rather than a fitted threshold: below
the beam value, monotonically increasing under refinement, and the gap closing.
None needs a tuned constant, and each fails for a different defect.

The exact two-dimensional answer slightly *exceeds* Euler-Bernoulli, by a shear
term of order `(H/L)^2`, so "below" is a statement about element stiffness
dominating at these resolutions rather than a universal bound.

## Energy consistency

Strain energy equals external work. A conservation identity, so it holds
**exactly** on any mesh at any resolution — unlike every closed-form comparison
here, which holds only in the limit.

That makes it the one oracle separating a solve that is merely *inaccurate*
from one that is *inconsistent*. A coarse mesh gives a poor cantilever
deflection while satisfying this identity to the solver tolerance; a defect in
assembly or in the constrained operator breaks it at any resolution.

## Scalar generality

The accuracy and identity oracles run at `f32` as well as `f64`, which is what
forces the whole solve path — assembly, constrained operator, Athena adapter,
Krylov solver — to monomorphise at a second scalar at all. A body pinned to a
concrete type does not compile against a second one.

The convergence-*rate* studies stay at `f64`, and the reason is stated rather
than the studies quietly omitted: at `f32` the representable relative precision
is about `1e-7`, and the discretisation error passes below that by the third
refinement. A rate measured there describes rounding. It would look like
evidence without being any.

## On mutation testing

Several claims above are stated as measured rather than argued. That is
deliberate: an oracle nobody has tried to break is a claim about a test, not
about the code.

Every headline oracle in this crate has been checked by injecting the defect it
claims to catch and confirming it fails — and, in three cases, by discovering
that it does *not*, which is how the false claims recorded above were found.
