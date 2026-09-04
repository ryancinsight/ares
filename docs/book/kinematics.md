# 1. Kinematics: measuring deformation

Before you can ask what force a deformed body carries, you need a number that
says *how deformed it is*. Getting that number right is most of the difficulty,
because the obvious candidates are all wrong in the same way: they report
deformation when nothing has deformed.

## The displacement field

A body occupies some region of space. Load it, and every point `x` moves to a
new position. The **displacement** `u(x)` is the vector from where a point was
to where it went.

Displacement alone is not deformation. Pick up a steel bar and carry it across
the room: every point has moved by metres, and the bar is not strained at all.
So whatever measures deformation cannot be `u` itself.

## What actually deforms is the *gradient*

The bar is unstrained because every point moved by the *same* amount. A body
deforms when nearby points move by *different* amounts — when the displacement
varies across the body. That variation is the **displacement gradient**:

```text
(grad u)_ij = du_i / dx_j
```

a `D × D` matrix in `D` dimensions. For a rigid translation it is zero
everywhere, which is the first thing we wanted.

But it is not enough on its own. Rotate the bar instead of translating it, and
the gradient is *not* zero — points on opposite sides of the bar move in
opposite directions. Yet the bar is still unstrained.

## Splitting rotation out

Any matrix splits uniquely into a symmetric and an antisymmetric part:

```text
grad u = sym(grad u) + skew(grad u)
```

For a small rigid rotation the gradient is *entirely* antisymmetric. So the
symmetric part is exactly what survives when rotation is removed. That is the
**small-strain tensor**:

```text
eps = (1/2) (grad u + grad u^T)
```

It is zero for translation, zero for infinitesimal rotation, and non-zero
precisely when the body's shape or volume has changed. In the crate:

```rust
# extern crate ares;
use ares::SmallStrain;

// A pure stretch along x: every point moves proportionally to its own x.
let gradient = [[1.0e-3, 0.0], [0.0, 0.0]];
let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&gradient);
assert_eq!(strain.volumetric(), 1.0e-3);
```

The word **small** is load-bearing. This measure is exact only in the limit of
infinitesimal displacement; a finite rotation still produces a small spurious
strain, of order the rotation angle squared. Phase 0 assumes displacements
small enough for that to be negligible, and the assumption is the boundary
between this phase and the finite-deformation one.

## Reading the tensor

Two decompositions of the strain matter enough to have their own methods.

**Volumetric strain** is the trace, `tr(eps)`, and to first order it is the
relative change in volume. A material that resists volume change resists this
part.

**Deviatoric strain** is what is left after removing the volumetric part:
`dev(eps) = eps - tr(eps) I / D`. It is shape change at constant volume —
shear. Metals yield in response to this part and are almost indifferent to the
volumetric one, which is why the split is not merely algebraic.

```rust
# extern crate ares;
use ares::SmallStrain;

// Pure shear: no volume change, all shape change.
let gradient = [[0.0, 2.0e-3], [2.0e-3, 0.0]];
let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&gradient);
assert_eq!(strain.volumetric(), 0.0);
```

## A note on shear that has cost people real money

There are two conventions for shear strain and they differ by a factor of two.
The **tensor** shear is `eps_xy`, the actual off-diagonal component. The
**engineering** shear is `gamma_xy = 2 eps_xy`, which is what appears in the
Voigt vector notation most finite element textbooks use.

Both are correct; mixing them is a factor-of-two error in every shear term,
and it is a bad one because the result stays smooth, stays symmetric, and
converges under refinement — to the wrong answer.

Ares avoids the trap structurally rather than by care: it carries strain and
stress as full tensors and never forms a Voigt vector at all, so there is no
place for the convention to be chosen wrongly. The hand-computed stiffness
check in [chapter 10](oracles.md) deliberately derives its reference through
the Voigt route *because* the two paths would disagree if the factor were
mishandled.

## Why the null space is asserted exactly

The property that rigid motion produces zero strain is not a nicety, and the
crate asserts it as an exact equality rather than a tolerance wherever it can.

The reason is the shape of the failure. A strain measure that returned a tiny
non-zero value under rigid motion would manufacture stress out of nothing, and
that stress would scale with **how far the body had moved**, not with how fine
the mesh was. Every convergence study would pass. The error would only appear
in a problem where something moved a long way — which, in a coupled simulation,
is exactly the interesting case.

[Chapter 5](element.md) returns to this, because making the property survive
the discretisation took more than making it true in the continuum.
