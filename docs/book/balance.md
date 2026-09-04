# 3. Balance: what equilibrium actually asserts

Chapters 1 and 2 gave a way to measure deformation and a rule turning it into
stress. The balance law is what ties them to the loads you actually apply.

## The statement

Take any sub-volume of a body at rest. The forces on it must sum to zero: the
tractions its neighbours exert across its surface, plus any body force such as
gravity acting throughout it. Apply the divergence theorem to convert the
surface term into a volume one, and since the sub-volume was arbitrary, the
integrand must vanish everywhere:

```text
div sigma + b = 0
```

That is static equilibrium. `b` is the body force per unit volume, and
`div sigma` is the vector whose `i`-th component is `d sigma_ij / d x_j`.

Everything else in this book is machinery for solving that equation on a real
geometry.

## Why it needs boundary conditions

The equation on its own does not determine the displacement. Two reasons, and
they are different.

First, it constrains stress, and stress comes from strain, and strain comes
from the *gradient* of displacement. Any rigid motion added to a solution is
still a solution — the body would be in equilibrium while floating away. So
the problem has a null space that must be removed by holding something.

Second, the equation says nothing about what happens at the surface. A bar
pulled at one end and a bar with that end free are the same differential
equation and different problems.

So the balance law is only half of a well-posed problem. The other half is the
boundary, and it comes in exactly two kinds:

- **Dirichlet**, prescribing the displacement — a clamp, a support, a
  symmetry plane;
- **Neumann**, prescribing the traction — a pressure, an applied load, or a
  free surface, which is the case `t = 0`.

Every point of the boundary needs one or the other, in each direction, and no
point may have both in the same direction. [Chapter 7](boundary.md) covers how
each is imposed.

## The complete problem

Putting it together, Phase 0 solves:

```text
div sigma + b = 0          in the body
sigma = lambda tr(eps) I + 2 mu eps
eps   = (1/2)(grad u + grad u^T)
u = g                      on the Dirichlet boundary
sigma . n = t              on the Neumann boundary
```

Three of those five lines are owned elsewhere or established earlier: the
constitutive line is Proteus's, the strain line is chapter 1. What Ares adds
is the first line and the last two.

## Plane strain, and why 2-D is a real problem rather than a toy

A two-dimensional analysis is not a cartoon of a three-dimensional one; it is a
three-dimensional problem with an assumption attached, and there are two
different assumptions in common use.

**Plane strain** assumes the body is long in the third direction and prevented
from stretching along it. A dam, a tunnel lining, a long pipe. The out-of-plane
strain is zero — but the out-of-plane *stress* is not, because holding a
material that wants to contract requires force.

**Plane stress** assumes the opposite: a thin plate, free to thin, so the
out-of-plane stress is zero and the strain is not.

Ares's 2-D case is **plane strain**, and it is so because of what
`isotropic_hooke` does rather than because of a flag. Applying the
three-dimensional Hooke's law with the out-of-plane strain set to zero gives
exactly the plane-strain relation. Nothing was configured; the assumption is
the consequence of writing the general law in two dimensions.

This is worth stating plainly because it changes the answers. A cantilever
computed under plane strain is stiffer than under plane stress by a factor of
`1 / (1 - nu^2)` — about 10% for a typical metal. The verification chapter's
beam comparison uses the plane-strain modulus for that reason, and using the
other one would have shown a 10% discrepancy that no amount of mesh refinement
would remove.

## What "static" excludes

There is no time in these equations, and no mass. The body is not accelerating;
it has reached equilibrium and stopped. This excludes vibration, wave
propagation, impact, and anything where inertia matters.

That is a real restriction and it is the boundary of Phase 0. It also shapes
the coupling in [chapter 9](coupling.md): with no time in the structural
problem there is no velocity continuity to impose at a fluid interface, which
is one of the conditions a full fluid-structure formulation would need.
