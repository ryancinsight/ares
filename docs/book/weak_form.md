# 4. From the balance law to a linear system

`div sigma + b = 0` is a statement about every point of a continuous body. A
computer has finitely many numbers. This chapter is the bridge, and it is worth
following once because the shape of the finished code — matrix-free assembly,
consistent loads, why a traction is not a force — all follows from it.

## The problem with the strong form

Solving `div sigma + b = 0` directly requires the displacement to be twice
differentiable, since stress involves one derivative and the divergence another.
That is more smoothness than real solutions have: at a corner, or where the
material changes, the stress is discontinuous and the second derivative does
not exist.

The formulation is also awkward to discretise. There is no obvious way to
enforce a pointwise differential equation with a finite set of unknowns.

## Multiplying by a test function

The standard move: instead of requiring the equation to hold at every point,
require its *weighted average* to vanish for every weight in a suitable family.

Take any **test function** `v` — think of it as a virtual displacement, an
imaginary wiggle you could apply to the body — that is zero wherever the
displacement is prescribed. Multiply and integrate:

```text
integral over Omega of v . (div sigma + b) dV = 0
```

If this holds for *every* admissible `v`, the original equation holds too. So
far nothing is gained.

The gain comes from integrating by parts, which moves one derivative off the
stress and onto the test function:

```text
integral of grad v : sigma dV = integral of v . b dV + integral over the
                                boundary of v . (sigma . n) dS
```

This is the **weak form**, and it is better in three distinct ways.

It needs only *one* derivative of the displacement, so the solutions it admits
are the ones real problems have.

The boundary term is `v . t` — the test function against the **traction**. So a
Neumann condition is not something imposed afterwards; it appears in the
equation naturally, as the boundary term. That is why chapter 7 can say a
traction condition is applied by adding to the right-hand side and nothing
else.

And since `sigma` is symmetric, `grad v : sigma` equals `eps(v) : sigma`. The
left side becomes strain against stress, which is the **virtual work** of the
internal forces. The whole equation reads: internal virtual work equals
external virtual work, for every virtual displacement.

## Discretising

Now make it finite. Chop the body into cells and pick a finite set of basis
functions `N_a`, one per node, each equal to one at its own node and zero at
every other. Approximate the displacement as a weighted sum:

```text
u(x) = sum over nodes a of N_a(x) u_a
```

The unknowns are now the nodal values `u_a` — a finite list of numbers.

**Galerkin's choice** is to use the same family for the test functions. Take
`v = N_b` for each node `b` in turn, and the single weak-form equation becomes
one equation per node:

```text
sum over a of K_ab u_a = f_b
```

with

```text
K_ab = integral of eps(N_b) : C : eps(N_a) dV        the stiffness
f_b  = integral of N_b . b dV + integral of N_b . t dS   the load
```

A linear system. That is the whole derivation.

## Three consequences worth carrying forward

**The load is an integral, not an assignment.** `f_b` is `integral of N_b t dS`,
not "the traction at node b". Chapter 7 calls this the *consistent* load, and
the distinction is invisible on a uniform mesh under uniform load and matters
everywhere else. It is also the reason interface work is exactly conserved in
[chapter 9](coupling.md).

**`K` is symmetric.** `K_ab = K_ba` follows from the symmetry of `C`. That is
what permits conjugate gradients in [chapter 8](solving.md), and it is why the
constrained operator there is built carefully enough to preserve it.

**`K` is never actually formed.** The system needs `K u`, not `K`. Chapter 6
computes that product directly from the mesh, and the sparse matrix that a
textbook would assemble is a structure nothing in this crate reads.
