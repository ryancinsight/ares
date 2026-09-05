# ares — Solid Momentum Balance for Atlas

`ares` computes how a solid deforms under load. Given a mesh, a material, and
the forces and restraints acting on a body, it answers: where does every point
of that body end up, and what stress does it carry there?

This book teaches the subject and the crate together. It assumes you can
program in Rust and remember some calculus, and assumes nothing about
continuum mechanics or the finite element method. Part I builds the physics
from the idea of measuring deformation. Part II turns the resulting equation
into something a computer can solve. Parts III and IV cover the solve, the
coupling to other physics, and — at length, because it is the part that
decides whether any of it is true — how the result is verified.

## What Ares owns, and what it deliberately does not

Ares owns the **balance** side of solid mechanics: kinematics, stress
measures, equilibrium, and the boundary conditions that close them.

It owns **no material data at all**. There is no steel in this crate, no
aluminium, no table of moduli. Material response belongs to
[Proteus](https://github.com/ryancinsight/proteus), and the division is
deliberate rather than tidy-minded: a material property is a claim about the
world that needs a source and a validity range, while a balance law is a
mathematical statement that is true regardless of what the body is made of.
Mixing them means every new alloy touches the equilibrium code.

Atlas ADR 0055 states the split in four words: **Proteus closes, Ares
balances.**

The same rule places the other neighbours. Gaia owns meshes and geometry.
Athena owns solver policy. Harmonia owns every coupling to another physics
domain, which is why Ares has no dependency on a fluid or acoustic package in
any phase — see [chapter 11](stack_position.md).

## Scope of this phase

Phase 0 is **small-strain linear elastostatics on an unstructured mesh**:

- displacements small enough that the geometry does not change appreciably;
- a linear relationship between strain and stress;
- no time dependence — the body is in equilibrium, not in motion.

Plasticity, finite deformation, contact, dynamics, fracture, fatigue, and
anisotropy are later phases. None of them is scaffolded here, and that is a
policy rather than an oversight: an empty module for a capability that does
not exist is a placeholder, and placeholders are how a codebase comes to
promise things it cannot do.

## How to read the verification chapter

If you read only one chapter, read [chapter 10](oracles.md). Ares was built
with no reference implementation to compare against, so every claim it makes
rests on analytical oracles — closed-form solutions and conservation
identities. That chapter explains what each one is blind to, which is the only
honest way to describe a test suite, and records the cases where an oracle
falsified a claim its own author had written down.

## The crates

| Crate | What it is |
| --- | --- |
| `ares-solid` (imported as `ares`) | The domain core. `no_std`, allocation-free, depends only on vocabulary crates. |
| `ares-operator` | The Athena linear-operator seam. |
| `ares-coupling` | The Harmonia coupling partition. |

The split exists because Athena's solver seam can only be implemented against
a named backend, and the only host backend links `std`. Keeping that out of
the domain core is what lets the balance code stay `no_std` and
allocation-free. [ADR 0001](https://github.com/ryancinsight/ares/blob/main/docs/adr/0001-athena-seam-as-a-separate-crate.md)
records the decision.
