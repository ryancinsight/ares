# Ares

Solid momentum balance for the [Atlas](https://github.com/ryancinsight/atlas)
stack. Ares owns the balance side of solid mechanics — kinematics, stress,
equilibrium, and the boundary conditions that close them — and no material
data at all: constitutive closure belongs to
[Proteus](https://github.com/ryancinsight/proteus). Proteus closes, Ares
balances (atlas ADR 0055).

## Crates

| Crate | Registry name | What it is |
| --- | --- | --- |
| [`crates/ares`](crates/ares) | `ares-solid` | The domain core. `#![no_std]`, allocation-free, depends only on vocabulary crates. |
| [`crates/ares-athena`](crates/ares-athena) | `ares-athena` | The Athena linear-operator seam. Links `std` through `leto`. |

The split keeps the domain core free of infrastructure — see
[ADR 0001](docs/adr/0001-athena-seam-as-a-separate-crate.md). Dependencies run
strictly inward: `ares-athena` depends on `ares`, never the reverse.

## Scope

Phase 0 is small-strain linear elastostatics on an unstructured mesh, charted
by atlas ADR 0057. Plasticity, finite deformation, contact, dynamics,
fracture, fatigue, and anisotropy are later phases, and none is scaffolded — a
module for a capability that does not exist is a placeholder.

## Governing decisions

[`docs/adr`](docs/adr) records decisions local to this repository. Decisions
that bind Ares from outside it live in the
[Atlas](https://github.com/ryancinsight/atlas) meta-repository:

| Record | Subject |
| --- | --- |
| atlas ADR 0055 | Continuum domain decomposition — Proteus closes, Ares balances |
| atlas ADR 0056 | New-construction promotion path |
| atlas ADR 0057 | Ares phase 0 charter |
| atlas ADR 0059 | Fluid–structure coupling phase 0 |

## The book

[`docs/book`](docs/book) teaches continuum solid mechanics and this crate
together, from measuring deformation through to the coupling interface. It
assumes Rust and some calculus and nothing about the finite element method.
Chapter 10 is the verification chapter, and it is the one to read if you read
only one: it says what each oracle catches, what it is blind to, and where an
oracle falsified a claim its own author had written down.

## Verification

Analytical oracles throughout; there is no reference implementation to
difference against, which is why the oracle breadth is the safety net. The
headline is the **patch test**: a constant-strain field on an arbitrary
distorted patch is reproduced to machine precision, in 2-D and 3-D, under pure
shear and pure dilation, at `f32` and `f64`. Element stiffness columns are
additionally checked against hand computation through the Voigt `B^T D B`
route, which shares no code with the tensor formulation the implementation
uses.

Oracles are mutation-measured rather than assumed: each records which defects
it was observed to catch, and where a mutation revealed an overclaim the claim
was corrected rather than the mutation dismissed.

## Licence

MIT OR Apache-2.0.
