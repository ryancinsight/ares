# 11. Position in the stack

Ares is one member of the [Atlas](https://github.com/ryancinsight/atlas) stack.
This chapter says what it depends on, what depends on it, and — the part that
took an ADR to settle — where the boundaries with its neighbours fall.

## The decomposition axis

Atlas ADR 0055 decomposes continuum physics along one axis. Every domain is a
**conserved quantity**, a **balance operator** for it, and a **constitutive
closure**:

| Layer | Owner | Conserved quantity |
| --- | --- | --- |
| Balance | `CFDrs` | fluid momentum and mass |
| Balance | `kwavers` | acoustic momentum |
| Balance | `helios` / `hyperion` | radiative energy |
| Balance | `asclepius` | bio-heat |
| Balance | **`ares`** | **solid momentum** |
| Balance | `prometheus` | species mass (chartered, not yet built) |
| Closure | `proteus` | material response — closes all of them |
| Coupling | `harmonia` | routes between balance domains |

The axis is what makes the boundaries decidable rather than negotiable. Ares
owns solid momentum balance; Proteus closes it; the pair is not a matter of
taste.

## What Ares depends on

Strictly inward, and only on vocabulary and infrastructure:

| Dependency | For |
| --- | --- |
| `aequitas` | typed physical quantities — a stress cannot be assigned to a length |
| `eunomia` | scalar and numeric-trait vocabulary |
| `proteus` | the constitutive closure |
| `athena`, `leto` | the solve (in `ares-operator` only) |
| `harmonia` | coupling (in `ares-coupling` only) |

The domain crate has the first three and nothing else. It is `no_std`, forbids
`unsafe`, and allocates nothing.

## What Ares must never depend on

**No other balance domain.** ADR 0055 rule R7 forbids a direct dependency edge
between two balance domains and routes coupling through Harmonia. So Ares has
no edge to CFDrs, kwavers, helios, or asclepius — in any phase, not just this
one.

The rule is mechanically enforced rather than reviewed: the Atlas conformance
scan classifies every direct runtime edge and counts violations. `ares` was the
first entry in that rule's boundary list, which is when the rule stopped being
preventive and started checking something.

**No mesh generation.** Gaia owns meshes, geometry, and proximity queries.
`SimplexMesh` is a borrowed view over coordinates and connectivity — the shape
assembly reads, not a mesh representation competing with Gaia's.

**No solver policy.** Athena owns iteration, preconditioning, and convergence.
`ares-operator` presents the assembled operator in the shape Athena consumes and
makes no decision about how it is solved.

**No time integration.** Horae owns time. Phase 0 is static and has none; when
dynamics arrive, they arrive through Horae rather than through a loop written
here.

## The three crates and why they are three

| Crate | Links `std`? | Depends on |
| --- | --- | --- |
| `ares-solid` | no | aequitas, eunomia, proteus |
| `ares-operator` | yes | ares, athena, leto |
| `ares-coupling` | yes | ares, ares-operator, harmonia |

The split follows from Athena's trait, as chapter 8 explains: its operator seam
can only be implemented against a *named* backend, and the only host backend
links `std`. Putting the seam in the domain core would push a solver dependency
into the balance code.

A cargo feature was the obvious alternative and was rejected on stronger
grounds than taste: it would make the shipped configuration the one CI does not
build by default, and a feature-gated solver path is an untested path.
[ADR 0001](https://github.com/ryancinsight/ares/blob/main/docs/adr/0001-athena-seam-as-a-separate-crate.md)
records the decision.

Dependencies run strictly inward. `ares-coupling` depends on `ares`; nothing
depends back.

## Registry naming

The crates.io name `ares` belongs to an unrelated third-party crate, so the
core publishes as **`ares-solid`** while the import path stays `ares` through
`[lib] name`:

```toml
[dependencies]
ares = { package = "ares-solid", version = "0.1.0" }
```

The convention `proteus-mat` and `gaia-mesh` follow. It is worth knowing about
because tooling that reads a dependency-table *key* rather than the resolved
package name gets this wrong — and in this stack, twice: once in the
architecture test's boundary set and once in the publish-order tool, where it
placed `ares-solid` a wave ahead of a crate it depends on.

## What comes next

Phase 1 and beyond, each with its own charter and none scaffolded here:

- **Dynamics** — mass, inertia, modal analysis, through Horae.
- **Plasticity** — through a Proteus closure that carries state.
- **Finite deformation** — where chapter 1's small-strain assumption ends.
- **Contact** — with Gaia supplying proximity.
- **Two-way FSI** — waiting on ALE in the fluid solver, per chapter 9.
