# 9. Coupling to a fluid

A pipe carrying pressurised flow deforms. The fluid pushes on the wall, the
wall moves, and if it moves enough the flow changes. Computing that is
**fluid-structure interaction**, and this chapter covers the half of it Phase 0
does.

## One way, and why that is a design choice rather than a limitation

Full FSI is two-way: the fluid loads the solid, the solid's motion changes the
fluid domain, and the two are iterated to consistency. That requires the fluid
solver to handle a **moving mesh** — an arbitrary Lagrangian-Eulerian
formulation. CFDrs has none.

Phase 0 therefore couples one way: traction flows fluid to solid, the solid
deforms, and the deformation is exported but does not move the fluid mesh. That
is a real physical assumption — it holds when the deformation is small enough
not to change the flow appreciably, which covers a stiff wall under moderate
pressure and does not cover a flapping membrane.

Stating it as a limitation is the point. Building a two-way interface that
silently ignored the return path would be a mock of the thing it claims to be.

## The partition

Harmonia drives coupled problems by asking each side to advance and to export
its interface values. `ares-harmonia` presents the structural solve in that
shape.

```rust,ignore
use ares_harmonia::{StructuralInterface, StructuralPartition};

let interface = StructuralInterface::try_new(&interface_nodes, &facets, &mesh)?;
let mut partition = StructuralPartition::try_new(mesh, operator, interface, policy)?;

partition.solve_for_traction(&mut state, &traction)?;
partition.export(&state, &mut displacement)?;
```

Neither side depends on the other. CFDrs computes the interface traction from
its own flow state; the structural side consumes it; a Harmonia driver joins
them. Atlas ADR 0055 forbids a direct dependency between two balance domains
and routes coupling through Harmonia, which is exactly this shape.

## The exchange orderings differ, and that is not an oversight

- **traction** is facet index major, component minor;
- **displacement** is node index major, component minor.

They differ because the quantities live on different entities. A traction is a
stress resolved on a **surface**, and the fluid side computes one per face from
the flow state either side of it. A displacement is a property of a **point**.

Forcing either onto the other's entity would mean interpolating, and Harmonia
explicitly does not interpolate — its transfers are index maps. Inventing an
interpolation scheme inside a coupling adapter would put a scheme with its own
conservation properties somewhere nobody owns it.

## Interface work is conserved, exactly

The headline property: the work the fluid does on the interface equals twice
the strain energy the structure stores.

It is exact, and the reason is worth following. With traction constant per
facet and displacement linear over it:

```text
integral(t . u) dS = sum over facets of t_f |A_f| (mean u over the facet nodes)
                   = sum over facets, nodes of t_f |A_f| u_a / D
                   = sum over nodes of u_a . f_a
```

and `f_a = t_f |A_f| / D` is precisely the consistent nodal load of chapter 7.
The identity holds **because** the load is the consistent one. A lumped load
carries the same resultant force and breaks it — which is not a hypothesis: the
test suite was mutation-checked against exactly that substitution, and the work
identity is the only test that fails.

The two sides are computed by routes sharing no arithmetic — a facet integral
over the traction exchange against a stiffness assembly over every cell in the
mesh — so their agreement is evidence rather than a restatement.

## Conformity

Phase 0 requires the two interface discretisations to **conform**: the fluid's
faces and the solid's facets must be the same surface with the same nodes. A
non-conforming interface is rejected with a typed error rather than transferred
approximately.

Non-conforming transfer is a genuine body of work with its own conservation
properties to prove, and it belongs to whoever owns interpolation. Today nobody
does. Inventing it inside Ares would put fluid-side knowledge in a solid
solver; inventing it inside Harmonia would contradict what Harmonia is.
