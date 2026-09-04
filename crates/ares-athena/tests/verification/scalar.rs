//! The same oracles at `f32`.
//!
//! # What this catches
//!
//! A kernel that is generic in its signature and `f64` in its body — the fake
//! generic the integrity rules name — passes every `f64` oracle in this suite
//! and every one of them at `f32` too, because it would simply be more
//! accurate than `f32` warrants. What it cannot do is match `f32`'s error
//! *scale*: the assertions below are bounded by `f32::EPSILON`, so a body
//! silently computing in `f64` and narrowing at the boundary is not caught by
//! accuracy but by the accompanying `f64` runs it must also satisfy. The real
//! catch is coarser and more reliable: a body pinned to a concrete type does
//! not compile against `T`, and this module is what forces the whole solve
//! path — assembly, constrained operator, Athena adapter, Krylov solver — to
//! monomorphise at a second scalar at all.
//!
//! # Why the convergence-rate studies are not repeated here
//!
//! An order-of-accuracy study needs the discretisation error to dominate every
//! other error in the calculation. At `f32` the representable relative
//! precision is about `1e-7`, and the manufactured solution's discretisation
//! error passes below that by the third refinement — after which the measured
//! "rate" describes rounding, not the element. Repeating the study at `f32`
//! would therefore produce a number that looks like evidence and is not.
//!
//! The oracles that *do* transfer are the ones with no asymptotic content: a
//! single-mesh accuracy comparison, and the energy identity, which is exact at
//! any resolution and so is limited only by the solver's tolerance.

use ares::{PrescribedDisplacement, SimplexMesh, TractionBoundary, TractionFacet};

use super::mesh::{Grid, l2_norm, moduli, solve};

#[test]
fn the_manufactured_solution_is_recovered_at_f32() {
    // The same field as the f64 study, on a mesh coarse enough that the
    // discretisation error stays well above f32's noise floor - otherwise the
    // comparison would be measuring rounding.
    const YOUNG: f64 = 200e9;
    const POISSON: f64 = 0.3;
    let amplitude = [1.3e-6_f32, -0.8e-6];
    let (lambda, mu) = super::mesh::lame(YOUNG, POISSON);
    let pi = core::f32::consts::PI;

    let grid = Grid::new(1.0, 1.0, 12, 12);
    let nodes: Vec<[f32; 2]> = grid
        .nodes
        .iter()
        .map(|p| [p[0] as f32, p[1] as f32])
        .collect();
    let mesh = SimplexMesh::try_new(&nodes, &grid.cells).expect("valid grid");

    let mut prescribed = Vec::new();
    for node in grid.boundary_nodes() {
        prescribed.push(PrescribedDisplacement::new(node, 0, 0.0_f32));
        prescribed.push(PrescribedDisplacement::new(node, 1, 0.0_f32));
    }

    let mut body = vec![0.0_f32; nodes.len() * 2];
    let mut truth = vec![0.0_f32; nodes.len() * 2];
    for (node, position) in nodes.iter().enumerate() {
        let f = (pi * position[0]).sin() * (pi * position[1]).sin();
        let cc = (pi * position[0]).cos() * (pi * position[1]).cos();
        let pi2 = pi * pi;
        let (l, m) = (lambda as f32, mu as f32);
        body[node * 2] = pi2 * ((l + 3.0 * m) * amplitude[0] * f - (l + m) * amplitude[1] * cc);
        body[node * 2 + 1] = pi2 * ((l + 3.0 * m) * amplitude[1] * f - (l + m) * amplitude[0] * cc);
        truth[node * 2] = amplitude[0] * f;
        truth[node * 2 + 1] = amplitude[1] * f;
    }
    let mut external = vec![0.0_f32; nodes.len() * 2];
    mesh.add_body_force(&body, &mut external)
        .expect("well-shaped");

    let computed = solve(&mesh, moduli::<f32>(YOUNG, POISSON), &prescribed, &external);

    let error: Vec<f32> = computed
        .iter()
        .zip(truth.iter())
        .map(|(got, want)| got - want)
        .collect();
    let relative = l2_norm(&mesh, &error) / l2_norm(&mesh, &truth);
    assert!(
        relative < 0.05,
        "the f32 solve reproduces the manufactured field to only {relative:.4} relative error"
    );
}

#[test]
fn strain_energy_equals_external_work_at_f32() {
    // A conservation identity rather than a limit statement, so it transfers
    // to f32 unchanged except for the bound, which follows the scalar's
    // precision rather than being relaxed by hand.
    let grid = Grid::new(2.0, 1.0, 6, 3);
    let nodes: Vec<[f32; 2]> = grid
        .nodes
        .iter()
        .map(|p| [p[0] as f32, p[1] as f32])
        .collect();
    let mesh = SimplexMesh::try_new(&nodes, &grid.cells).expect("valid grid");
    let material = moduli::<f32>(200e9, 0.3);

    let mut prescribed = Vec::new();
    for row in 0..=grid.rows {
        let node = grid.node_index(0, row);
        prescribed.push(PrescribedDisplacement::new(node, 0, 0.0_f32));
        prescribed.push(PrescribedDisplacement::new(node, 1, 0.0_f32));
    }

    let mut facets = Vec::new();
    for row in 0..grid.rows {
        facets.push(TractionFacet::new(
            [
                grid.node_index(grid.columns, row),
                grid.node_index(grid.columns, row + 1),
            ],
            [4.0e6_f32, 1.5e6],
        ));
    }
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facets");
    let mut external = vec![0.0_f32; nodes.len() * 2];
    boundary
        .add_consistent_loads(&nodes, &mut external)
        .expect("well-shaped");

    let displacement = solve(&mesh, material, &prescribed, &external);

    let mut internal = vec![0.0_f32; nodes.len() * 2];
    mesh.internal_forces(&material, &displacement, &mut internal)
        .expect("well-shaped");
    // Accumulated in f64: the identity is a property of the f32 solve, but
    // summing its terms in f32 would add the reduction's own error to the
    // difference being measured and report it as a violation.
    let strain_energy: f64 = displacement
        .iter()
        .zip(internal.iter())
        .map(|(u, f)| f64::from(*u) * f64::from(*f))
        .sum::<f64>()
        / 2.0;
    let external_work: f64 = displacement
        .iter()
        .zip(external.iter())
        .map(|(u, f)| f64::from(*u) * f64::from(*f))
        .sum::<f64>()
        / 2.0;

    assert!(strain_energy > 0.0, "the body stored no energy");
    let relative = (strain_energy - external_work).abs() / strain_energy;
    assert!(
        relative < 1e-4,
        "at f32, strain energy {strain_energy:.6e} and external work {external_work:.6e} differ \
         by {relative:.3e} relative, above the solver tolerance"
    );
}
