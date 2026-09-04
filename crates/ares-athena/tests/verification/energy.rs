//! Strain energy equals external work.
//!
//! # Why this oracle earns its place beside the others
//!
//! It is a conservation identity, so it holds **exactly** on any mesh at any
//! resolution — unlike every closed-form comparison here, which holds only in
//! the limit. That makes it the one oracle that separates a solve which is
//! merely inaccurate from one which is inconsistent: a coarse mesh gives a
//! poor cantilever deflection while satisfying this identity to rounding, and
//! a defect in assembly or in the constrained operator breaks the identity at
//! any resolution.
//!
//! For a linear elastic static solve, the stored strain energy is
//! `U = (1/2) u . K u` and the work done by the applied load is
//! `W = (1/2) u . f`. Equilibrium is `K u = f`, so `U = W` follows directly —
//! which is exactly why it is a check on the *assembly* rather than on
//! elasticity: it fails when the operator that was solved is not the operator
//! the energy is computed from.

use ares::{PrescribedDisplacement, SimplexMesh, TractionBoundary, TractionFacet};

use super::mesh::{Grid, moduli, solve};

#[test]
fn strain_energy_equals_external_work() {
    let grid = Grid::new(2.0, 1.0, 6, 3);
    let mesh = SimplexMesh::try_new(&grid.nodes, &grid.cells).expect("valid grid");
    let material = moduli(200e9, 0.3);

    // Clamp x = 0, pull the far edge in +x and +y at once so the field is not
    // a pure extension.
    let mut prescribed = Vec::new();
    for row in 0..=grid.rows {
        let node = grid.node_index(0, row);
        prescribed.push(PrescribedDisplacement::new(node, 0, 0.0));
        prescribed.push(PrescribedDisplacement::new(node, 1, 0.0));
    }

    let mut facets = Vec::new();
    for row in 0..grid.rows {
        facets.push(TractionFacet::new(
            [
                grid.node_index(grid.columns, row),
                grid.node_index(grid.columns, row + 1),
            ],
            [4.0e6_f64, 1.5e6],
        ));
    }
    let boundary = TractionBoundary::try_new(&facets, &grid.nodes).expect("valid facets");
    let mut external = vec![0.0_f64; grid.degrees_of_freedom()];
    boundary
        .add_consistent_loads(&grid.nodes, &mut external)
        .expect("well-shaped");

    let displacement = solve(&mesh, material, &prescribed, &external);

    // U = (1/2) u . K u, from the unconstrained operator: the stored energy is
    // a property of the body, not of how it was held.
    let mut internal = vec![0.0_f64; grid.degrees_of_freedom()];
    mesh.internal_forces(&material, &displacement, &mut internal)
        .expect("well-shaped");
    let strain_energy: f64 = displacement
        .iter()
        .zip(internal.iter())
        .map(|(u, f)| u * f)
        .sum::<f64>()
        / 2.0;

    // W = (1/2) u . f_ext. The constrained nodes contribute nothing because
    // they did not move, so the reactions there do no work — which is the
    // physical content of a rigid support.
    let external_work: f64 = displacement
        .iter()
        .zip(external.iter())
        .map(|(u, f)| u * f)
        .sum::<f64>()
        / 2.0;

    assert!(
        strain_energy > 0.0,
        "the body stored no energy, so this identity is being satisfied trivially"
    );
    // The bound is the solver's, not the discretisation's: the identity is
    // exact for the exact solution of the discrete system, so the residual
    // tolerance is what separates the two sides.
    let relative = (strain_energy - external_work).abs() / strain_energy;
    assert!(
        relative < 1e-9,
        "strain energy {strain_energy:.6e} and external work {external_work:.6e} differ by \
         {relative:.3e} relative, far above the solver tolerance"
    );
}
