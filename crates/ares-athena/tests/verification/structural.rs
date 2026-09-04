//! Cantilever tip deflection against beam theory.
//!
//! # Why beam theory rather than another elasticity result
//!
//! Every other oracle here is elasticity checking elasticity. Euler-Bernoulli
//! beam theory is an independent structural model with its own assumptions, so
//! agreement is evidence across theories rather than internal consistency.
//!
//! # Approached from the known direction
//!
//! Linear triangles are stiff in bending. Their displacement field is linear,
//! so their strain is constant per element, and a bending state — whose strain
//! varies linearly through the depth — cannot be represented within one
//! element. The computed deflection is therefore **below** the beam value and
//! rises toward it under refinement.
//!
//! That direction is the assertion. A result that matched the beam value on a
//! coarse mesh, or overshot it, would be evidence of a defect rather than of
//! accuracy — the most likely being an element that is too soft because its
//! stiffness lost a factor.
//!
//! The assertions are therefore structural rather than a tuned threshold: the
//! deflection is below the beam value, it increases monotonically under
//! refinement, and the gap closes. None of those needs a fitted constant, and
//! each fails for a different defect.
//!
//! The exact two-dimensional answer slightly *exceeds* Euler-Bernoulli, by a
//! shear term of order `(H/L)^2` — about a percent at the aspect ratio used
//! here — so "below" is a statement about element stiffness dominating at
//! these resolutions, not a universal bound.

use ares::{PrescribedDisplacement, SimplexMesh, TractionBoundary, TractionFacet};

use super::mesh::{Grid, moduli, solve};

const LENGTH: f64 = 10.0;
const DEPTH: f64 = 1.0;
const YOUNG: f64 = 200e9;
const POISSON: f64 = 0.3;
const LOAD: f64 = 1.0e5;

/// `delta = P L^3 / (3 E' I)`, with `I = H^3 / 12` per unit thickness and
/// `E' = E / (1 - nu^2)` for plane strain.
fn beam_deflection() -> f64 {
    let inertia = DEPTH.powi(3) / 12.0;
    let plane_strain_modulus = YOUNG / (1.0 - POISSON * POISSON);
    LOAD * LENGTH.powi(3) / (3.0 * plane_strain_modulus * inertia)
}

/// Mean tip deflection for a mesh of `columns x rows` cells.
fn tip_deflection(columns: usize, rows: usize) -> f64 {
    let grid = Grid::new(LENGTH, DEPTH, columns, rows);
    let mesh = SimplexMesh::try_new(&grid.nodes, &grid.cells).expect("valid grid");

    // Clamp the root.
    let mut prescribed = Vec::new();
    for row in 0..=grid.rows {
        let node = grid.node_index(0, row);
        prescribed.push(PrescribedDisplacement::new(node, 0, 0.0));
        prescribed.push(PrescribedDisplacement::new(node, 1, 0.0));
    }

    // The tip load as a uniform shear traction over the end face, so the total
    // is `LOAD` regardless of how finely the face is divided.
    let traction = [0.0_f64, -LOAD / DEPTH];
    let mut facets = Vec::new();
    for row in 0..grid.rows {
        facets.push(TractionFacet::new(
            [
                grid.node_index(grid.columns, row),
                grid.node_index(grid.columns, row + 1),
            ],
            traction,
        ));
    }
    let boundary = TractionBoundary::try_new(&facets, &grid.nodes).expect("valid facets");
    let mut external = vec![0.0_f64; grid.degrees_of_freedom()];
    boundary
        .add_consistent_loads(&grid.nodes, &mut external)
        .expect("well-shaped");

    let displacement = solve(&mesh, moduli(YOUNG, POISSON), &prescribed, &external);

    // Averaged over the tip face: the end section rotates, so a single corner
    // node would report a deflection contaminated by that rotation.
    let mut total = 0.0_f64;
    for row in 0..=grid.rows {
        total -= displacement[grid.node_index(grid.columns, row) * 2 + 1];
    }
    total / (grid.rows + 1) as f64
}

#[test]
fn the_cantilever_approaches_beam_theory_from_below() {
    let beam = beam_deflection();
    let deflections: Vec<f64> = [(20_usize, 2_usize), (40, 4), (80, 8)]
        .into_iter()
        .map(|(columns, rows)| tip_deflection(columns, rows))
        .collect();

    for (index, deflection) in deflections.iter().enumerate() {
        assert!(
            *deflection > 0.0,
            "mesh {index} deflected {deflection:.6e}, the wrong way or not at all"
        );
        assert!(
            *deflection < beam,
            "mesh {index} deflected {deflection:.6e}, above the beam value {beam:.6e}; linear \
             triangles are stiff in bending, so an element this soft has lost stiffness"
        );
    }

    for pair in deflections.windows(2) {
        assert!(
            pair[1] > pair[0],
            "refinement made the beam stiffer ({:.6e} then {:.6e}), the opposite of convergence",
            pair[0],
            pair[1]
        );
    }

    // The gap closes. Without this, a monotone sequence converging to the
    // wrong value would pass everything above.
    let first_gap = beam - deflections[0];
    let last_gap = beam - deflections[deflections.len() - 1];
    assert!(
        last_gap < first_gap / 2.0,
        "the gap to beam theory closed only from {first_gap:.6e} to {last_gap:.6e}"
    );

    // A loose floor against gross error, well below where any of the meshes
    // above land. This is the one number here that is not structural, and it
    // is deliberately far from the measured values so it catches an order-of
    // -magnitude defect rather than encoding a result.
    let finest = deflections[deflections.len() - 1];
    assert!(
        finest > beam * 0.5,
        "the finest mesh reached only {:.1}% of beam theory",
        100.0 * finest / beam
    );
}
