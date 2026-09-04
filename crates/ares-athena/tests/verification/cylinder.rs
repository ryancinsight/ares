//! Lame's thick-walled cylinder under internal pressure.
//!
//! # What this oracle adds that the others cannot
//!
//! It is the only one here with a **curved boundary**. Every other fixture is
//! a rectangle whose edges the mesh represents exactly, so none of them can
//! detect a defect that appears only when element edges approximate a surface
//! — a traction resolved along the wrong normal, or a facet measure that is
//! right for an axis-aligned edge and wrong for an oblique one.
//!
//! It is also axisymmetric, so it exercises coupling between normal components
//! that a uniaxial fixture leaves at zero.
//!
//! # The closed form
//!
//! For a cylinder of inner radius `a` and outer radius `b` under internal
//! pressure `p`, in plane strain:
//!
//! ```text
//! sigma_r(r)  = (a^2 p / (b^2 - a^2)) (1 - b^2 / r^2)
//! sigma_th(r) = (a^2 p / (b^2 - a^2)) (1 + b^2 / r^2)
//! u_r(r)      = (1 + nu) a^2 p / (E (b^2 - a^2)) [ (1 - 2 nu) r + b^2 / r ]
//! ```
//!
//! # Why the headline assertion is a rate rather than a tolerance
//!
//! The mesh's straight edges approximate the circular boundary, so the
//! geometry carries an error of its own on top of the discretisation's. A
//! fixed tolerance on a single mesh would be a fitted number standing in for
//! two effects at once. Refinement separates them: both are second order, so
//! the total must fall as `O(h^2)`, and that is a statement about the method
//! rather than about the mesh that happened to be chosen.

use ares::{PrescribedDisplacement, SimplexMesh, TractionBoundary, TractionFacet};

use super::mesh::{moduli, solve};

const INNER: f64 = 1.0;
const OUTER: f64 = 2.0;
const PRESSURE: f64 = 5.0e7;
const YOUNG: f64 = 200e9;
const POISSON: f64 = 0.3;

/// Exact radial displacement at radius `r`.
fn exact_radial_displacement(radius: f64) -> f64 {
    let squares = OUTER * OUTER - INNER * INNER;
    (1.0 + POISSON) * INNER * INNER * PRESSURE / (YOUNG * squares)
        * ((1.0 - 2.0 * POISSON) * radius + OUTER * OUTER / radius)
}

/// A quarter annulus meshed on a polar grid over `0 <= theta <= pi/2`.
struct Annulus {
    nodes: Vec<[f64; 2]>,
    cells: Vec<[usize; 3]>,
    radial: usize,
    angular: usize,
}

impl Annulus {
    fn new(radial: usize, angular: usize) -> Self {
        let mut nodes = Vec::with_capacity((radial + 1) * (angular + 1));
        for ring in 0..=radial {
            let r = INNER + (OUTER - INNER) * ring as f64 / radial as f64;
            for spoke in 0..=angular {
                let theta = core::f64::consts::FRAC_PI_2 * spoke as f64 / angular as f64;
                nodes.push([r * theta.cos(), r * theta.sin()]);
            }
        }
        let index = |ring: usize, spoke: usize| ring * (angular + 1) + spoke;
        let mut cells = Vec::with_capacity(radial * angular * 2);
        for ring in 0..radial {
            for spoke in 0..angular {
                let (a, b) = (index(ring, spoke), index(ring, spoke + 1));
                let (c, d) = (index(ring + 1, spoke + 1), index(ring + 1, spoke));
                // In (x, y) the quad a -> b -> c -> d runs clockwise, because
                // increasing theta moves along the arc while increasing radius
                // moves outward. The triangles therefore take the reverse
                // order. The mesh constructor rejects an inverted cell, so a
                // mistake here fails loudly rather than quietly negating a
                // cell's stiffness.
                cells.push([a, d, c]);
                cells.push([a, c, b]);
            }
        }
        Self {
            nodes,
            cells,
            radial,
            angular,
        }
    }

    fn index(&self, ring: usize, spoke: usize) -> usize {
        ring * (self.angular + 1) + spoke
    }

    fn degrees_of_freedom(&self) -> usize {
        self.nodes.len() * 2
    }

    /// Element size: the larger of the radial and circumferential spacings.
    fn cell_size(&self) -> f64 {
        let radial = (OUTER - INNER) / self.radial as f64;
        let circumferential = OUTER * core::f64::consts::FRAC_PI_2 / self.angular as f64;
        radial.max(circumferential)
    }
}

/// Solve the quarter cylinder and return the displacement field.
fn solve_field(annulus: &Annulus) -> Vec<f64> {
    let mesh = SimplexMesh::try_new(&annulus.nodes, &annulus.cells).expect("valid annulus");

    // Symmetry. The theta = 0 cut lies on the x-axis and cannot move in y; the
    // theta = pi/2 cut lies on the y-axis and cannot move in x. Together they
    // remove the rigid-body modes while prescribing no radial motion at all,
    // so the displacement compared below is entirely the solve's.
    let mut constrained: Vec<(usize, f64)> = Vec::new();
    for ring in 0..=annulus.radial {
        constrained.push((annulus.index(ring, 0) * 2 + 1, 0.0));
        constrained.push((annulus.index(ring, annulus.angular) * 2, 0.0));
    }
    constrained.sort_unstable_by_key(|(dof, _)| *dof);
    let prescribed: Vec<PrescribedDisplacement<f64>> = constrained
        .into_iter()
        .map(|(dof, value)| PrescribedDisplacement::new(dof / 2, dof % 2, value))
        .collect();

    // Internal pressure. The solid's outward normal at r = a points toward the
    // origin, so the traction on the solid is `+p r_hat`, directed outward.
    let mut facets = Vec::new();
    for spoke in 0..annulus.angular {
        let (first, second) = (annulus.index(0, spoke), annulus.index(0, spoke + 1));
        let midpoint = [
            f64::midpoint(annulus.nodes[first][0], annulus.nodes[second][0]),
            f64::midpoint(annulus.nodes[first][1], annulus.nodes[second][1]),
        ];
        let length = midpoint[0].hypot(midpoint[1]);
        let outward = [midpoint[0] / length, midpoint[1] / length];
        facets.push(TractionFacet::new(
            [first, second],
            [PRESSURE * outward[0], PRESSURE * outward[1]],
        ));
    }
    let boundary = TractionBoundary::try_new(&facets, &annulus.nodes).expect("valid facets");
    let mut external = vec![0.0_f64; annulus.degrees_of_freedom()];
    boundary
        .add_consistent_loads(&annulus.nodes, &mut external)
        .expect("well-shaped");

    solve(&mesh, moduli(YOUNG, POISSON), &prescribed, &external)
}

/// Solve at a resolution; return `(h, relative error in u_r, mean u_r)`.
fn solve_at(radial: usize, angular: usize) -> (f64, f64, f64) {
    let annulus = Annulus::new(radial, angular);
    let displacement = solve_field(&annulus);

    let mut error = 0.0_f64;
    let mut truth = 0.0_f64;
    let mut mean = 0.0_f64;
    for (node, position) in annulus.nodes.iter().enumerate() {
        let radius = position[0].hypot(position[1]);
        let radial_component = (displacement[node * 2] * position[0]
            + displacement[node * 2 + 1] * position[1])
            / radius;
        let want = exact_radial_displacement(radius);
        error += (radial_component - want).powi(2);
        truth += want * want;
        mean += radial_component;
    }
    (
        annulus.cell_size(),
        (error / truth).sqrt(),
        mean / annulus.nodes.len() as f64,
    )
}

/// Mean radial displacement on the inner and outer walls.
fn wall_displacements(radial: usize, angular: usize) -> (f64, f64) {
    let annulus = Annulus::new(radial, angular);
    let displacement = solve_field(&annulus);
    let mean_on = |ring: usize| -> f64 {
        let mut total = 0.0_f64;
        for spoke in 0..=annulus.angular {
            let node = annulus.index(ring, spoke);
            let position = annulus.nodes[node];
            let radius = position[0].hypot(position[1]);
            total += (displacement[node * 2] * position[0]
                + displacement[node * 2 + 1] * position[1])
                / radius;
        }
        total / (annulus.angular + 1) as f64
    };
    (mean_on(0), mean_on(annulus.radial))
}

#[test]
fn the_thick_walled_cylinder_matches_lame() {
    let (_, error, _) = solve_at(12, 18);
    assert!(
        error < 0.02,
        "the radial displacement differs from Lame by {error:.4} relative"
    );
}

#[test]
fn the_cylinder_expands_under_internal_pressure() {
    // Direction, asserted on its own. An inward traction would still converge,
    // still be smooth, and still refine cleanly; it differs from the closed
    // form only by a sign, which a relative error norm reports as roughly two
    // rather than as the reversal it is.
    let (_, _, mean) = solve_at(8, 12);
    assert!(
        mean > 0.0,
        "internal pressure moved the wall inward by {mean:.6e}"
    );

    // And the inner wall moves further than the outer one. That is a real
    // feature of the closed form rather than a restatement of the sign:
    // `du_r/dr = C[(1 - 2 nu) - b^2 / r^2]` is negative across the whole wall
    // for these radii, so a solve that got the pressure onto the wrong surface
    // would reverse this gradient while keeping every displacement positive.
    assert!(
        exact_radial_displacement(INNER) > exact_radial_displacement(OUTER),
        "the closed form does not decay outward, so this check is vacuous"
    );
    let (inner, outer) = wall_displacements(8, 12);
    assert!(
        inner > outer,
        "the outer wall ({outer:.6e}) moved further than the inner ({inner:.6e}), which means          the pressure landed on the wrong surface"
    );
}

#[test]
fn the_cylinder_converges_at_second_order() {
    let mut previous: Option<(f64, f64)> = None;
    let mut rates = Vec::new();
    for (radial, angular) in [(4_usize, 6_usize), (8, 12), (16, 24)] {
        let (size, error, _) = solve_at(radial, angular);
        if let Some((last_size, last_error)) = previous {
            rates.push((last_error / error).ln() / (last_size / size).ln());
        }
        previous = Some((size, error));
    }
    for (index, rate) in rates.iter().enumerate() {
        assert!(
            *rate > 1.7,
            "refinement step {index} converged at order {rate:.3}, below second order"
        );
    }
}
