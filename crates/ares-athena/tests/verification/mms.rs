//! Manufactured solution and the convergence-order study.
//!
//! # Why a manufactured solution
//!
//! Every other oracle here depends on a special geometry — a slender beam, an
//! axisymmetric annulus — and so tests the solver only where that geometry's
//! assumptions hold. A manufactured solution has no such dependence: any
//! smooth field is admissible, and the body force that produces it follows
//! from the governing equation by differentiation.
//!
//! # The field and its body force
//!
//! On the unit square, with `f = sin(pi x) sin(pi y)`:
//!
//! ```text
//! u = (A f, B f)
//! ```
//!
//! which vanishes on the whole boundary, so the Dirichlet data is homogeneous
//! and no boundary term contaminates the interior error.
//!
//! For linear isotropic elasticity `div sigma = (l + m) grad(div u) + m lap u`,
//! so `b = -div sigma` works out to
//!
//! ```text
//! b_x = pi^2 [ (l + 3m) A f - (l + m) B cc ]
//! b_y = pi^2 [ (l + 3m) B f - (l + m) A cc ]
//! ```
//!
//! with `cc = cos(pi x) cos(pi y)`. Derived here rather than quoted: `div u =
//! pi (A cos(pi x) sin(pi y) + B sin(pi x) cos(pi y))`, whose gradient
//! contributes `pi^2 (-A f + B cc, -B f + A cc)`, and `lap u = -2 pi^2 u`.

use ares::SimplexMesh;

use super::mesh::{Grid, l2_norm, lame, moduli, solve};
use ares::PrescribedDisplacement;

const YOUNG: f64 = 200e9;
const POISSON: f64 = 0.3;
const AMPLITUDE: [f64; 2] = [1.3e-6, -0.8e-6];

fn exact(position: &[f64; 2]) -> [f64; 2] {
    let f =
        (core::f64::consts::PI * position[0]).sin() * (core::f64::consts::PI * position[1]).sin();
    [AMPLITUDE[0] * f, AMPLITUDE[1] * f]
}

fn body_force(position: &[f64; 2]) -> [f64; 2] {
    let (lambda, mu) = lame(YOUNG, POISSON);
    let pi = core::f64::consts::PI;
    let f = (pi * position[0]).sin() * (pi * position[1]).sin();
    let cc = (pi * position[0]).cos() * (pi * position[1]).cos();
    let pi2 = pi * pi;
    [
        pi2 * ((lambda + 3.0 * mu) * AMPLITUDE[0] * f - (lambda + mu) * AMPLITUDE[1] * cc),
        pi2 * ((lambda + 3.0 * mu) * AMPLITUDE[1] * f - (lambda + mu) * AMPLITUDE[0] * cc),
    ]
}

/// Solve the manufactured problem on an `n x n` grid and return
/// `(cell size, L2 error, L2 norm of the exact field)`.
fn solve_at(divisions: usize) -> (f64, f64, f64) {
    let grid = Grid::new(1.0, 1.0, divisions, divisions);
    let mesh = SimplexMesh::try_new(&grid.nodes, &grid.cells).expect("valid grid");

    // The field vanishes on the boundary, so every boundary node is held at
    // zero. Conditions must be strictly increasing in degree of freedom.
    let mut prescribed = Vec::new();
    for node in grid.boundary_nodes() {
        prescribed.push(PrescribedDisplacement::new(node, 0, 0.0));
        prescribed.push(PrescribedDisplacement::new(node, 1, 0.0));
    }

    let mut nodal_body_force = vec![0.0_f64; grid.degrees_of_freedom()];
    for (node, position) in grid.nodes.iter().enumerate() {
        let value = body_force(position);
        nodal_body_force[node * 2] = value[0];
        nodal_body_force[node * 2 + 1] = value[1];
    }
    let mut external = vec![0.0_f64; grid.degrees_of_freedom()];
    mesh.add_body_force(&nodal_body_force, &mut external)
        .expect("well-shaped");

    let computed = solve(&mesh, moduli(YOUNG, POISSON), &prescribed, &external);

    let mut error = vec![0.0_f64; grid.degrees_of_freedom()];
    let mut truth = vec![0.0_f64; grid.degrees_of_freedom()];
    for (node, position) in grid.nodes.iter().enumerate() {
        let want = exact(position);
        for component in 0..2 {
            truth[node * 2 + component] = want[component];
            error[node * 2 + component] = computed[node * 2 + component] - want[component];
        }
    }
    (
        grid.cell_size(),
        l2_norm(&mesh, &error),
        l2_norm(&mesh, &truth),
    )
}

#[test]
fn the_manufactured_solution_is_recovered() {
    // A single mesh, checked for relative accuracy before any rate is claimed.
    // A rate computed between two equally wrong fields is still a rate.
    let (_, error, truth) = solve_at(16);
    let relative = error / truth;
    assert!(
        relative < 0.02,
        "the manufactured field is reproduced to only {relative:.4} relative error"
    );
    assert!(truth > 0.0, "the exact field is identically zero");
}

#[test]
fn h_refinement_recovers_second_order_convergence() {
    // The rate that certifies the element. A defect leaving the solve
    // plausible but inconsistent - a mis-scaled shape gradient, a wrong
    // element measure - caps this at first order while every single-mesh
    // check above still passes.
    //
    // Substituting a lumped body force does *not* cap it, which was measured
    // rather than assumed: the vertex rule is also exact for a linear
    // integrand, so it is second order too. The consistent load is chosen for
    // being the Galerkin integral, not for the order.
    //
    // # Why the assertion is the shape of the approach, not a floor
    //
    // Second-order convergence is an asymptotic statement, and the coarsest
    // mesh here is not in the asymptotic regime: at `h = 0.35` a single sine
    // hump spans four elements. The measured rates are 1.65, 1.88, 1.97 -
    // monotonically approaching 2 from below, which is the signature of a
    // second-order method observed from too far out, not of a first-order one.
    //
    // A floor of 1.8 on every step would have failed on that first pair and
    // said nothing about why. Asserting instead that the rates rise toward 2
    // and that the finest exceeds 1.9 uses every data point and is strictly
    // harder to satisfy by accident: a genuinely first-order method gives
    // rates flat near 1.0, which fails both halves, and no amount of starting
    // finer would rescue it.
    let mut previous: Option<(f64, f64)> = None;
    let mut rates = Vec::new();
    let mut trace = Vec::new();
    for divisions in [4_usize, 8, 16, 32] {
        let (size, error, _) = solve_at(divisions);
        if let Some((last_size, last_error)) = previous {
            rates.push((last_error / error).ln() / (last_size / size).ln());
        }
        trace.push((divisions, size, error));
        previous = Some((size, error));
    }
    let report = format!("{trace:?} rates {rates:?}");

    for pair in rates.windows(2) {
        assert!(
            pair[1] > pair[0],
            "the convergence rate fell under refinement, so it is not approaching a limit:              {report}"
        );
    }

    let finest = *rates.last().expect("four meshes give three rates");
    assert!(
        finest > 1.9,
        "the finest refinement converged at order {finest:.3}, short of second order: {report}"
    );
    assert!(
        finest < 2.5,
        "order {finest:.3} is too high to be second-order convergence; the error has probably          reached the solver tolerance rather than the discretisation: {report}"
    );
}
