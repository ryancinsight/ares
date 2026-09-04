//! Fixtures and helpers shared by the oracle modules.
//!
//! One definition each. These are the geometries every oracle is measured on,
//! so a copy per module would be three chances for them to drift apart while
//! every test still passed.

use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use eunomia::RealField;
use proteus::IsotropicModuli;

pub fn moduli<T: RealField>(young: f64, poisson: f64) -> IsotropicModuli<T> {
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(T::from_f64(young)),
        Dimensionless::from_base(T::from_f64(poisson)),
    )
    .expect("inside the positive-definite domain")
}

/// Lame parameters, for deriving bounds that scale with stiffness.
pub fn lame(young: f64, poisson: f64) -> (f64, f64) {
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    (lambda, mu)
}

/// The unit triangle: (0,0), (1,0), (0,1).
pub fn unit_triangle() -> [[f64; 2]; 3] {
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
}

/// The unit tetrahedron: origin plus the three axis points.
pub fn unit_tetrahedron() -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

/// A distorted unit square: four corners and one off-centre interior node,
/// triangulated into four cells that meet at the interior node.
///
/// The interior node sits at `(0.37, 0.61)` rather than the centre so that no
/// cell is a mirror of another. A symmetric patch lets sign errors cancel in
/// pairs, which is exactly the failure the patch test exists to catch.
pub fn square_patch() -> ([[f64; 2]; 5], [[usize; 3]; 4], usize) {
    let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.37, 0.61]];
    let cells = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
    (nodes, cells, 4)
}

/// A distorted hexahedral patch: eight perturbed corners and one interior
/// node, coned into twelve tetrahedra over a triangulated boundary.
///
/// The corners are perturbed off the unit cube so the cells are genuinely
/// irregular; an undistorted cube has cells related by symmetry.
pub fn cube_patch() -> ([[f64; 3]; 9], [[usize; 4]; 12], usize) {
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.07, 0.0, 0.0],
        [1.0, 0.93, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.11],
        [1.0, 0.0, 0.96],
        [1.04, 1.0, 1.0],
        [0.0, 1.09, 1.0],
        [0.46, 0.53, 0.58],
    ];
    // The six faces, each split into two triangles wound so that the cone to
    // the interior node 8 has positive volume. Verified by construction: the
    // mesh constructor rejects an inverted cell, so a wrong winding here is a
    // test failure rather than a silently wrong oracle.
    let cells = [
        // Each face is split into two triangles wound so the cone to the
        // interior node has positive volume â€” that is, wound clockwise seen
        // from outside, so the winding normal points inward at the apex.
        // z = 0
        [0, 1, 2, 8],
        [0, 2, 3, 8],
        // z = 1
        [4, 7, 6, 8],
        [4, 6, 5, 8],
        // y = 0
        [0, 4, 5, 8],
        [0, 5, 1, 8],
        // y = 1
        [3, 2, 6, 8],
        [3, 6, 7, 8],
        // x = 0
        [0, 3, 7, 8],
        [0, 7, 4, 8],
        // x = 1
        [1, 5, 6, 8],
        [1, 6, 2, 8],
    ];
    (nodes, cells, 8)
}
