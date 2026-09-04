//! Neumann conditions: consistent nodal loads from a surface traction.
//!
//! A consistent nodal load is not "the traction split between the nodes" -
//! it is `integral(N_a t) dS`, and the two coincide only for linear elements
//! under a uniform traction. Two independent statics identities pin it: the
//! loads must carry the exact resultant **force** `t * A`, and the exact
//! resultant **moment** about any origin. Force alone does not distinguish an
//! equal split from an unequal one with the same total, so the moment is what
//! fixes the distribution.

use ares::{InvalidBoundary, TractionBoundary, TractionFacet};

// ---------------------------------------------------------------------------
// Neumann: consistent nodal loads
// ---------------------------------------------------------------------------

#[test]
fn a_facet_measure_is_its_length_in_two_dimensions() {
    // A 3-4-5 edge, so the closed form is an exact integer.
    let nodes = [[1.0_f64, 2.0], [4.0, 6.0]];
    assert_eq!(ares::boundary::facet_measure(&nodes), 5.0);
}

#[test]
fn a_facet_measure_is_its_area_in_three_dimensions() {
    // A right triangle with legs 3 and 4 in a plane oblique to every axis, so
    // the Gram determinant is doing real work rather than reading off a
    // coordinate difference.
    let nodes = [[0.0_f64, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 4.0]];
    assert_eq!(ares::boundary::facet_measure(&nodes), 6.0);
}

#[test]
fn a_traction_carries_its_exact_resultant_force() {
    // The first statics identity: sum of nodal loads == t * A.
    let nodes = [[0.0_f64, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]];
    // The right-hand edge, length 3.
    let traction = [1.7e6_f64, -0.4e6];
    let facets = [TractionFacet::new([1, 2], traction)];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");

    let mut loads = [0.0_f64; 8];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    for component in 0..2 {
        let total: f64 = (0..4).map(|node| loads[node * 2 + component]).sum();
        let expected = traction[component] * 3.0;
        let tolerance = expected.abs() * f64::EPSILON * 8.0;
        assert!(
            (total - expected).abs() <= tolerance,
            "component {component}: resultant {total} != {expected}"
        );
    }
}

#[test]
fn a_traction_carries_its_exact_resultant_moment() {
    // The second, independent identity. A uniform traction acts through the
    // facet centroid, so `sum_a x_a (x) f_a` must equal `A * x_centroid (x) t`.
    // The force check above passes for any distribution with the right total;
    // this one fixes the distribution.
    let nodes = [[0.0_f64, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]];
    let traction = [1.7e6_f64, -0.4e6];
    let facets = [TractionFacet::new([1, 2], traction)];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");

    let mut loads = [0.0_f64; 8];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    // Scalar moment about the origin: sum of (x f_y - y f_x).
    let moment: f64 = (0..4)
        .map(|node| nodes[node][0] * loads[node * 2 + 1] - nodes[node][1] * loads[node * 2])
        .sum();
    let centroid = [
        f64::midpoint(nodes[1][0], nodes[2][0]),
        f64::midpoint(nodes[1][1], nodes[2][1]),
    ];
    let expected = 3.0 * (centroid[0] * traction[1] - centroid[1] * traction[0]);
    let tolerance = expected.abs() * f64::EPSILON * 16.0;
    assert!(
        (moment - expected).abs() <= tolerance,
        "resultant moment {moment} != {expected}"
    );
}

#[test]
fn a_three_dimensional_traction_carries_its_exact_resultant() {
    let nodes = [
        [0.0_f64, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 0.0, 4.0],
        [1.0, 5.0, 1.0],
    ];
    let traction = [0.0_f64, 2.5e6, 0.0];
    let facets = [TractionFacet::new([0, 1, 2], traction)];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");

    let mut loads = [0.0_f64; 12];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    let total: f64 = (0..4).map(|node| loads[node * 3 + 1]).sum();
    let expected = traction[1] * 6.0;
    assert!((total - expected).abs() <= expected.abs() * f64::EPSILON * 8.0);
    // Node 3 is not on the facet and must carry nothing.
    for component in 0..3 {
        assert_eq!(loads[3 * 3 + component], 0.0);
    }
}

#[test]
fn a_traction_scales_with_the_facet_rather_than_the_node_count() {
    // The distinction between a traction and a force: doubling the facet
    // doubles the load. A load applied per node instead would not move.
    let short = [[0.0_f64, 0.0], [0.0, 1.0]];
    let long = [[0.0_f64, 0.0], [0.0, 2.0]];
    let traction = [3.0e5_f64, 0.0];
    let facets = [TractionFacet::new([0, 1], traction)];

    let mut short_loads = [0.0_f64; 4];
    TractionBoundary::try_new(&facets, &short)
        .expect("valid")
        .add_consistent_loads(&short, &mut short_loads)
        .expect("well-shaped");
    let mut long_loads = [0.0_f64; 4];
    TractionBoundary::try_new(&facets, &long)
        .expect("valid")
        .add_consistent_loads(&long, &mut long_loads)
        .expect("well-shaped");

    let short_total: f64 = (0..2).map(|node| short_loads[node * 2]).sum();
    let long_total: f64 = (0..2).map(|node| long_loads[node * 2]).sum();
    assert!((long_total - 2.0 * short_total).abs() <= short_total.abs() * f64::EPSILON * 8.0);
}

#[test]
fn tractions_accumulate_rather_than_overwrite() {
    // Two facets sharing a node, plus a pre-existing body force. Assignment
    // instead of accumulation would drop every load but the last.
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [2.0, 0.0]];
    let traction = [0.0_f64, 1.0e5];
    let facets = [
        TractionFacet::new([0, 1], traction),
        TractionFacet::new([1, 2], traction),
    ];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid");

    let mut loads = [0.0_f64, 7.0, 0.0, 0.0, 0.0, 0.0];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    // Node 1 is shared, so it takes half of each facet: a full facet's worth.
    let expected_shared = traction[1] * 1.0;
    assert!((loads[3] - expected_shared).abs() <= expected_shared * f64::EPSILON * 8.0);
    // The pre-existing entry survived.
    assert!((loads[1] - (7.0 + traction[1] / 2.0)).abs() <= traction[1] * f64::EPSILON * 8.0);
}

#[test]
fn an_invalid_traction_boundary_is_rejected() {
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0]];
    let facets = [TractionFacet::new([0, 5], [1.0_f64, 0.0])];
    assert_eq!(
        TractionBoundary::try_new(&facets, &nodes).expect_err("node 5 is beyond the mesh"),
        InvalidBoundary::NodeOutOfRange {
            facet: 0,
            position: 1,
            node: 5,
            nodes: 2,
        }
    );

    let facets = [TractionFacet::new([0, 0], [1.0_f64, 0.0])];
    assert_eq!(
        TractionBoundary::try_new(&facets, &nodes).expect_err("the facet has no length"),
        InvalidBoundary::DegenerateFacet { facet: 0 }
    );
}
