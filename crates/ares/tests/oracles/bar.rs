//! A one-dimensional mesh, the smallest spatial dimension assembly admits
//! (ADR 0004).
//!
//! A zero-dimensional mesh does not compile, and the `compile_fail` example on
//! `SimplexMesh` pins that side. These oracles pin the adjacent valid side: a
//! unit bar `[0, 1]` as one two-node cell, whose fields are read as one-wide
//! chunks. Every expected value below is derived from the bar's geometry and
//! stays exact at `f32` and `f64`, since each is a small dyadic rational.

use super::support::moduli;
use ares::SimplexMesh;
use eunomia::RealField;

fn unit_bar<T: RealField>() -> ([[T; 1]; 2], [[usize; 2]; 1]) {
    ([[T::from_f64(0.0)], [T::from_f64(1.0)]], [[0, 1]])
}

fn assert_a_unit_bar_assembles_its_exact_body_load<T: RealField>() {
    let (nodes, cells) = unit_bar::<T>();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("invariant: the bar has length 1");
    assert_eq!(mesh.degrees_of_freedom(), 2);

    // The consistent load is `f_a = |Om| (sum_b b_b + b_a) / ((D + 1)(D + 2))`.
    // With `|Om| = 1`, `D = 1` and a body force of 2 at both nodes it is
    // `(4 + 2) / 6 = 1` at each node.
    let two = T::from_f64(2.0);
    let mut loads = [T::from_f64(0.0); 2];
    mesh.add_body_force(&[two, two], &mut loads)
        .expect("invariant: two entries hold two 1-D nodes");
    assert_eq!(loads, [T::from_f64(1.0); 2]);
}

#[test]
fn a_unit_bar_assembles_its_exact_body_load() {
    assert_a_unit_bar_assembles_its_exact_body_load::<f32>();
    assert_a_unit_bar_assembles_its_exact_body_load::<f64>();
}

fn assert_a_stretched_unit_bar_recovers_its_exact_strain<T: RealField>() {
    let (nodes, cells) = unit_bar::<T>();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("invariant: the bar has length 1");

    // Node 1 displaced by 1/2 over a bar of length 1 is the strain 1/2.
    let displacements = [T::from_f64(0.0), T::from_f64(0.5)];
    let strains: Vec<_> = mesh
        .cell_strains(&displacements)
        .expect("invariant: two entries hold two 1-D nodes")
        .collect();
    assert_eq!(strains.len(), 1);
    assert_eq!(*strains[0].tensor().components(), [[T::from_f64(0.5)]]);
}

#[test]
fn a_stretched_unit_bar_recovers_its_exact_strain() {
    assert_a_stretched_unit_bar_recovers_its_exact_strain::<f32>();
    assert_a_stretched_unit_bar_recovers_its_exact_strain::<f64>();
}

fn assert_a_stretched_unit_bar_carries_equal_and_opposite_end_forces<T: RealField>() {
    let (nodes, cells) = unit_bar::<T>();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("invariant: the bar has length 1");

    // With Poisson's ratio 0 the Lame parameters are `lambda = 0` and
    // `mu = E / 2`, so the one-dimensional stress is `(lambda + 2 mu) eps =
    // E eps = 4 * 1/2 = 2`. The shape gradients of the two nodes are -1 and
    // +1, so the end forces are `|Om| sigma grad N_a = -2` and `+2`.
    let displacements = [T::from_f64(0.0), T::from_f64(0.5)];
    let mut forces = [T::from_f64(7.0); 2];
    mesh.internal_forces(&moduli::<T>(4.0, 0.0), &displacements, &mut forces)
        .expect("invariant: two entries hold two 1-D nodes");
    assert_eq!(forces, [T::from_f64(-2.0), T::from_f64(2.0)]);
}

#[test]
fn a_stretched_unit_bar_carries_equal_and_opposite_end_forces() {
    assert_a_stretched_unit_bar_carries_equal_and_opposite_end_forces::<f32>();
    assert_a_stretched_unit_bar_carries_equal_and_opposite_end_forces::<f64>();
}
