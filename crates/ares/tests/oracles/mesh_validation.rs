//! Mesh and field validation: every condition assembly refuses to run on.

use super::support::{moduli, square_patch};
use ares::{InvalidMesh, MisshapedField, SimplexMesh};

// ---------------------------------------------------------------------------
// Mesh validation
// ---------------------------------------------------------------------------

#[test]
fn an_out_of_range_node_index_is_rejected() {
    let (nodes, _, _) = square_patch();
    let cells = [[0, 1, 9]];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("node 9 is beyond the patch"),
        InvalidMesh::NodeIndexOutOfRange {
            cell: 0,
            position: 2,
            node: 9,
            nodes: 5,
        }
    );
}

#[test]
fn a_degenerate_cell_is_rejected() {
    // Three collinear nodes: no area, so no shape gradients. Left in the mesh
    // it would put infinities into every assembled force.
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0]];
    let cells = [[0, 1, 3], [0, 1, 2]];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("the collinear cell has no area"),
        InvalidMesh::DegenerateCell { cell: 1 }
    );
}

#[test]
fn a_repeated_node_within_a_cell_is_rejected() {
    // Degenerate by a different route: the same node twice collapses the
    // element without any coordinate being unusual.
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let cells = [[0, 1, 1]];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("the repeated node collapses the cell"),
        InvalidMesh::DegenerateCell { cell: 0 }
    );
}

#[test]
fn an_inverted_cell_is_rejected() {
    // A negative measure negates that cell's stiffness, so a mixed-winding
    // mesh assembles an indefinite operator. Rejecting it here is why the
    // operator handed to conjugate gradients is positive definite.
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let cells = [[0, 2, 1]];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("the cell is wound clockwise"),
        InvalidMesh::InvertedCell { cell: 0 }
    );
}

#[test]
fn non_finite_coordinates_are_rejected() {
    // NaN passes both the sign tests, so it needs its own guard: without one
    // it reaches assembly and poisons every node the cell touches.
    let nodes = [[0.0_f64, 0.0], [f64::NAN, 0.0], [0.0, 1.0]];
    let cells = [[0, 1, 2]];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("a NaN coordinate has no measure"),
        InvalidMesh::NonFiniteCell { cell: 0 }
    );
}

#[test]
fn an_empty_mesh_is_rejected() {
    let nodes: [[f64; 2]; 0] = [];
    let cells = [[0, 1, 2]];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("the mesh has no nodes"),
        InvalidMesh::NoNodes
    );

    let (nodes, _, _) = square_patch();
    let cells: [[usize; 3]; 0] = [];
    assert_eq!(
        SimplexMesh::try_new(&nodes, &cells).expect_err("the mesh has no cells"),
        InvalidMesh::NoCells
    );
}

#[test]
fn a_misshaped_field_is_rejected() {
    let (nodes, cells, _) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let m = moduli::<f64>(200e9, 0.3);
    assert_eq!(mesh.degrees_of_freedom(), 10);

    let mut forces = [0.0_f64; 10];
    assert_eq!(
        mesh.internal_forces(&m, &[0.0; 8], &mut forces)
            .expect_err("the displacement field is two entries short"),
        MisshapedField::Displacement {
            expected: 10,
            found: 8
        }
    );

    let mut short = [0.0_f64; 6];
    assert_eq!(
        mesh.internal_forces(&m, &[0.0; 10], &mut short)
            .expect_err("the force field is four entries short"),
        MisshapedField::Force {
            expected: 10,
            found: 6
        }
    );
}
