//! Global assembly over a simplex mesh.
//!
//! The element kernel gives the stiffness action on one cell; this walks the
//! mesh, gathers each cell's displacements, and scatters its nodal forces back.
//! It is matrix-free throughout: no global sparse matrix is formed, because
//! nothing downstream reads one.

mod action;
mod body;
mod constrained;
mod mesh;

pub use action::MisshapedField;
pub use mesh::{InvalidMesh, SimplexMesh};
