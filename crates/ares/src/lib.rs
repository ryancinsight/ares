#![doc = include_str!("../README.md")]
#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// Global assembly of the stiffness action over a mesh.
pub mod assembly;
/// Dirichlet and Neumann conditions that close the balance.
pub mod boundary;
/// Constitutive coupling between strain and stress.
pub mod constitutive;
/// Linear simplex elements and their stiffness action.
pub mod element;
/// Deformation measures derived from a displacement field.
pub mod kinematics;

pub use assembly::{InvalidMesh, MisshapedField, SimplexMesh};
pub use boundary::{
    DirichletConditions, InvalidBoundary, InvalidConditions, MisshapedLoad, PrescribedDisplacement,
    TractionBoundary, TractionFacet,
};
pub use constitutive::{CauchyStress, PressureSquared, isotropic_hooke};
pub use element::{DegenerateElement, Simplex, stiffness_action};
pub use kinematics::{AsymmetricInput, SmallStrain, SymmetricTensor};
