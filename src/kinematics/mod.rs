//! Kinematics: deformation measures derived from a displacement field.
//!
//! Ares owns kinematics because they derive from the displacement it solves
//! for (atlas ADR 0055 R5: kinematics belong to the owner of the primal
//! field). Nothing here consumes a material property; strain is geometry.

mod strain;
mod tensor;

pub use strain::SmallStrain;
pub use tensor::{AsymmetricInput, SymmetricTensor};
