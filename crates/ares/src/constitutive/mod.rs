//! Constitutive coupling: strain to stress through a Proteus closure.
//!
//! Ares balances and Proteus closes (atlas ADR 0055). This module is the seam
//! between them: it applies a material law it is handed and stores no material
//! data of its own. No alloy is named anywhere in this crate.

mod hooke;
mod stress;

pub use hooke::isotropic_hooke;
pub use stress::CauchyStress;
