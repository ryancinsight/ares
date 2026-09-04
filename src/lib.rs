//! Solid momentum balance for the Atlas stack.
//!
//! Ares owns the balance side of solid mechanics: kinematics, stress measures,
//! equilibrium, and the boundary conditions that close them. It owns no
//! material data. Constitutive closure is
//! [`proteus`](https://github.com/ryancinsight/proteus)\'s, per atlas ADR 0055:
//! Proteus closes, Ares balances.
//!
//! # Scope
//!
//! Phase 0 is small-strain linear elastostatics on an unstructured mesh
//! (atlas ADR 0057). Plasticity, finite deformation, contact, dynamics,
//! fracture, and fatigue are later phases and are not scaffolded here — an
//! empty module for a capability that does not exist is a placeholder, and the
//! gate that admitted this package refuses those.
//!
//! # Boundaries
//!
//! Gaia owns mesh, geometry, and proximity queries. Athena owns solver policy.
//! Horae owns time when dynamics arrive. Harmonia owns every coupling to
//! another balance domain, so Ares carries no dependency on a fluid, acoustic,
//! or transport package in any phase.

#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

extern crate alloc;

/// Deformation measures derived from a displacement field.
pub mod kinematics;

pub use kinematics::{AsymmetricInput, SmallStrain, SymmetricTensor};
