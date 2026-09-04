//! Linear simplex elements and their stiffness action.
//!
//! Elements are borrowed views over mesh-owned node coordinates. Gaia owns the
//! mesh (atlas ADR 0055); this module interprets a cell's nodes as an
//! element and never stores geometry of its own.

mod gradient;
mod simplex;
mod stiffness;

pub use simplex::{DegenerateElement, Simplex};
pub use stiffness::stiffness_action;
