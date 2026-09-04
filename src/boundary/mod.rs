//! Boundary conditions that close the momentum balance.
//!
//! Both kinds are typed rather than expressed as index manipulation on an
//! assembled system. A Dirichlet condition applied by striking rows out of a
//! matrix, or a Neumann load applied by adding a force where a traction was
//! meant, produces a system that still solves — so the error surfaces as a
//! wrong displacement field rather than as a failure.

mod dirichlet;
mod traction;

pub use dirichlet::{DirichletConditions, InvalidConditions, PrescribedDisplacement};
pub use traction::{
    InvalidBoundary, MisshapedLoad, TractionBoundary, TractionFacet, facet_measure,
};
