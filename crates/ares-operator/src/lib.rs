//! Athena's linear-operator seam for the Ares solid momentum balance.
//!
//! Ares assembles `f = K u` matrix-free over a mesh; Athena's Krylov solvers
//! consume an [`athena_core::LinearOperator`] that answers exactly that question. This crate
//! is the join, and it is deliberately thin — it owns no physics, no
//! discretisation, and no solver policy.
//!
//! # Why this is a separate crate
//!
//! `ares` is `#![no_std]` and depends on nothing but vocabulary crates.
//! Athena's operator trait fixes the error type to its backend's, so an
//! implementation must name a concrete backend, and the only host backend
//! links `std` through `leto`. Implementing the seam inside `ares` would push
//! that dependency into the domain core.
//!
//! The alternative — a cargo feature — is worse: it makes the shipped
//! configuration the one CI does not build by default, and a feature-gated
//! solver path is an untested path. Two crates keep the dependency direction
//! inward and leave both build configurations real.
//!
//! # Why the operator carries no fallible path of its own
//!
//! [`LinearOperator::apply`](athena_core::LinearOperator::apply) must return
//! `B::Error`, which for the host
//! backend is `LetoBackendError` — a closed, non-exhaustive enum with no
//! variant for "this element is degenerate". That is not a limitation to work
//! around but a constraint that shaped the design upstream:
//! `SimplexMesh::try_new` establishes that every cell integrates, and the
//! Dirichlet conditions are validated against the mesh, so by the time an
//! operator exists the only failures left are shape mismatches — which
//! `LengthMismatch` names exactly.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// The Dirichlet-constrained stiffness as an Athena linear operator.
mod operator;

pub use operator::ConstrainedStiffness;
