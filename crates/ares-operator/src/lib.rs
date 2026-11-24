#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// The rounding floor of the constrained operator's residual evaluation.
mod floor;
/// The Dirichlet-constrained stiffness as an Athena linear operator.
mod operator;

pub use operator::ConstrainedStiffness;
