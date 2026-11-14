//! Analytical oracles for the Phase 0 solid momentum balance (atlas ADR 0057).
//!
//! One test binary rather than one per area. Each integration-test file is an
//! independent full link that re-instantiates and re-codegens the whole crate,
//! so a file per concern buys separation at the cost of multiplying build and
//! link work; modules give the same separation for one binary.
//!
//! There is no reference implementation to difference against, which ADR 0057
//! records as this package's central risk. Oracle breadth is the mitigation:
//! each module below is blind to something another catches, and the load
//! -bearing ones are exact rather than approximate.

#![expect(
    clippy::float_cmp,
    reason = "the exact comparisons are exact by construction rather than by luck: unit-simplex measures are ratios of small integers; a rigid translation differences to identically zero relative displacements, so the forces it produces are sums of exact zeros; and a projected field's constrained entries are assigned or copied rather than computed. Everything with genuine rounding - the patch-test residual, the recovered strain, the statics identities, equilibrium, and rotation - carries a derived bound instead."
)]
#![expect(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "the f32 fixtures are built by narrowing the f64 ones, which is the point of the scalar-generality oracle rather than an accident: a bound computed in f64 and applied to an f32 result would not be an f32 bound. The remaining casts turn small loop and node counts into the reals of a derived bound, far inside f64's exact-integer range."
)]

mod assembled_operator;
mod constitutive;
mod dirichlet;
mod element_geometry;
mod element_stiffness;
mod kinematics;
mod mesh_validation;
mod neumann;
mod patch_test;
mod support;
