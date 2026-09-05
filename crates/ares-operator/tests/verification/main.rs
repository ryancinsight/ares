//! End-to-end verification of the solid momentum balance (atlas ADR 0057, A6).
//!
//! Every oracle here is closed-form or a conservation statement. No reference
//! implementation exists to difference against, which is the risk ADR 0057
//! records and the reason the oracles are broad rather than deep: each is
//! blind to something the others catch.
//!
//! - **Manufactured solution** is the general one — any smooth field, no
//!   geometric assumption — and its refinement study is what certifies the
//!   element's order.
//! - **Cantilever tip deflection** is an independent structural theory rather
//!   than a restatement of elasticity, and it is approached from the known
//!   direction: linear triangles are stiff in bending, so the computed
//!   deflection sits below the beam value and rises toward it under
//!   refinement. A result that overshot would be evidence of a defect, not
//!   accuracy.
//! - **Energy consistency** is a conservation identity and holds on any mesh
//!   at any resolution, so it separates a solve that is inaccurate from one
//!   that is inconsistent.
//! - **Lame's thick-walled cylinder** is the axisymmetric closed form, and the
//!   only oracle here that exercises a curved boundary.
//!
//! These are integration tests of the whole composition, so they live with the
//! solver rather than with the domain crate: an oracle needs a solve, and the
//! domain crate deliberately has none.

#![expect(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "grid indices and division counts become the reals of a mesh coordinate or a refinement ratio, all far inside f64's exact-integer range. The f64-to-f32 narrowing builds the f32 fixtures from the f64 ones, which is the point of the scalar-generality oracle rather than an accident."
)]

mod cylinder;
mod energy;
mod mesh;
mod mms;
mod scalar;
mod structural;
