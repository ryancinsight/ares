//! Athena's linear-operator seam for the Ares solid momentum balance.
//!
//! Ares assembles `f = K u` matrix-free over a mesh; Athena's Krylov solvers
//! consume a [`LinearOperator`] that answers exactly that question. This crate
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
//! [`LinearOperator::apply`] must return `B::Error`, which for the host
//! backend is `LetoBackendError` — a closed, non-exhaustive enum with no
//! variant for "this element is degenerate". That is not a limitation to work
//! around but a constraint that shaped the design upstream:
//! `SimplexMesh::try_new` establishes that every cell integrates, and the
//! Dirichlet conditions are validated against the mesh, so by the time an
//! operator exists the only failures left are shape mismatches — which
//! `LengthMismatch` names exactly.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

use ares::{DirichletConditions, SimplexMesh};
use athena_core::LinearOperator;
use athena_leto::{LetoBackend, LetoBackendError};
use eunomia::RealField;
use leto::{ArrayView1, ArrayViewMut1};
use leto_ops::RealScalar;
use proteus::IsotropicModuli;

/// The Dirichlet-constrained stiffness of one mesh, as an Athena operator.
///
/// Borrows the mesh, conditions, and material; owns only the scratch buffer
/// the constrained action needs. A Krylov solve applies this once per
/// iteration, so the buffer is allocated at construction and reused — the
/// domain crate refuses to allocate, and this is the layer that is allowed to.
///
/// # Interior mutability, and why it is here rather than in `ares`
///
/// `LinearOperator::apply` takes `&self`, because a solver holds the operator
/// immutably while iterating. The scratch buffer must nonetheless be written
/// on every application, so it lives behind a [`RefCell`](core::cell::RefCell)
/// — the one place in this stack where that is the right answer, precisely
/// because it is confined to a buffer no caller can observe. `ares` keeps the
/// scratch as an explicit parameter, so the domain crate stays free of
/// interior mutability and a caller with a `&mut` buffer never pays for a
/// borrow flag.
pub struct ConstrainedStiffness<'system, T, const D: usize, const N: usize> {
    mesh: SimplexMesh<'system, T, D, N>,
    moduli: IsotropicModuli<T>,
    conditions: DirichletConditions<'system, T, D>,
    scratch: core::cell::RefCell<Vec<T>>,
}

impl<'system, T: RealField, const D: usize, const N: usize> ConstrainedStiffness<'system, T, D, N> {
    /// Build the operator for a mesh, material, and condition set.
    ///
    /// Allocates the scratch buffer once. The mesh and conditions were already
    /// validated against each other by their own constructors, so this cannot
    /// fail.
    #[must_use]
    pub fn new(
        mesh: SimplexMesh<'system, T, D, N>,
        moduli: IsotropicModuli<T>,
        conditions: DirichletConditions<'system, T, D>,
    ) -> Self {
        let scratch = vec![<T as eunomia::NumericElement>::ZERO; mesh.degrees_of_freedom()];
        Self {
            mesh,
            moduli,
            conditions,
            scratch: core::cell::RefCell::new(scratch),
        }
    }

    /// The mesh this operator acts on.
    #[must_use]
    pub const fn mesh(&self) -> &SimplexMesh<'system, T, D, N> {
        &self.mesh
    }

    /// The conditions constraining it.
    #[must_use]
    pub const fn conditions(&self) -> &DirichletConditions<'system, T, D> {
        &self.conditions
    }

    /// The material closing it.
    #[must_use]
    pub const fn moduli(&self) -> &IsotropicModuli<T> {
        &self.moduli
    }

    /// Build the right-hand side for an applied load.
    ///
    /// Wraps [`SimplexMesh::constrained_load`] with this operator's scratch, so
    /// the load and the operator cannot disagree about the conditions they
    /// were built from.
    ///
    /// # Errors
    ///
    /// Returns [`LetoBackendError::LengthMismatch`] when the load lengths do
    /// not match the mesh.
    pub fn load(&self, external: &[T], right_hand_side: &mut [T]) -> Result<(), LetoBackendError> {
        let mut scratch = self.scratch.borrow_mut();
        self.mesh
            .constrained_load(
                &self.moduli,
                &self.conditions,
                external,
                right_hand_side,
                &mut scratch,
            )
            .map_err(|_| LetoBackendError::LengthMismatch {
                left: self.mesh.degrees_of_freedom(),
                right: external.len().min(right_hand_side.len()),
            })
    }
}

impl<T: RealScalar + RealField, const D: usize, const N: usize> LinearOperator<LetoBackend<T>>
    for ConstrainedStiffness<'_, T, D, N>
{
    fn dimension(&self) -> usize {
        self.mesh.degrees_of_freedom()
    }

    fn apply(
        &self,
        _backend: &LetoBackend<T>,
        input: ArrayView1<'_, T>,
        mut output: ArrayViewMut1<'_, T>,
    ) -> Result<(), LetoBackendError> {
        // Athena's host vectors are contiguous by construction, but a strided
        // view is representable, and assembly indexes by degree of freedom. A
        // silent stride would read the wrong entries rather than fail.
        let input = input
            .as_slice()
            .ok_or(LetoBackendError::NonContiguousVector)?;
        let output = output
            .as_mut_slice()
            .ok_or(LetoBackendError::NonContiguousVector)?;
        let mut scratch = self.scratch.borrow_mut();

        self.mesh
            .constrained_action(&self.moduli, &self.conditions, input, output, &mut scratch)
            .map_err(|_| LetoBackendError::LengthMismatch {
                left: self.mesh.degrees_of_freedom(),
                right: input.len(),
            })
    }
}
