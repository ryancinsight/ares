use ares::{SimplexMesh, TractionBoundary, TractionFacet};
use ares_athena::ConstrainedStiffness;
use athena_core::{Cg, CgWorkspace, ConvergencePolicy, Identity, KrylovBackend, Termination};
use athena_leto::LetoBackend;
use eunomia::{NumericElement, RealField};
use harmonia::{Partition, Substep};
use leto::Array1;
use leto_ops::RealScalar;

use crate::interface::{InvalidInterface, StructuralInterface};

/// A structural solve, presented to Harmonia as a coupling partition.
///
/// `state` is the full nodal displacement field, `input` an interface traction
/// exchange, and `output` an interface displacement exchange. Both exchange
/// orderings are fixed by [`StructuralInterface`].
///
/// # Quasi-static, and what that means for the substep
///
/// `advance` solves equilibrium for the traction it was handed and does not
/// integrate in time. Phase 0 of Ares is static (atlas ADR 0057), so there is
/// no time derivative and no velocity continuity to impose; Phase 0 of the
/// coupling is one-way (atlas ADR 0059), so the exported displacement does not
/// move the fluid mesh.
///
/// The substep is therefore unused rather than partly used. It is not
/// validated either, and deliberately: `StepSize::new` already refuses a
/// non-positive step, so a check here would be unreachable code asserting an
/// invariant its own type system enforces upstream.
//
// The bound sits on the definition rather than only on the impls, which the
// bound-placement rule otherwise forbids. It is the sanctioned exception: the
// `CgWorkspace<LetoBackend<T>>` field cannot be named without it.
pub struct StructuralPartition<'system, T: RealScalar + RealField, const D: usize, const N: usize> {
    operator: ConstrainedStiffness<'system, T, D, N>,
    mesh: SimplexMesh<'system, T, D, N>,
    interface: StructuralInterface<'system, D>,
    facets: Vec<TractionFacet<T, D>>,
    external: Vec<T>,
    load: Vec<T>,
    policy: ConvergencePolicy<T>,
    backend: LetoBackend<T>,
    workspace: CgWorkspace<LetoBackend<T>>,
    solution: Array1<T>,
}

impl<'system, T: RealScalar + RealField, const D: usize, const N: usize>
    StructuralPartition<'system, T, D, N>
{
    /// Build the partition for a mesh, its constrained operator, and an
    /// interface.
    ///
    /// Allocates every buffer the coupling loop needs once, including the
    /// Krylov workspace. Harmonia drives `advance` repeatedly, so allocating
    /// per step would be a heap round trip inside the coupling iteration.
    ///
    /// # Errors
    ///
    /// Returns [`PartitionError::Solve`] when the Krylov workspace cannot be
    /// allocated.
    pub fn try_new(
        mesh: SimplexMesh<'system, T, D, N>,
        operator: ConstrainedStiffness<'system, T, D, N>,
        interface: StructuralInterface<'system, D>,
        policy: ConvergencePolicy<T>,
    ) -> Result<Self, PartitionError> {
        let dofs = mesh.degrees_of_freedom();
        let backend = LetoBackend::<T>::default();
        let workspace = CgWorkspace::new(&backend, dofs).map_err(|_| PartitionError::Solve)?;
        Ok(Self {
            operator,
            mesh,
            interface,
            facets: vec![
                TractionFacet::new([0; D], [<T as NumericElement>::ZERO; D]);
                interface.facets().len()
            ],
            external: vec![<T as NumericElement>::ZERO; dofs],
            load: vec![<T as NumericElement>::ZERO; dofs],
            policy,
            backend,
            workspace,
            solution: Array1::zeros([dofs]),
        })
    }

    /// The interface this partition exchanges over.
    #[must_use]
    pub const fn interface(&self) -> &StructuralInterface<'system, D> {
        &self.interface
    }

    /// The mesh it solves on.
    #[must_use]
    pub const fn mesh(&self) -> &SimplexMesh<'system, T, D, N> {
        &self.mesh
    }

    /// The work `integral(t . u) dS` the interface traction does on a state.
    ///
    /// # Errors
    ///
    /// Returns [`PartitionError::Interface`] on a length mismatch.
    pub fn interface_work(&self, traction: &[T], state: &[T]) -> Result<T, PartitionError> {
        self.interface
            .interface_work(&self.mesh, traction, state)
            .map_err(PartitionError::Interface)
    }

    /// The strain energy `(1/2) u . K u` stored by a displacement state.
    ///
    /// # Errors
    ///
    /// Returns [`PartitionError::Solve`] when the state is misshaped.
    pub fn strain_energy(&self, state: &[T]) -> Result<T, PartitionError> {
        let mut forces = vec![<T as NumericElement>::ZERO; state.len()];
        self.mesh
            .internal_forces(self.operator.moduli(), state, &mut forces)
            .map_err(|_| PartitionError::Solve)?;
        let mut total = <T as NumericElement>::ZERO;
        for (displacement, force) in state.iter().zip(forces.iter()) {
            total += *displacement * *force;
        }
        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        Ok(total / two)
    }

    /// Solve equilibrium for one interface traction, writing the nodal
    /// displacement into `state`.
    ///
    /// This is the whole of the coupling step; [`Partition::advance`] is a
    /// delegation to it that discards the substep.
    ///
    /// It exists separately because `harmonia::Substep` has no public
    /// constructor, so `advance` can only be reached through Harmonia's own
    /// coupling driver. Routing the work through an inherent method keeps the
    /// physics reachable from a test that is about the physics, rather than
    /// forcing every such test to stand up a two-partition driver first.
    ///
    /// # Errors
    ///
    /// Returns [`PartitionError`] when the exchange is misshaped, the traction
    /// is not a valid boundary, or the solve fails to converge.
    pub fn solve_for_traction(
        &mut self,
        state: &mut [T],
        input: &[T],
    ) -> Result<(), PartitionError> {
        if state.len() != self.state_dimension() {
            return Err(PartitionError::Interface(
                InvalidInterface::MisshapedExchange {
                    expected: self.state_dimension(),
                    found: state.len(),
                },
            ));
        }

        self.interface
            .read_traction(input, &mut self.facets)
            .map_err(PartitionError::Interface)?;
        let boundary = TractionBoundary::try_new(&self.facets, self.mesh.nodes())
            .map_err(|_| PartitionError::Boundary)?;

        for slot in &mut self.external {
            *slot = <T as NumericElement>::ZERO;
        }
        boundary
            .add_consistent_loads(self.mesh.nodes(), &mut self.external)
            .map_err(|_| PartitionError::Boundary)?;
        self.operator
            .load(&self.external, &mut self.load)
            .map_err(|_| PartitionError::Solve)?;

        let right_hand_side = Array1::from_shape_vec([self.load.len()], self.load.clone())
            .map_err(|_| PartitionError::Solve)?;
        // Warm start from the incoming state: successive coupling iterations
        // differ by whatever the traction changed, so the previous solution is
        // a better guess than zero and the Krylov count falls accordingly.
        {
            let mut view = self.backend.view_mut(&mut self.solution);
            let slots = view.as_mut_slice().ok_or(PartitionError::Solve)?;
            slots.copy_from_slice(state);
        }

        let report = Cg::<LetoBackend<T>>::solve_into(
            &self.backend,
            &self.operator,
            &Identity,
            &right_hand_side,
            &mut self.solution,
            &mut self.workspace,
            self.policy,
        )
        .map_err(|_| PartitionError::Solve)?;
        if !report.converged() {
            return Err(PartitionError::NotConverged(report.termination));
        }

        let solved = self.backend.view(&self.solution);
        state.copy_from_slice(solved.as_slice().ok_or(PartitionError::Solve)?);
        Ok(())
    }
}

impl<T: RealScalar + RealField, const D: usize, const N: usize> Partition<T>
    for StructuralPartition<'_, T, D, N>
{
    type Error = PartitionError;

    fn state_dimension(&self) -> usize {
        self.mesh.degrees_of_freedom()
    }

    fn input_dimension(&self) -> usize {
        self.interface.input_dimension()
    }

    fn output_dimension(&self) -> usize {
        self.interface.output_dimension()
    }

    fn advance(
        &mut self,
        _substep: Substep<T>,
        state: &mut [T],
        input: &[T],
    ) -> Result<(), Self::Error> {
        self.solve_for_traction(state, input)
    }

    fn export(&self, state: &[T], output: &mut [T]) -> Result<(), Self::Error> {
        self.interface
            .write_displacement(state, output)
            .map_err(PartitionError::Interface)
    }
}

/// A failure of the structural partition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum PartitionError {
    /// The exchange does not match the interface.
    Interface(InvalidInterface),
    /// The interface traction is not a valid traction boundary.
    Boundary,
    /// The structural system could not be assembled or its workspace built.
    Solve,
    /// The Krylov solve did not converge.
    ///
    /// Surfaced rather than absorbed: an unconverged displacement field is
    /// smooth and plausible, so a coupling driver handed one would iterate
    /// against a wrong structural response without any sign of it.
    NotConverged(Termination),
}

impl core::fmt::Display for PartitionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Interface(inner) => write!(formatter, "{inner}"),
            Self::Boundary => write!(formatter, "the interface traction is not a valid boundary"),
            Self::Solve => write!(formatter, "the structural system could not be assembled"),
            Self::NotConverged(termination) => write!(
                formatter,
                "the structural solve did not converge: {termination:?}"
            ),
        }
    }
}

impl core::error::Error for PartitionError {}
