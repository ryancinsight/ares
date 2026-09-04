use eunomia::{NumericElement, RealField};
use proteus::IsotropicModuli;

use super::{DegenerateElement, Simplex, gradient::accumulate_gradient};
use crate::constitutive::isotropic_hooke;
use crate::kinematics::SmallStrain;

/// Apply the element stiffness to nodal displacements: `f = K u`.
///
/// Writes nodal forces into `forces`, one `D`-vector per node.
///
/// # Matrix-free, and why
///
/// Athena's solvers take a `LinearOperator` — they ask for `K u`, never for
/// `K`. So this computes the action directly and never forms an element
/// matrix, which also means no global sparse matrix is assembled or stored.
///
/// It is allocation-free: `N` is a compile-time node count, so the shape
/// gradients live in a stack array. Assembly calls this once per cell, and an
/// allocation there would be a per-element heap round trip in the hottest loop
/// the crate has.
///
/// # Formulation
///
/// Rather than building a `B` matrix and evaluating `B^T D B`, the action goes
/// through the physical quantities that already have owners:
///
/// ```text
/// eps   = sum_a sym(u_a (x) grad N_a)    kinematics
/// sigma = C : eps                        Proteus closure, via isotropic_hooke
/// f_a   = measure * sigma . grad N_a     balance
/// ```
///
/// That is the same operator, expressed in the three quantities the ADR 0055
/// decomposition names. It reuses the strain measure and the constitutive law
/// unchanged, so a defect in either surfaces here rather than being duplicated
/// in a separate `B`-matrix path — and the `B` matrix is where Voigt's shear
/// factor is normally introduced.
///
/// # Theorem: rigid-body motions are in the null space
///
/// A rigid displacement gives an antisymmetric (or zero) gradient, so `eps` is
/// zero, so `sigma` is zero, so every nodal force is. Translation is exact —
/// see the differencing note in the body — while rotation holds to rounding,
/// because its relative displacements do not vanish and the reconstructed
/// gradient is antisymmetric only to rounding.
///
/// # Errors
///
/// Returns [`DegenerateElement`] when the element has collapsed. The buffers
/// cannot be misshaped: `N` fixes their length at compile time.
pub fn stiffness_action<T: RealField, const D: usize, const N: usize>(
    element: &Simplex<'_, T, D, N>,
    moduli: &IsotropicModuli<T>,
    displacements: &[[T; D]; N],
    forces: &mut [[T; D]; N],
) -> Result<(), DegenerateElement> {
    let shape_gradients = element.shape_gradients()?;

    // The differencing that makes a rigid translation exactly force-free
    // lives in `accumulate_gradient`, shared with the strain recovery so the
    // two cannot drift apart.
    let gradient = accumulate_gradient(&shape_gradients, displacements);

    let strain = SmallStrain::<T, D>::from_displacement_gradient(&gradient);
    let stress = isotropic_hooke(moduli, &strain);
    let measure = element.signed_measure();
    let components = *stress.tensor().components();

    // f_a = measure * sigma . grad N_a
    for (force, shape) in forces.iter_mut().zip(shape_gradients.iter()) {
        for (i, entry) in force.iter_mut().enumerate() {
            let mut sum = <T as NumericElement>::ZERO;
            for (j, shape_component) in shape.iter().enumerate() {
                sum += components[i][j] * *shape_component;
            }
            *entry = measure * sum;
        }
    }
    Ok(())
}
