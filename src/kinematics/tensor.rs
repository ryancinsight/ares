use aequitas::systems::si::quantities::Dimensionless;
use eunomia::{NumericElement, RealField};

/// A symmetric second-order tensor in `D` dimensions.
///
/// # Full storage rather than Voigt
///
/// The `D(D+1)/2` independent components would fit a Voigt vector, and this
/// stores the full `D x D` array instead — nine scalars for a 3-D tensor
/// rather than six.
///
/// The three saved scalars are not worth the bug class they buy. Voigt
/// notation carries a factor of two on its shear entries that differs between
/// strain and stress, and every mixed operation has to remember which
/// convention each operand is in. That factor is the classic silent error in
/// elasticity code: it produces a stiffness that is wrong only in shear, which
/// a uniaxial test never exercises and a converged solve never reveals.
///
/// Full storage makes the tensor algebra literal. `D` is at most three here,
/// so the storage difference is three scalars per tensor.
///
/// # Symmetry
///
/// Symmetry is a construction invariant, not a checked property: the
/// constructors either symmetrise their input or build from components that
/// cannot be asymmetric. There is no way to write an off-diagonal pair
/// independently.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymmetricTensor<T, const D: usize> {
    components: [[T; D]; D],
}

impl<T: RealField, const D: usize> SymmetricTensor<T, D> {
    /// The zero tensor.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            components: [[<T as NumericElement>::ZERO; D]; D],
        }
    }

    /// Symmetrise a general second-order tensor: `(A + A^T) / 2`.
    ///
    /// # Theorem
    ///
    /// For antisymmetric `A`, where `A[j][i] == -A[i][j]`, every entry of the
    /// result is `(a + (-a)) / 2`. In IEEE-754 arithmetic `a + (-a)` is
    /// exactly `+0.0` for every finite `a`, so the result is **exactly** zero,
    /// not merely small.
    ///
    /// That exactness is what makes an infinitesimal rigid rotation produce
    /// exactly zero strain, which is asserted rather than approximated in the
    /// rigid-body test.
    #[must_use]
    pub fn from_symmetrised(general: &[[T; D]; D]) -> Self {
        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        let mut components = [[<T as NumericElement>::ZERO; D]; D];
        for (i, row) in components.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                // Both operands are read from the input; `general` is indexed
                // only at `i`/`j` already proven in range by the iterators.
                let upper = general[i][j];
                let lower = general[j][i];
                *entry = (upper + lower) / two;
            }
        }
        Self { components }
    }

    /// Build directly from an array that is already symmetric.
    ///
    /// # Errors
    ///
    /// Returns [`AsymmetricInput`] naming the first `(i, j)` whose transpose
    /// entry differs. Rejecting rather than silently symmetrising keeps a
    /// caller from believing an asymmetric tensor round-trips.
    pub fn try_from_components(components: [[T; D]; D]) -> Result<Self, AsymmetricInput> {
        for (i, row) in components.iter().enumerate() {
            for (j, entry) in row.iter().enumerate() {
                if *entry != components[j][i] {
                    return Err(AsymmetricInput { row: i, column: j });
                }
            }
        }
        Ok(Self { components })
    }

    /// The spatial dimension as a scalar, built by repeated addition.
    ///
    /// `D as f64` would be a lossy-cast lint and a needless one: summing `ONE`
    /// `D` times is exact for every dimension a tensor is written in, and does
    /// not depend on the target's `usize` width.
    #[must_use]
    fn dimension() -> T {
        let mut count = <T as NumericElement>::ZERO;
        for _ in 0..D {
            count += <T as NumericElement>::ONE;
        }
        count
    }

    /// Trace, the first invariant `tr(A) = sum A[i][i]`.
    #[must_use]
    pub fn trace(&self) -> T {
        let mut sum = <T as NumericElement>::ZERO;
        for (i, row) in self.components.iter().enumerate() {
            sum += row[i];
        }
        sum
    }

    /// Deviatoric part `A - tr(A)/D * I`, which is trace-free.
    #[must_use]
    pub fn deviator(&self) -> Self {
        let mut components = self.components;
        let mean = self.trace() / Self::dimension();
        for (i, row) in components.iter_mut().enumerate() {
            row[i] -= mean;
        }
        Self { components }
    }

    /// Double contraction `A : B = sum_ij A[i][j] B[i][j]`.
    ///
    /// The energy pairing between a stress and a strain.
    #[must_use]
    pub fn double_dot(&self, other: &Self) -> T {
        let mut sum = <T as NumericElement>::ZERO;
        for (row, other_row) in self.components.iter().zip(other.components.iter()) {
            for (entry, other_entry) in row.iter().zip(other_row.iter()) {
                sum += *entry * *other_entry;
            }
        }
        sum
    }

    /// Borrow the full component array.
    #[must_use]
    pub const fn components(&self) -> &[[T; D]; D] {
        &self.components
    }
}

impl<T: RealField, const D: usize> SymmetricTensor<T, D> {
    /// One component, as a dimensionless quantity.
    ///
    /// Returns `None` when either index is outside the tensor.
    #[must_use]
    pub fn component(&self, row: usize, column: usize) -> Option<Dimensionless<T>> {
        self.components
            .get(row)
            .and_then(|entries| entries.get(column))
            .map(|value| Dimensionless::from_base(*value))
    }
}

/// A component array offered to [`SymmetricTensor::try_from_components`] that
/// is not symmetric.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AsymmetricInput {
    /// Row of the first entry differing from its transpose.
    pub row: usize,
    /// Column of that entry.
    pub column: usize,
}

impl core::fmt::Display for AsymmetricInput {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "component ({}, {}) differs from its transpose",
            self.row, self.column
        )
    }
}

impl core::error::Error for AsymmetricInput {}
