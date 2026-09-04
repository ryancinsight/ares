use eunomia::{NumericElement, RealField};

/// A linear simplex element: `D + 1` nodes in `D` dimensions.
///
/// A borrowed view over caller-owned node coordinates rather than an owning
/// type. Assembly walks millions of elements, and the node positions already
/// live in the mesh; copying `D + 1` coordinates per element to describe one
/// would be the allocation churn the hot-path rules exist to avoid.
///
/// # Why linear simplices need no quadrature loop
///
/// Shape functions on a linear simplex are affine, so their gradients are
/// **constant** over the element. Every integrand in the stiffness action is
/// therefore constant, and `integral = value * measure` is exact — not a
/// one-point quadrature rule that happens to be accurate enough. Higher-order
/// elements need a quadrature loop; these do not, and pretending otherwise
/// would add a rule whose error term is identically zero.
#[derive(Clone, Copy, Debug)]
pub struct Simplex<'nodes, T, const D: usize> {
    nodes: &'nodes [[T; D]],
}

impl<'nodes, T: RealField, const D: usize> Simplex<'nodes, T, D> {
    /// Borrow `D + 1` node coordinates as an element.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidElement::NodeCount`] when the slice is not exactly
    /// `D + 1` long. A simplex is defined by that count; accepting any other
    /// would silently interpret the extra or missing nodes.
    pub fn try_new(nodes: &'nodes [[T; D]]) -> Result<Self, InvalidElement> {
        if nodes.len() == D + 1 {
            Ok(Self { nodes })
        } else {
            Err(InvalidElement::NodeCount {
                expected: D + 1,
                found: nodes.len(),
            })
        }
    }

    /// The element's node coordinates.
    #[must_use]
    pub const fn nodes(&self) -> &'nodes [[T; D]] {
        self.nodes
    }

    /// Edge matrix `J`, whose columns are `x_i - x_0` for `i` in `1..=D`.
    fn edge_matrix(&self) -> [[T; D]; D] {
        let mut jacobian = [[<T as NumericElement>::ZERO; D]; D];
        // `nodes` has exactly D + 1 entries, checked at construction.
        for (column, node) in self.nodes.iter().skip(1).enumerate() {
            for (row, target) in jacobian.iter_mut().enumerate() {
                target[column] = node[row] - self.nodes[0][row];
            }
        }
        jacobian
    }

    /// Signed measure: area in 2-D, volume in 3-D.
    ///
    /// `det(J) / D!`. The sign follows node ordering, so a negative measure
    /// means the element is inverted rather than merely small.
    #[must_use]
    pub fn signed_measure(&self) -> T {
        let mut factorial = <T as NumericElement>::ONE;
        let mut term = <T as NumericElement>::ONE;
        for _ in 1..=D {
            factorial *= term;
            term += <T as NumericElement>::ONE;
        }
        determinant(&self.edge_matrix()) / factorial
    }

    /// Constant gradients of the `D + 1` shape functions.
    ///
    /// Writes `grad N_i` into `out[i]`. `out` must hold `D + 1` entries.
    ///
    /// # Partition of unity, and how exactly it holds
    ///
    /// The shape functions sum to one everywhere, so their gradients sum to
    /// zero. `grad N_0` is computed as the negated sum of the others rather
    /// than independently, which is the closest this can get to that identity.
    ///
    /// It is **not** identically zero downstream. The cancellation is exact
    /// only when the sum is re-accumulated in the same order; a consumer
    /// summing in another order sees a residual, measured at `1.4e-17` for an
    /// ordinary triangle. Some geometries — axis-aligned ones, where the edge
    /// matrix is diagonal — do cancel exactly, which makes the difference easy
    /// to miss in a test that happens to pick one.
    ///
    /// Consumers that need translation invariance must therefore not rely on
    /// this cancellation.
    /// [`stiffness_action`](crate::stiffness_action) does not: it differences
    /// displacements against a reference node, so a uniform translation gives
    /// an exactly zero gradient whatever these gradients rounded to.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidElement::NodeCount`] when `out` is misshaped, or
    /// [`InvalidElement::Degenerate`] when the element has no measure — a
    /// collapsed element has no well-defined gradient, and returning infinities
    /// would push a plausible-looking `NaN` into the assembled system.
    pub fn shape_gradients(&self, out: &mut [[T; D]]) -> Result<(), InvalidElement> {
        if out.len() != D + 1 {
            return Err(InvalidElement::NodeCount {
                expected: D + 1,
                found: out.len(),
            });
        }

        // grad N_i for i in 1..=D are the rows of J^{-1}; solving J^T G = I
        // gives them directly.
        let inverse = invert(&self.edge_matrix()).ok_or(InvalidElement::Degenerate)?;

        let mut sum = [<T as NumericElement>::ZERO; D];
        for (i, row) in inverse.iter().enumerate() {
            for (component, value) in row.iter().enumerate() {
                out[i + 1][component] = *value;
                sum[component] += *value;
            }
        }
        for (component, total) in sum.iter().enumerate() {
            out[0][component] = -*total;
        }
        Ok(())
    }
}

/// Determinant of a small square matrix by Gaussian elimination.
fn determinant<T: RealField, const D: usize>(matrix: &[[T; D]; D]) -> T {
    let mut work = *matrix;
    let mut result = <T as NumericElement>::ONE;

    for pivot in 0..D {
        let mut best = pivot;
        for row in (pivot + 1)..D {
            if work[row][pivot].abs() > work[best][pivot].abs() {
                best = row;
            }
        }
        if work[best][pivot] == <T as NumericElement>::ZERO {
            return <T as NumericElement>::ZERO;
        }
        if best != pivot {
            work.swap(best, pivot);
            result = -result;
        }
        result *= work[pivot][pivot];
        // The pivot row is copied out so the target row can be iterated
        // rather than indexed; `D` is the spatial dimension, so this is a
        // two- or three-element copy.
        let pivot_row = work[pivot];
        for row in work.iter_mut().skip(pivot + 1) {
            let factor = row[pivot] / pivot_row[pivot];
            for (column, target) in row.iter_mut().enumerate().skip(pivot) {
                *target -= factor * pivot_row[column];
            }
        }
    }
    result
}

/// Inverse of a small square matrix, or `None` when singular.
///
/// Gauss-Jordan with partial pivoting. `D` is the spatial dimension, so this
/// is a two- or three-row solve; a factorisation cache would cost more than it
/// saves.
fn invert<T: RealField, const D: usize>(matrix: &[[T; D]; D]) -> Option<[[T; D]; D]> {
    let mut work = *matrix;
    let mut inverse = [[<T as NumericElement>::ZERO; D]; D];
    for (i, row) in inverse.iter_mut().enumerate() {
        row[i] = <T as NumericElement>::ONE;
    }

    for pivot in 0..D {
        let mut best = pivot;
        for row in (pivot + 1)..D {
            if work[row][pivot].abs() > work[best][pivot].abs() {
                best = row;
            }
        }
        if !work[best][pivot].is_finite() || work[best][pivot] == <T as NumericElement>::ZERO {
            return None;
        }
        work.swap(best, pivot);
        inverse.swap(best, pivot);

        let diagonal = work[pivot][pivot];
        for column in 0..D {
            work[pivot][column] = work[pivot][column] / diagonal;
            inverse[pivot][column] = inverse[pivot][column] / diagonal;
        }
        for row in 0..D {
            if row == pivot {
                continue;
            }
            let factor = work[row][pivot];
            for column in 0..D {
                let a = factor * work[pivot][column];
                let b = factor * inverse[pivot][column];
                work[row][column] -= a;
                inverse[row][column] -= b;
            }
        }
    }
    Some(inverse)
}

/// An element that cannot be interpreted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InvalidElement {
    /// A simplex in `D` dimensions has exactly `D + 1` nodes.
    NodeCount {
        /// Nodes a `D`-simplex requires.
        expected: usize,
        /// Nodes supplied.
        found: usize,
    },
    /// The element has no measure, so its shape gradients are undefined.
    Degenerate,
}

impl core::fmt::Display for InvalidElement {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NodeCount { expected, found } => {
                write!(formatter, "a simplex needs {expected} nodes, got {found}")
            }
            Self::Degenerate => {
                write!(formatter, "the element is degenerate and has no measure")
            }
        }
    }
}

impl core::error::Error for InvalidElement {}
