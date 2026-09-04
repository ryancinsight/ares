use eunomia::{NumericElement, RealField};

/// A linear simplex element: `D + 1` nodes in `D` dimensions.
///
/// A borrowed view over caller-owned node coordinates rather than an owning
/// type. Assembly walks millions of elements, and the node positions already
/// live in the mesh; copying `D + 1` coordinates per element to describe one
/// would be the allocation churn the hot-path rules exist to avoid.
///
/// # Why the node count is a second const parameter
///
/// `N` is `D + 1`, and stable Rust cannot spell that as an array length. The
/// alternative — a slice checked at construction — makes the count a runtime
/// error that assembly can never actually trigger, since it gathers into a
/// buffer it sized itself. Carrying `N` instead makes the count a compile-time
/// fact: [`Simplex::new`] is infallible, the displacement and force buffers of
/// [`stiffness_action`](super::stiffness_action) cannot be misshaped, and the
/// element's only remaining failure is geometric.
///
/// Both parameters infer from the argument at every call site, so the pair
/// costs no annotation:
///
/// ```
/// use ares::Simplex;
/// // `D = 2` and `N = 3` both come from the argument's type.
/// let triangle = Simplex::new(&[[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]]);
/// assert_eq!(triangle.signed_measure(), 0.5);
/// ```
///
/// A mismatched pair is not a runtime error but a compilation failure — the
/// executable form of the claim that the node count cannot be wrong:
///
/// ```compile_fail
/// use ares::Simplex;
/// // Four nodes in 2-D: `N = 4`, `D + 1 = 3`. The `const` block rejects it.
/// let bad = Simplex::new(&[[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
/// let _ = bad.signed_measure();
/// ```
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
pub struct Simplex<'nodes, T, const D: usize, const N: usize> {
    nodes: &'nodes [[T; D]; N],
}

impl<'nodes, T: RealField, const D: usize, const N: usize> Simplex<'nodes, T, D, N> {
    /// Borrow `D + 1` node coordinates as an element.
    ///
    /// # Panics
    ///
    /// Does not panic at runtime. The `const` block fails compilation when
    /// `N != D + 1`, so an ill-shaped instantiation never links.
    #[must_use]
    pub fn new(nodes: &'nodes [[T; D]; N]) -> Self {
        const { assert!(N == D + 1, "a simplex in D dimensions has D + 1 nodes") }
        Self { nodes }
    }

    /// The element's node coordinates.
    #[must_use]
    pub const fn nodes(&self) -> &'nodes [[T; D]; N] {
        self.nodes
    }

    /// Edge matrix `J`, whose columns are `x_i - x_0` for `i` in `1..=D`.
    ///
    /// Taken by value: a `Simplex` is one reference wide, so borrowing it
    /// costs more than copying it.
    fn edge_matrix(self) -> [[T; D]; D] {
        let mut jacobian = [[<T as NumericElement>::ZERO; D]; D];
        // `N == D + 1`, so skipping node 0 leaves exactly `D` columns.
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

    /// Constant gradients of the `D + 1` shape functions, `grad N_i` at index `i`.
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
    /// [`stiffness_action`](super::stiffness_action) does not: it differences
    /// displacements against a reference node, so a uniform translation gives
    /// an exactly zero gradient whatever these gradients rounded to.
    ///
    /// # Errors
    ///
    /// Returns [`DegenerateElement`] when the element has no measure — a
    /// collapsed element has no well-defined gradient, and returning
    /// infinities would push a plausible-looking `NaN` into the assembled
    /// system, where it is far harder to attribute.
    pub fn shape_gradients(&self) -> Result<[[T; D]; N], DegenerateElement> {
        // grad N_i for i in 1..=D are the rows of J^{-1}; solving J^T G = I
        // gives them directly.
        let inverse = invert(&self.edge_matrix()).ok_or(DegenerateElement)?;

        let mut out = [[<T as NumericElement>::ZERO; D]; N];
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
        Ok(out)
    }
}

/// Determinant of a small square matrix by Gaussian elimination.
fn determinant<T: RealField, const D: usize>(matrix: &[[T; D]; D]) -> T {
    leading_determinant(matrix, D)
}

/// Determinant of the leading `size` by `size` block of a matrix.
///
/// The general form, because a facet's Gram matrix is `(D - 1)` square inside
/// a `D`-square buffer — stable Rust cannot spell `D - 1` as an array length,
/// and a second elimination routine for the smaller block would be the same
/// algorithm twice.
pub(crate) fn leading_determinant<T: RealField, const D: usize>(
    matrix: &[[T; D]; D],
    size: usize,
) -> T {
    let mut work = *matrix;
    let mut result = <T as NumericElement>::ONE;

    for pivot in 0..size {
        let mut best = pivot;
        for row in (pivot + 1)..size {
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
        for row in work.iter_mut().take(size).skip(pivot + 1) {
            let factor = row[pivot] / pivot_row[pivot];
            for (column, target) in row.iter_mut().enumerate().take(size).skip(pivot) {
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

/// The element has no measure, so its shape gradients are undefined.
///
/// The only way an element can fail, now that its node count is structural.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DegenerateElement;

impl core::fmt::Display for DegenerateElement {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "the element is degenerate and has no measure")
    }
}

impl core::error::Error for DegenerateElement {}
