use eunomia::{NumericElement, RealField};

use crate::element::Simplex;

/// A conforming simplex mesh: node coordinates plus cell connectivity.
///
/// A borrowed view, like [`Simplex`]. Gaia owns mesh generation, geometry, and
/// proximity queries (atlas ADR 0055); this is the shape assembly reads, not a
/// mesh representation competing with Gaia's. Any producer that can lend node
/// coordinates and connectivity satisfies it.
///
/// # Why coordinates are nodal and fields are flat
///
/// Node coordinates arrive as `&[[T; D]]` and displacement or force fields as
/// flat `&[T]`, and the asymmetry is deliberate. Geometry is built once and
/// never leaves the crate. A field crosses the solver boundary on **every**
/// Krylov iteration, and Athena's vector views are flat, so a nodal field type
/// would impose a copy per iteration. Flat storage makes that boundary
/// zero-copy, and `as_chunks` recovers the nodal view inside the loop for
/// free.
///
/// # Validation, and what it buys
///
/// [`SimplexMesh::try_new`] rejects every cell that assembly could not
/// integrate: an out-of-range node index, a cell with no measure, a cell wound
/// the wrong way, and a cell whose coordinates are not finite. It establishes
/// those by running the same [`Simplex::shape_gradients`] call assembly will
/// run, so the invariant is not merely implied — it is the recorded outcome of
/// the identical computation on the identical data.
///
/// That is what makes
/// [`internal_forces`](SimplexMesh::internal_forces) fail only on a misshaped
/// field, which in turn is what lets the Athena operator report the errors its
/// backend defines rather than needing one of its own.
///
/// Inversion is rejected because it is a correctness failure, not a
/// convention. A negative measure negates that cell's stiffness contribution,
/// so a mesh with mixed winding assembles an indefinite operator — one that
/// stores negative energy under deformation, and on which conjugate gradients
/// has no reason to converge.
#[derive(Clone, Copy, Debug)]
pub struct SimplexMesh<'mesh, T, const D: usize, const N: usize> {
    nodes: &'mesh [[T; D]],
    cells: &'mesh [[usize; N]],
}

impl<'mesh, T: RealField, const D: usize, const N: usize> SimplexMesh<'mesh, T, D, N> {
    /// Borrow node coordinates and connectivity as a mesh.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidMesh`] naming the offending cell. Every variant is a
    /// condition under which assembly would produce a plausible wrong answer
    /// rather than an obvious one, which is why they are rejected here instead
    /// of being tolerated in the loop.
    ///
    /// # Panics
    ///
    /// Does not panic at runtime. The `const` block fails compilation when
    /// `N != D + 1`.
    pub fn try_new(
        nodes: &'mesh [[T; D]],
        cells: &'mesh [[usize; N]],
    ) -> Result<Self, InvalidMesh> {
        const { assert!(N == D + 1, "a simplex in D dimensions has D + 1 nodes") }

        if nodes.is_empty() {
            return Err(InvalidMesh::NoNodes);
        }
        if cells.is_empty() {
            return Err(InvalidMesh::NoCells);
        }

        let mesh = Self { nodes, cells };
        for (cell, connectivity) in cells.iter().enumerate() {
            for (position, node) in connectivity.iter().enumerate() {
                if *node >= nodes.len() {
                    return Err(InvalidMesh::NodeIndexOutOfRange {
                        cell,
                        position,
                        node: *node,
                        nodes: nodes.len(),
                    });
                }
            }

            let coordinates = mesh.gather(connectivity);
            let element = Simplex::new(&coordinates);
            let measure = element.signed_measure();
            if !measure.is_finite() {
                return Err(InvalidMesh::NonFiniteCell { cell });
            } else if measure < <T as NumericElement>::ZERO {
                return Err(InvalidMesh::InvertedCell { cell });
            } else if measure <= <T as NumericElement>::ZERO {
                // Finite and not negative, so this is exactly zero. Written as
                // an ordering rather than an equality because comparing a
                // float to zero for equality reads as the mistake it is not.
                return Err(InvalidMesh::DegenerateCell { cell });
            }

            // Run the exact call assembly will run, rather than inferring its
            // success from the measure. The two share a determinant but not a
            // pivot sequence, so a measure test would leave the inference
            // one step short of a proof — and that step is the one the
            // `expect` in `internal_forces` discharges.
            if element.shape_gradients().is_err() {
                return Err(InvalidMesh::DegenerateCell { cell });
            }
        }
        Ok(mesh)
    }

    /// The mesh's node coordinates.
    #[must_use]
    pub const fn nodes(&self) -> &'mesh [[T; D]] {
        self.nodes
    }

    /// The mesh's cell connectivity, `D + 1` node indices per cell.
    #[must_use]
    pub const fn cells(&self) -> &'mesh [[usize; N]] {
        self.cells
    }

    /// Node count.
    #[must_use]
    pub const fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Degrees of freedom: one displacement component per node per dimension.
    #[must_use]
    pub const fn degrees_of_freedom(&self) -> usize {
        self.nodes.len() * D
    }

    /// Copy a cell's node coordinates into a contiguous element buffer.
    ///
    /// Indices are checked once in [`Self::try_new`], so every access here is
    /// in range for the life of the borrow.
    pub(crate) fn gather(&self, connectivity: &[usize; N]) -> [[T; D]; N] {
        let mut coordinates = [[<T as NumericElement>::ZERO; D]; N];
        for (slot, node) in coordinates.iter_mut().zip(connectivity.iter()) {
            *slot = self.nodes[*node];
        }
        coordinates
    }
}

/// A mesh assembly could not integrate.
///
/// Each variant names a condition that yields a plausible wrong answer rather
/// than an obvious failure, which is why construction rejects it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidMesh {
    /// The mesh has no nodes.
    NoNodes,
    /// The mesh has no cells.
    NoCells,
    /// A cell names a node the mesh does not have.
    NodeIndexOutOfRange {
        /// Offending cell.
        cell: usize,
        /// Position within that cell's connectivity.
        position: usize,
        /// The out-of-range index.
        node: usize,
        /// Nodes the mesh actually has.
        nodes: usize,
    },
    /// A cell has no measure, so its shape gradients do not exist.
    DegenerateCell {
        /// Offending cell.
        cell: usize,
    },
    /// A cell is wound the wrong way, so its stiffness enters with the wrong
    /// sign and the assembled operator is indefinite.
    InvertedCell {
        /// Offending cell.
        cell: usize,
    },
    /// A cell's coordinates are not finite, so its measure is not a number.
    NonFiniteCell {
        /// Offending cell.
        cell: usize,
    },
}

impl core::fmt::Display for InvalidMesh {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoNodes => write!(formatter, "the mesh has no nodes"),
            Self::NoCells => write!(formatter, "the mesh has no cells"),
            Self::NodeIndexOutOfRange {
                cell,
                position,
                node,
                nodes,
            } => write!(
                formatter,
                "cell {cell} position {position} names node {node}, but the mesh has {nodes}"
            ),
            Self::DegenerateCell { cell } => {
                write!(formatter, "cell {cell} is degenerate and has no measure")
            }
            Self::InvertedCell { cell } => write!(
                formatter,
                "cell {cell} is inverted; its measure is negative"
            ),
            Self::NonFiniteCell { cell } => {
                write!(formatter, "cell {cell} has a measure that is not finite")
            }
        }
    }
}

impl core::error::Error for InvalidMesh {}
