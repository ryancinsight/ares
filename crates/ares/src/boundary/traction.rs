use eunomia::{NumericElement, RealField};

use crate::element::leading_determinant;

/// A uniform traction applied over one boundary facet.
///
/// A facet of a `D`-simplex is a `(D - 1)`-simplex with exactly `D` nodes — an
/// edge in 2-D, a triangle in 3-D — so both arrays are `D` long and the shape
/// needs no further parameter.
///
/// `traction` is force per unit area (per unit length in 2-D): a stress, not a
/// force. The distinction is the one Neumann conditions are usually got wrong
/// on, because the two differ by the facet measure and a small mesh makes the
/// error look like a modelling choice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TractionFacet<T, const D: usize> {
    nodes: [usize; D],
    traction: [T; D],
}

impl<T, const D: usize> TractionFacet<T, D> {
    /// Apply `traction` uniformly over the facet spanned by `nodes`.
    #[must_use]
    pub const fn new(nodes: [usize; D], traction: [T; D]) -> Self {
        Self { nodes, traction }
    }

    /// The facet's nodes.
    #[must_use]
    pub const fn nodes(&self) -> &[usize; D] {
        &self.nodes
    }

    /// The applied traction.
    #[must_use]
    pub const fn traction(&self) -> &[T; D] {
        &self.traction
    }
}

/// The `(D - 1)`-measure of a facet: length in 2-D, area in 3-D.
///
/// Computed as `sqrt(det(E^T E)) / (D - 1)!` from the facet's edge matrix `E`,
/// whose `D - 1` columns are `x_i - x_0`. The Gram determinant is the general
/// form of the cross-product magnitude that gives a triangle's area and the
/// edge-length that gives a segment's, so one expression covers both
/// dimensions instead of a match on `D`.
#[must_use]
pub fn facet_measure<T: RealField, const D: usize>(nodes: &[[T; D]; D]) -> T {
    let edges = D - 1;

    // Column `i` is `x_{i+1} - x_0`; only the leading `edges` columns are used.
    let mut edge = [[<T as NumericElement>::ZERO; D]; D];
    for (column, node) in nodes.iter().skip(1).enumerate() {
        for (row, target) in edge.iter_mut().enumerate() {
            target[column] = node[row] - nodes[0][row];
        }
    }

    // Gram matrix G[i][j] = e_i . e_j, again in the leading block.
    let mut gram = [[<T as NumericElement>::ZERO; D]; D];
    for i in 0..edges {
        for j in 0..edges {
            let mut sum = <T as NumericElement>::ZERO;
            for row in &edge {
                sum += row[i] * row[j];
            }
            gram[i][j] = sum;
        }
    }

    let mut factorial = <T as NumericElement>::ONE;
    let mut term = <T as NumericElement>::ONE;
    for _ in 1..=edges {
        factorial *= term;
        term += <T as NumericElement>::ONE;
    }
    let squared = leading_determinant(&gram, edges);
    if squared > <T as NumericElement>::ZERO {
        squared.sqrt() / factorial
    } else {
        // A Gram determinant is non-negative in exact arithmetic, so a
        // non-positive value is a collapsed facet rather than a signed one.
        // Returning zero keeps the caller's rejection test a comparison
        // instead of a `NaN` check.
        <T as NumericElement>::ZERO
    }
}

/// A validated set of traction facets on one mesh.
#[derive(Clone, Copy, Debug)]
pub struct TractionBoundary<'bc, T, const D: usize> {
    facets: &'bc [TractionFacet<T, D>],
}

impl<'bc, T: RealField, const D: usize> TractionBoundary<'bc, T, D> {
    /// Validate facets against a mesh's nodes.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidBoundary`] for a node outside the mesh, a repeated
    /// node within one facet, or a facet with no measure — each of which would
    /// otherwise contribute a load that is silently zero or infinite.
    pub fn try_new(
        facets: &'bc [TractionFacet<T, D>],
        nodes: &[[T; D]],
    ) -> Result<Self, InvalidBoundary> {
        for (facet, description) in facets.iter().enumerate() {
            for (position, node) in description.nodes.iter().enumerate() {
                if *node >= nodes.len() {
                    return Err(InvalidBoundary::NodeOutOfRange {
                        facet,
                        position,
                        node: *node,
                        nodes: nodes.len(),
                    });
                }
            }

            let mut coordinates = [[<T as NumericElement>::ZERO; D]; D];
            for (slot, node) in coordinates.iter_mut().zip(description.nodes.iter()) {
                *slot = nodes[*node];
            }
            let measure = facet_measure(&coordinates);
            if !measure.is_finite() {
                return Err(InvalidBoundary::NonFiniteFacet { facet });
            }
            if measure <= <T as NumericElement>::ZERO {
                return Err(InvalidBoundary::DegenerateFacet { facet });
            }
        }
        Ok(Self { facets })
    }

    /// The facets, in the order given.
    #[must_use]
    pub const fn facets(&self) -> &'bc [TractionFacet<T, D>] {
        self.facets
    }

    /// Add the consistent nodal loads for every facet into `loads`.
    ///
    /// Adds rather than assigns, so tractions compose with a body force and
    /// with each other.
    ///
    /// # The load a facet contributes
    ///
    /// For linear shape functions on a `(D - 1)`-simplex,
    /// `integral(N_a) dS = measure / D` for each of the `D` facet nodes, so a
    /// uniform traction distributes equally:
    ///
    /// ```text
    /// f_a = t * measure / D
    /// ```
    ///
    /// Equal distribution is a property of *linear* elements and a uniform
    /// traction, not a simplification. A quadratic element's consistent load
    /// is famously unequal — and negative at the corners — which is why this
    /// is derived from the shape-function integral rather than assumed to be
    /// "the traction split between the nodes".
    ///
    /// # Errors
    ///
    /// Returns [`MisshapedLoad`] when `loads` is not `nodes.len() * D` long.
    pub fn add_consistent_loads(
        &self,
        nodes: &[[T; D]],
        loads: &mut [T],
    ) -> Result<(), MisshapedLoad> {
        let expected = nodes.len() * D;
        if loads.len() != expected {
            return Err(MisshapedLoad {
                expected,
                found: loads.len(),
            });
        }

        let mut divisor = <T as NumericElement>::ZERO;
        for _ in 0..D {
            divisor += <T as NumericElement>::ONE;
        }

        let (nodal, _) = loads.as_chunks_mut::<D>();
        for facet in self.facets {
            let mut coordinates = [[<T as NumericElement>::ZERO; D]; D];
            for (slot, node) in coordinates.iter_mut().zip(facet.nodes.iter()) {
                *slot = nodes[*node];
            }
            let share = facet_measure(&coordinates) / divisor;
            for node in &facet.nodes {
                for (slot, component) in nodal[*node].iter_mut().zip(facet.traction.iter()) {
                    *slot += *component * share;
                }
            }
        }
        Ok(())
    }
}

/// A traction boundary that does not describe the mesh it loads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidBoundary {
    /// A facet names a node the mesh does not have.
    NodeOutOfRange {
        /// Offending facet.
        facet: usize,
        /// Position within that facet.
        position: usize,
        /// The out-of-range node.
        node: usize,
        /// Nodes the mesh has.
        nodes: usize,
    },
    /// A facet has no measure, so it would carry no load however large its
    /// traction.
    DegenerateFacet {
        /// Offending facet.
        facet: usize,
    },
    /// A facet's coordinates are not finite.
    NonFiniteFacet {
        /// Offending facet.
        facet: usize,
    },
}

impl core::fmt::Display for InvalidBoundary {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NodeOutOfRange {
                facet,
                position,
                node,
                nodes,
            } => write!(
                formatter,
                "facet {facet} position {position} names node {node}, but the mesh has {nodes}"
            ),
            Self::DegenerateFacet { facet } => {
                write!(formatter, "facet {facet} has no measure")
            }
            Self::NonFiniteFacet { facet } => {
                write!(formatter, "facet {facet} has a measure that is not finite")
            }
        }
    }
}

impl core::error::Error for InvalidBoundary {}

/// A load vector whose length does not match the mesh.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MisshapedLoad {
    /// Degrees of freedom the mesh has.
    pub expected: usize,
    /// Length supplied.
    pub found: usize,
}

impl core::fmt::Display for MisshapedLoad {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "the load vector has {} entries, but the mesh has {} degrees of freedom",
            self.found, self.expected
        )
    }
}

impl core::error::Error for MisshapedLoad {}
