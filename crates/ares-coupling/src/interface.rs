use ares::{SimplexMesh, TractionFacet};
use eunomia::{NumericElement, RealField};

/// The structural side of a coupling interface: which facets receive traction
/// and which nodes report displacement.
///
/// # The exchange ordering contract
///
/// Harmonia's `Partition` exchanges flat `&[T]`, so the two sides must agree
/// on what position `k` means. This type fixes it:
///
/// - **input** — traction, **facet index major, component minor**:
///   `input[f * D + c]`.
/// - **output** — displacement, **node index major, component minor**:
///   `output[k * D + c]`, where `k` indexes [`Self::nodes`], not the mesh.
///
/// The two orderings differ because the two quantities live on different
/// entities, and that asymmetry is deliberate rather than an oversight. A
/// traction is a stress resolved on a *surface*, and the fluid side computes
/// one per face from the flow state either side of it. A displacement is a
/// property of a *point*, and the structural solve carries one per node.
/// Forcing either onto the other's entity would mean interpolating, which
/// atlas ADR 0050 explicitly places outside Harmonia.
///
/// Atlas ADR 0059 states the contract as "interface node index major,
/// component minor". That describes the displacement side exactly and the
/// traction side not at all; the record is worth correcting rather than the
/// code bent to match it.
///
/// # Why per-facet traction makes the work balance exact
///
/// With traction piecewise constant per facet and displacement linear over it,
///
/// ```text
/// integral(t . u) dS = sum_f t_f |A_f| (mean of u over f's nodes)
///                    = sum_f sum_a t_f |A_f| u_a / D
///                    = sum_a u_a . f_a
/// ```
///
/// where `f_a = t_f |A_f| / D` is exactly the consistent nodal load Ares
/// already computes. The interface work identity is therefore an equality
/// rather than an approximation — and it is an equality *because* the load is
/// the consistent one. A lumped load would carry the same total force and
/// break this.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInterface<'mesh, const D: usize> {
    nodes: &'mesh [usize],
    facets: &'mesh [[usize; D]],
}

impl<'mesh, const D: usize> StructuralInterface<'mesh, D> {
    /// Validate an interface against the mesh it couples.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidInterface`] when the interface is not conforming: a
    /// node outside the mesh, a duplicate in the exchange ordering, or a facet
    /// naming a node the interface does not carry. Phase 0 requires a
    /// conforming interface (atlas ADR 0059), so a non-conforming one is
    /// rejected here rather than silently transferred — the failure the
    /// conformity oracle exists to force.
    pub fn try_new<T: RealField, const N: usize>(
        nodes: &'mesh [usize],
        facets: &'mesh [[usize; D]],
        mesh: &SimplexMesh<'_, T, D, N>,
    ) -> Result<Self, InvalidInterface> {
        if nodes.is_empty() {
            return Err(InvalidInterface::NoNodes);
        }
        if facets.is_empty() {
            return Err(InvalidInterface::NoFacets);
        }

        for (position, node) in nodes.iter().enumerate() {
            if *node >= mesh.node_count() {
                return Err(InvalidInterface::NodeOutOfRange {
                    position,
                    node: *node,
                    nodes: mesh.node_count(),
                });
            }
            if nodes[..position].contains(node) {
                return Err(InvalidInterface::DuplicateNode {
                    position,
                    node: *node,
                });
            }
        }

        for (facet, connectivity) in facets.iter().enumerate() {
            for (position, node) in connectivity.iter().enumerate() {
                if !nodes.contains(node) {
                    return Err(InvalidInterface::FacetNodeNotOnInterface {
                        facet,
                        position,
                        node: *node,
                    });
                }
            }
        }

        Ok(Self { nodes, facets })
    }

    /// Interface nodes, in exchange order.
    #[must_use]
    pub const fn nodes(&self) -> &'mesh [usize] {
        self.nodes
    }

    /// Interface facets, in exchange order.
    #[must_use]
    pub const fn facets(&self) -> &'mesh [[usize; D]] {
        self.facets
    }

    /// Length of an incoming traction exchange: `facets * D`.
    #[must_use]
    pub const fn input_dimension(&self) -> usize {
        self.facets.len() * D
    }

    /// Length of an outgoing displacement exchange: `nodes * D`.
    #[must_use]
    pub const fn output_dimension(&self) -> usize {
        self.nodes.len() * D
    }

    /// Read a flat traction exchange into Ares traction facets.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidInterface::MisshapedExchange`] when `traction` is not
    /// [`input_dimension`](Self::input_dimension) long.
    pub fn read_traction<T: RealField>(
        &self,
        traction: &[T],
        out: &mut [TractionFacet<T, D>],
    ) -> Result<(), InvalidInterface> {
        if traction.len() != self.input_dimension() {
            return Err(InvalidInterface::MisshapedExchange {
                expected: self.input_dimension(),
                found: traction.len(),
            });
        }
        if out.len() != self.facets.len() {
            return Err(InvalidInterface::MisshapedExchange {
                expected: self.facets.len(),
                found: out.len(),
            });
        }
        let (per_facet, _) = traction.as_chunks::<D>();
        for ((slot, connectivity), value) in
            out.iter_mut().zip(self.facets.iter()).zip(per_facet.iter())
        {
            *slot = TractionFacet::new(*connectivity, *value);
        }
        Ok(())
    }

    /// Write interface displacement out of a full nodal state.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidInterface::MisshapedExchange`] when either slice is
    /// the wrong length.
    pub fn write_displacement<T: RealField>(
        &self,
        state: &[T],
        out: &mut [T],
    ) -> Result<(), InvalidInterface> {
        if out.len() != self.output_dimension() {
            return Err(InvalidInterface::MisshapedExchange {
                expected: self.output_dimension(),
                found: out.len(),
            });
        }
        let (nodal_state, remainder) = state.as_chunks::<D>();
        if !remainder.is_empty() {
            return Err(InvalidInterface::MisshapedExchange {
                expected: nodal_state.len() * D,
                found: state.len(),
            });
        }
        let (nodal_out, _) = out.as_chunks_mut::<D>();
        for (slot, node) in nodal_out.iter_mut().zip(self.nodes.iter()) {
            let Some(value) = nodal_state.get(*node) else {
                return Err(InvalidInterface::NodeOutOfRange {
                    position: *node,
                    node: *node,
                    nodes: nodal_state.len(),
                });
            };
            *slot = *value;
        }
        Ok(())
    }

    /// The work `integral(t . u) dS` the interface traction does on a nodal
    /// displacement state.
    ///
    /// Computed from the facet integral rather than from the assembled load
    /// vector, so comparing it against the structural strain energy tests the
    /// coupling instead of restating it. Reading it off the load vector would
    /// make the interface-work oracle an identity about one array.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidInterface::MisshapedExchange`] on a length mismatch.
    pub fn interface_work<T: RealField, const N: usize>(
        &self,
        mesh: &SimplexMesh<'_, T, D, N>,
        traction: &[T],
        state: &[T],
    ) -> Result<T, InvalidInterface> {
        if traction.len() != self.input_dimension() {
            return Err(InvalidInterface::MisshapedExchange {
                expected: self.input_dimension(),
                found: traction.len(),
            });
        }
        let (nodal_state, _) = state.as_chunks::<D>();
        let (per_facet, _) = traction.as_chunks::<D>();

        let mut divisor = <T as NumericElement>::ZERO;
        for _ in 0..D {
            divisor += <T as NumericElement>::ONE;
        }

        let mut work = <T as NumericElement>::ZERO;
        for (connectivity, value) in self.facets.iter().zip(per_facet.iter()) {
            let mut coordinates = [[<T as NumericElement>::ZERO; D]; D];
            for (slot, node) in coordinates.iter_mut().zip(connectivity.iter()) {
                *slot = mesh.nodes()[*node];
            }
            let measure = ares::boundary::facet_measure(&coordinates);

            // The facet's mean displacement, which is what a constant traction
            // integrates against over a linear element.
            let mut mean = [<T as NumericElement>::ZERO; D];
            for node in connectivity {
                for (slot, component) in mean.iter_mut().zip(nodal_state[*node].iter()) {
                    *slot += *component;
                }
            }
            for (component, average) in value.iter().zip(mean.iter()) {
                work += *component * *average * measure / divisor;
            }
        }
        Ok(work)
    }
}

/// An interface that does not conform to the mesh it couples.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidInterface {
    /// The interface carries no nodes.
    NoNodes,
    /// The interface carries no facets.
    NoFacets,
    /// An interface node is not a node of the mesh.
    NodeOutOfRange {
        /// Position in the exchange ordering.
        position: usize,
        /// The out-of-range node.
        node: usize,
        /// Nodes the mesh has.
        nodes: usize,
    },
    /// A node appears twice in the exchange ordering, so an exchange position
    /// would be ambiguous.
    DuplicateNode {
        /// The later position.
        position: usize,
        /// The repeated node.
        node: usize,
    },
    /// A facet names a node the interface does not carry, so the two sides do
    /// not describe the same surface — the non-conforming case ADR 0059
    /// excludes from Phase 0.
    FacetNodeNotOnInterface {
        /// Offending facet.
        facet: usize,
        /// Position within it.
        position: usize,
        /// The node that is not on the interface.
        node: usize,
    },
    /// An exchange buffer is the wrong length.
    MisshapedExchange {
        /// Entries the interface expects.
        expected: usize,
        /// Entries supplied.
        found: usize,
    },
}

impl core::fmt::Display for InvalidInterface {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoNodes => write!(formatter, "the interface carries no nodes"),
            Self::NoFacets => write!(formatter, "the interface carries no facets"),
            Self::NodeOutOfRange {
                position,
                node,
                nodes,
            } => write!(
                formatter,
                "interface position {position} names node {node}, but the mesh has {nodes}"
            ),
            Self::DuplicateNode { position, node } => write!(
                formatter,
                "node {node} appears again at interface position {position}"
            ),
            Self::FacetNodeNotOnInterface {
                facet,
                position,
                node,
            } => write!(
                formatter,
                "facet {facet} position {position} names node {node}, which is not on the interface"
            ),
            Self::MisshapedExchange { expected, found } => write!(
                formatter,
                "the exchange buffer has {found} entries, but the interface expects {expected}"
            ),
        }
    }
}

impl core::error::Error for InvalidInterface {}
