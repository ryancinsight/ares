use eunomia::{NumericElement, RealField};

/// A prescribed displacement component at one node.
///
/// The typed form of a Dirichlet condition. The alternative — callers zeroing
/// rows and columns of a matrix by index — is where sign errors and off-by-one
/// mistakes live, and it produces a system that still solves, so nothing
/// reports the error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedDisplacement<T> {
    node: usize,
    component: usize,
    value: T,
}

impl<T> PrescribedDisplacement<T> {
    /// Prescribe `value` for `component` of the displacement at `node`.
    ///
    /// Range is checked against the mesh by
    /// [`DirichletConditions::try_new`], not here: a condition is meaningful
    /// only against the mesh it constrains.
    #[must_use]
    pub const fn new(node: usize, component: usize, value: T) -> Self {
        Self {
            node,
            component,
            value,
        }
    }

    /// The constrained node.
    #[must_use]
    pub const fn node(&self) -> usize {
        self.node
    }

    /// The constrained component.
    #[must_use]
    pub const fn component(&self) -> usize {
        self.component
    }

    /// The prescribed value.
    #[must_use]
    pub const fn value(&self) -> &T {
        &self.value
    }

    /// The degree of freedom this condition constrains, in the flat field
    /// layout `node * D + component`.
    #[must_use]
    pub const fn degree_of_freedom<const D: usize>(&self) -> usize {
        self.node * D + self.component
    }
}

/// A validated set of Dirichlet conditions on one mesh.
///
/// # Ordering is a precondition, and why
///
/// The conditions must be strictly increasing in degree of freedom. That makes
/// duplicate detection a single pass rather than a quadratic scan, and it
/// fixes the order in which conditions apply — two conditions on one degree of
/// freedom would otherwise resolve by whichever the caller happened to list
/// last, which is a silent dependence on input order. Out-of-order input is
/// rejected rather than sorted, because sorting would need to allocate and the
/// caller already knows the order it built them in.
#[derive(Clone, Copy, Debug)]
pub struct DirichletConditions<'bc, T, const D: usize> {
    prescribed: &'bc [PrescribedDisplacement<T>],
}

impl<'bc, T: RealField, const D: usize> DirichletConditions<'bc, T, D> {
    /// Validate conditions against a mesh of `node_count` nodes.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidConditions`] for a node or component outside the mesh,
    /// or for conditions that are not strictly increasing in degree of
    /// freedom — which covers both duplicates and unordered input.
    pub fn try_new(
        prescribed: &'bc [PrescribedDisplacement<T>],
        node_count: usize,
    ) -> Result<Self, InvalidConditions> {
        let mut previous: Option<usize> = None;
        for (position, condition) in prescribed.iter().enumerate() {
            if condition.node >= node_count {
                return Err(InvalidConditions::NodeOutOfRange {
                    position,
                    node: condition.node,
                    nodes: node_count,
                });
            }
            if condition.component >= D {
                return Err(InvalidConditions::ComponentOutOfRange {
                    position,
                    component: condition.component,
                    dimensions: D,
                });
            }
            let dof = condition.degree_of_freedom::<D>();
            if let Some(last) = previous
                && dof <= last
            {
                return Err(InvalidConditions::NotStrictlyIncreasing {
                    position,
                    degree_of_freedom: dof,
                    previous: last,
                });
            }
            previous = Some(dof);
        }
        Ok(Self { prescribed })
    }

    /// The conditions, in degree-of-freedom order.
    #[must_use]
    pub const fn prescribed(&self) -> &'bc [PrescribedDisplacement<T>] {
        self.prescribed
    }

    /// How many degrees of freedom are constrained.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.prescribed.len()
    }

    /// Whether the mesh is unconstrained.
    ///
    /// An unconstrained elastostatic system is singular — the rigid-body modes
    /// are in its null space — so this answering `true` is the diagnosis for a
    /// solve that will not converge.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.prescribed.is_empty()
    }

    /// Zero the constrained entries of a field: the projection `P`.
    ///
    /// Costs one pass over the conditions, not over the field.
    pub fn project(&self, field: &mut [T]) {
        for condition in self.prescribed {
            if let Some(slot) = field.get_mut(condition.degree_of_freedom::<D>()) {
                *slot = <T as NumericElement>::ZERO;
            }
        }
    }

    /// Write the prescribed values into the constrained entries of a field.
    pub fn impose(&self, field: &mut [T]) {
        for condition in self.prescribed {
            if let Some(slot) = field.get_mut(condition.degree_of_freedom::<D>()) {
                *slot = condition.value;
            }
        }
    }

    /// Copy the constrained entries of `source` into `target`.
    ///
    /// The identity rows of the constrained operator: the constrained part of
    /// the output is the constrained part of the input, so those rows read
    /// `u_c = u_c` and the operator stays symmetric and positive definite on
    /// the whole space rather than only on the free subspace.
    pub fn carry(&self, source: &[T], target: &mut [T]) {
        for condition in self.prescribed {
            let dof = condition.degree_of_freedom::<D>();
            if let Some(value) = source.get(dof)
                && let Some(slot) = target.get_mut(dof)
            {
                *slot = *value;
            }
        }
    }
}

/// A Dirichlet condition set that does not describe the mesh it constrains.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidConditions {
    /// A condition names a node the mesh does not have.
    NodeOutOfRange {
        /// Position within the condition list.
        position: usize,
        /// The out-of-range node.
        node: usize,
        /// Nodes the mesh has.
        nodes: usize,
    },
    /// A condition names a component the space does not have.
    ComponentOutOfRange {
        /// Position within the condition list.
        position: usize,
        /// The out-of-range component.
        component: usize,
        /// Spatial dimensions.
        dimensions: usize,
    },
    /// The conditions are not strictly increasing in degree of freedom, so
    /// either one is duplicated or the list is unordered.
    NotStrictlyIncreasing {
        /// Position within the condition list.
        position: usize,
        /// This condition's degree of freedom.
        degree_of_freedom: usize,
        /// The previous condition's, which it does not exceed.
        previous: usize,
    },
}

impl core::fmt::Display for InvalidConditions {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NodeOutOfRange {
                position,
                node,
                nodes,
            } => write!(
                formatter,
                "condition {position} names node {node}, but the mesh has {nodes}"
            ),
            Self::ComponentOutOfRange {
                position,
                component,
                dimensions,
            } => write!(
                formatter,
                "condition {position} names component {component} of a {dimensions}-dimensional \
                 displacement"
            ),
            Self::NotStrictlyIncreasing {
                position,
                degree_of_freedom,
                previous,
            } => write!(
                formatter,
                "condition {position} constrains degree of freedom {degree_of_freedom}, which \
                 does not exceed the previous {previous}"
            ),
        }
    }
}

impl core::error::Error for InvalidConditions {}
