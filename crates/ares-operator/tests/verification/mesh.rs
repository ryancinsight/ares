//! Structured fixtures and the solve harness shared by the oracles.

use ares::{DirichletConditions, PrescribedDisplacement, SimplexMesh};
use ares_operator::ConstrainedStiffness;
use athena_core::{Cg, CgWorkspace, ConvergencePolicy, Identity};
use athena_leto::LetoBackend;
use eunomia::RealField;
use leto::Array1;
use leto_ops::RealScalar;
use proteus::IsotropicModuli;

/// Lame parameters from the engineering pair.
pub fn lame(young: f64, poisson: f64) -> (f64, f64) {
    (
        young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
        young / (2.0 * (1.0 + poisson)),
    )
}

pub fn moduli<T: RealField>(young: f64, poisson: f64) -> IsotropicModuli<T> {
    use aequitas::systems::si::quantities::{Dimensionless, Pressure};
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(T::from_f64(young)),
        Dimensionless::from_base(T::from_f64(poisson)),
    )
    .expect("inside the positive-definite domain")
}

/// A right-triangulated grid on `[0, width] x [0, height]`.
///
/// Each quad splits along one diagonal, giving a mesh that is structured but
/// not symmetric under reflection — a symmetric split would let a sign error
/// in one triangle cancel against its mirror.
pub struct Grid {
    pub nodes: Vec<[f64; 2]>,
    pub cells: Vec<[usize; 3]>,
    pub columns: usize,
    pub rows: usize,
    pub width: f64,
    pub height: f64,
}

impl Grid {
    pub fn new(width: f64, height: f64, columns: usize, rows: usize) -> Self {
        let mut nodes = Vec::with_capacity((columns + 1) * (rows + 1));
        for row in 0..=rows {
            for column in 0..=columns {
                nodes.push([
                    width * column as f64 / columns as f64,
                    height * row as f64 / rows as f64,
                ]);
            }
        }
        let mut cells = Vec::with_capacity(columns * rows * 2);
        for row in 0..rows {
            for column in 0..columns {
                let index = |c: usize, r: usize| r * (columns + 1) + c;
                let (a, b) = (index(column, row), index(column + 1, row));
                let (c, d) = (index(column + 1, row + 1), index(column, row + 1));
                cells.push([a, b, c]);
                cells.push([a, c, d]);
            }
        }
        Self {
            nodes,
            cells,
            columns,
            rows,
            width,
            height,
        }
    }

    pub fn node_index(&self, column: usize, row: usize) -> usize {
        row * (self.columns + 1) + column
    }

    pub fn degrees_of_freedom(&self) -> usize {
        self.nodes.len() * 2
    }

    /// The largest cell diameter, which is the `h` of a convergence study.
    pub fn cell_size(&self) -> f64 {
        let dx = self.width / self.columns as f64;
        let dy = self.height / self.rows as f64;
        dx.hypot(dy)
    }

    /// Nodes on the boundary of the rectangle.
    pub fn boundary_nodes(&self) -> Vec<usize> {
        let mut nodes = Vec::new();
        for row in 0..=self.rows {
            for column in 0..=self.columns {
                if row == 0 || row == self.rows || column == 0 || column == self.columns {
                    nodes.push(self.node_index(column, row));
                }
            }
        }
        nodes
    }
}

/// Solve `A u = b` through Athena's conjugate gradients.
///
/// Panics rather than returning on non-convergence: an oracle that silently
/// compared an unconverged field against a closed form would be measuring the
/// solver's failure as if it were the discretisation's error.
pub fn solve<T: RealScalar + RealField + core::fmt::LowerExp>(
    mesh: &SimplexMesh<'_, T, 2, 3>,
    material: IsotropicModuli<T>,
    prescribed: &[PrescribedDisplacement<T>],
    external: &[T],
) -> Vec<T> {
    let conditions =
        DirichletConditions::try_new(prescribed, mesh.node_count()).expect("valid conditions");
    let operator = ConstrainedStiffness::new(*mesh, material, conditions);
    let dofs = mesh.degrees_of_freedom();

    let mut load = vec![<T as eunomia::NumericElement>::ZERO; dofs];
    operator.load(external, &mut load).expect("well-shaped");

    let backend = LetoBackend::<T>::default();
    let right_hand_side = Array1::from_shape_vec([dofs], load).expect("valid vector");
    let mut solution = Array1::zeros([dofs]);
    let mut workspace = CgWorkspace::new(&backend, dofs).expect("workspace");
    let zero = <T as eunomia::NumericElement>::ZERO;

    // The estimate only sizes `|u|` for the threshold below, so its tolerance
    // carries no accuracy claim; `sqrt(EPSILON)` keeps it attainable far
    // beyond the conditioning of these meshes.
    let estimate =
        ConvergencePolicy::<T>::new(zero, T::EPSILON.sqrt(), SOLVE_BUDGET).expect("valid policy");
    converge(
        &operator,
        &right_hand_side,
        &mut solution,
        &mut workspace,
        estimate,
    );

    let threshold = attainable_threshold(
        row_terms(mesh),
        infinity_norm(&operator),
        euclidean_norm(solution.as_slice().expect("contiguous")),
    );
    let policy = ConvergencePolicy::<T>::new(threshold, zero, SOLVE_BUDGET).expect("valid policy");
    converge(
        &operator,
        &right_hand_side,
        &mut solution,
        &mut workspace,
        policy,
    );
    solution.as_slice().expect("contiguous").to_vec()
}

/// Iteration budget of one verification solve.
///
/// It bounds a stalled solve; a converging one stops on its threshold first.
const SOLVE_BUDGET: usize = 20_000;

/// Run CG from the current `solution` and panic unless it converged.
fn converge<T: RealScalar + RealField + core::fmt::LowerExp>(
    operator: &ConstrainedStiffness<'_, T, 2, 3>,
    right_hand_side: &Array1<T>,
    solution: &mut Array1<T>,
    workspace: &mut CgWorkspace<LetoBackend<T>>,
    policy: ConvergencePolicy<T>,
) {
    let report = Cg::<LetoBackend<T>>::solve_into(
        &LetoBackend::<T>::default(),
        operator,
        &Identity,
        right_hand_side,
        solution,
        workspace,
        policy,
    )
    .expect("well-shaped system");
    assert!(
        report.converged(),
        "conjugate gradients did not converge: {:?} after {} iterations, residual {:.3e} against \
         threshold {:.3e}",
        report.termination,
        report.iterations,
        report.final_residual_norm,
        report.threshold
    );
}

/// Smallest residual threshold that one evaluation of `b - A u` can certify.
///
/// Athena accepts convergence only when the recomputed `|b - A u|_2` meets the
/// threshold, so a threshold below the rounding error of that evaluation is
/// unattainable whatever the iteration does. Forming `A u` accumulates
/// `row_terms - 1` products per entry and the subtraction from `b` adds one
/// more operation, so the computed residual differs from the exact one by
/// `|dr| <= gamma_k (|b| + |A| |u|)` componentwise, with
/// `gamma_k = k u / (1 - k u)`, `k = row_terms` and unit roundoff
/// `u = EPSILON / 2` (Higham, *Accuracy and Stability of Numerical
/// Algorithms*, 2nd ed., Lemma 3.1 and §3.5).
///
/// Taking 2-norms, `|A|` is symmetric with `||A|||_1 = |A|_inf`, so
/// `||A| |u||_2 <= |A|_inf |u|_2`, and `|b|_2 = |A u|_2 <= |A|_2 |u|_2 <=
/// |A|_inf |u|_2`. The floor is therefore at most `2 gamma_k |A|_inf |u|_2`,
/// which is what this returns. Relative to `|b|` it is `2 gamma_k |A|_inf
/// |u| / |b|`, a factor bounded by `2 gamma_k kappa(A)`-scale terms: a
/// relative tolerance fixed at a multiple of `EPSILON` is below it whenever
/// `|A| |u| / |b|` exceeds that multiple, which the cantilever (about `1e5`)
/// and the cylinders (about `1e2` to `1e3`) do.
///
/// Meeting this threshold makes the normwise backward error
/// `|b - A u| / (|A| |u|)` at most `2 gamma_k` (Rigal and Gaches; Higham,
/// Theorem 7.1), so the solution error is bounded by
/// `2 gamma_k kappa(A)` relative, far below every discretisation error these
/// oracles measure.
///
/// # Limits
///
/// `u` is the estimate solve's iterate, whose norm is below the exact
/// solution's because CG norms grow monotonically from a zero start
/// (Hestenes and Stiefel 1952, Theorem 6:3). An underestimate only tightens
/// the threshold, so a poor estimate fails the solve loudly and never passes
/// it falsely. The bound treats the matrix-free product as an assembled one;
/// element-matrix formation rounding is outside `row_terms`.
fn attainable_threshold<T: RealField>(row_terms: usize, operator_norm: T, solution_norm: T) -> T {
    let terms = T::from_f64(f64::from(
        u32::try_from(row_terms)
            .expect("invariant: a 2D mesh row couples far fewer than 2^32 terms"),
    ));
    let one = <T as eunomia::NumericElement>::ONE;
    let unit_roundoff = T::EPSILON / (one + one);
    let gamma = terms * unit_roundoff / (one - terms * unit_roundoff);
    (one + one) * gamma * operator_norm * solution_norm
}

/// Terms accumulated into one entry of `b - A u`.
///
/// Each cell sharing a node contributes one row of its `D N` element matrix,
/// and the subtraction from `b` is one more operation.
fn row_terms<T: RealField>(mesh: &SimplexMesh<'_, T, 2, 3>) -> usize {
    let mut cells_at_node = vec![0_usize; mesh.node_count()];
    for connectivity in mesh.cells() {
        for &node in connectivity {
            cells_at_node[node] += 1;
        }
    }
    let busiest = cells_at_node.into_iter().max().unwrap_or(0);
    busiest * 2 * 3 + 1
}

/// `|A|_inf`, the largest absolute row sum, from `A`'s columns.
///
/// The constrained stiffness is symmetric, so its largest absolute column sum
/// equals its largest absolute row sum. One application per column is exact
/// in the entries the solve itself uses, where a power iteration would
/// approach `|A|_2` from below.
fn infinity_norm<T: RealScalar + RealField>(operator: &ConstrainedStiffness<'_, T, 2, 3>) -> T {
    use athena_core::{KrylovBackend, LinearOperator};
    let backend = LetoBackend::<T>::default();
    let dofs = operator.dimension();
    let zero = <T as eunomia::NumericElement>::ZERO;
    let mut unit = Array1::zeros([dofs]);
    let mut column = Array1::zeros([dofs]);
    let mut largest = zero;
    for index in 0..dofs {
        {
            let mut view = backend.view_mut(&mut unit);
            let slots = view.as_mut_slice().expect("contiguous");
            slots.fill(zero);
            slots[index] = <T as eunomia::NumericElement>::ONE;
        }
        operator
            .apply(&backend, backend.view(&unit), backend.view_mut(&mut column))
            .expect("well-shaped");
        let sum = column
            .as_slice()
            .expect("contiguous")
            .iter()
            .fold(zero, |total, entry| total + entry.abs());
        if sum > largest {
            largest = sum;
        }
    }
    largest
}

fn euclidean_norm<T: RealField>(values: &[T]) -> T {
    values
        .iter()
        .fold(<T as eunomia::NumericElement>::ZERO, |total, value| {
            total + *value * *value
        })
        .sqrt()
}

/// The mass-weighted `L2` norm of a nodal field over a mesh.
///
/// Uses the consistent mass matrix `integral(N_a N_b) dV`, so this is the
/// exact `L2` norm of the field's linear interpolant rather than a discrete
/// sum of nodal values — a plain root-mean-square would weight a coarse mesh's
/// nodes the same as a fine one's and report a rate that depends on the node
/// count instead of the element.
pub fn l2_norm<T: RealField>(mesh: &SimplexMesh<'_, T, 2, 3>, field: &[T]) -> f64
where
    f64: From<T>,
{
    let mut total = 0.0_f64;
    for connectivity in mesh.cells() {
        let mut coordinates = [[<T as eunomia::NumericElement>::ZERO; 2]; 3];
        for (slot, node) in coordinates.iter_mut().zip(connectivity.iter()) {
            *slot = mesh.nodes()[*node];
        }
        let measure = f64::from(ares::Simplex::new(&coordinates).signed_measure());
        for component in 0..2 {
            let values: Vec<f64> = connectivity
                .iter()
                .map(|node| f64::from(field[node * 2 + component]))
                .collect();
            let sum: f64 = values.iter().sum();
            let squares: f64 = values.iter().map(|v| v * v).sum();
            // sum_ab (1 + delta_ab) v_a v_b = (sum v)^2 + sum v^2, over 12.
            total += measure * (sum * sum + squares) / 12.0;
        }
    }
    total.sqrt()
}
