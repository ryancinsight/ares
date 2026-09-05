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
    // Tight enough that the solver's residual is far below the discretisation
    // error every oracle here measures, so a convergence rate reflects the
    // element rather than the stopping rule.
    // Scaled to the scalar rather than fixed: an f64 tolerance is unreachable
    // in f32 and the solve would exhaust its budget instead of converging.
    let relative = T::EPSILON * T::from_f64(64.0);
    let policy =
        ConvergencePolicy::<T>::new(T::from_f64(0.0), relative, 20_000).expect("valid policy");

    let report = Cg::<LetoBackend<T>>::solve_into(
        &backend,
        &operator,
        &Identity,
        &right_hand_side,
        &mut solution,
        &mut workspace,
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
    solution.as_slice().expect("contiguous").to_vec()
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
