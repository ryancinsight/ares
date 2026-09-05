//! Oracles for the structural coupling partition (atlas ADR 0059, phase A8).
//!
//! # The headline: interface work is conserved
//!
//! The work the fluid does on the interface must equal the strain energy the
//! structure stores. It is a conservation identity, so it holds on any mesh at
//! any resolution — unlike the closed-form comparisons, which hold only in the
//! limit — and it is the one oracle that tests the *coupling* rather than
//! either side of it.
//!
//! The two quantities are computed by routes that share no arithmetic. The
//! work comes from the facet integral `sum_f t_f |A_f| (mean u over f)`,
//! evaluated from the traction exchange and the node coordinates. The strain
//! energy comes from `(1/2) u . K u`, evaluated by assembling the stiffness
//! over every cell in the mesh. Nothing links them but the physics.
//!
//! # Why it is an equality rather than a bound
//!
//! With traction constant per facet and displacement linear over it,
//! `integral(t . u) dS` collapses to `sum_a u_a . f_a` where `f_a` is exactly
//! the consistent nodal load. So the identity is exact — *because* the load is
//! the consistent one. A lumped load carries the same resultant force and
//! breaks it, which is the mutation this oracle is measured against.

#![expect(
    clippy::cast_precision_loss,
    reason = "grid indices and division counts become the reals of a mesh coordinate, far inside f64's exact-integer range."
)]
#![expect(
    clippy::float_cmp,
    reason = "the exact comparisons are exact by derivation. A zero traction assembles an identically zero load, so every Krylov iterate stays exactly zero and the displacement is bit-zero rather than small. The exported values are copied out of the state rather than computed, so they are bit-identical to it."
)]

use ares::{DirichletConditions, PrescribedDisplacement, SimplexMesh, TractionFacet};
use ares_coupling::{InvalidInterface, StructuralInterface, StructuralPartition};
use ares_operator::ConstrainedStiffness;
use athena_core::ConvergencePolicy;
use harmonia::Partition as _;
use proteus::IsotropicModuli;

fn moduli(young: f64, poisson: f64) -> IsotropicModuli<f64> {
    use aequitas::systems::si::quantities::{Dimensionless, Pressure};
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(young),
        Dimensionless::from_base(poisson),
    )
    .expect("inside the positive-definite domain")
}

/// A 2-by-1 strip of quads split into triangles, clamped at `x = 0`, with the
/// `x = 2` face as the coupling interface.
struct Fixture {
    nodes: Vec<[f64; 2]>,
    cells: Vec<[usize; 3]>,
    interface_nodes: Vec<usize>,
    interface_facets: Vec<[usize; 2]>,
    clamped: Vec<PrescribedDisplacement<f64>>,
}

impl Fixture {
    fn new(columns: usize, rows: usize) -> Self {
        let (width, height) = (2.0_f64, 1.0_f64);
        let index = |c: usize, r: usize| r * (columns + 1) + c;
        let mut nodes = Vec::new();
        for row in 0..=rows {
            for column in 0..=columns {
                nodes.push([
                    width * column as f64 / columns as f64,
                    height * row as f64 / rows as f64,
                ]);
            }
        }
        let mut cells = Vec::new();
        for row in 0..rows {
            for column in 0..columns {
                let (a, b) = (index(column, row), index(column + 1, row));
                let (c, d) = (index(column + 1, row + 1), index(column, row + 1));
                cells.push([a, b, c]);
                cells.push([a, c, d]);
            }
        }
        let interface_nodes: Vec<usize> = (0..=rows).map(|r| index(columns, r)).collect();
        let interface_facets: Vec<[usize; 2]> = (0..rows)
            .map(|r| [index(columns, r), index(columns, r + 1)])
            .collect();
        let mut clamped = Vec::new();
        for row in 0..=rows {
            let node = index(0, row);
            clamped.push(PrescribedDisplacement::new(node, 0, 0.0));
            clamped.push(PrescribedDisplacement::new(node, 1, 0.0));
        }
        clamped.sort_by_key(ares::PrescribedDisplacement::degree_of_freedom::<2>);
        Self {
            nodes,
            cells,
            interface_nodes,
            interface_facets,
            clamped,
        }
    }
}

fn policy() -> ConvergencePolicy<f64> {
    ConvergencePolicy::new(0.0, 1e-14, 20_000).expect("valid policy")
}

/// Run one coupling step and return `(partition, state, traction)`.
fn couple(
    fixture: &Fixture,
    traction_per_facet: [f64; 2],
) -> (StructuralPartition<'_, f64, 2, 3>, Vec<f64>, Vec<f64>) {
    let mesh = SimplexMesh::try_new(&fixture.nodes, &fixture.cells).expect("valid mesh");
    let conditions =
        DirichletConditions::try_new(&fixture.clamped, mesh.node_count()).expect("valid");
    let operator = ConstrainedStiffness::new(mesh, moduli(200e9, 0.3), conditions);
    let interface =
        StructuralInterface::try_new(&fixture.interface_nodes, &fixture.interface_facets, &mesh)
            .expect("conforming interface");
    let mut partition =
        StructuralPartition::try_new(mesh, operator, interface, policy()).expect("workspace");

    let mut traction = vec![0.0_f64; partition.input_dimension()];
    for (facet, slot) in traction.as_chunks_mut::<2>().0.iter_mut().enumerate() {
        let _ = facet;
        *slot = traction_per_facet;
    }
    let mut state = vec![0.0_f64; partition.state_dimension()];
    partition
        .solve_for_traction(&mut state, &traction)
        .expect("the coupling step solves");
    (partition, state, traction)
}

#[test]
fn interface_work_equals_the_strain_energy() {
    let fixture = Fixture::new(6, 3);
    let (partition, state, traction) = couple(&fixture, [3.0e6, 1.2e6]);

    let work = partition
        .interface_work(&traction, &state)
        .expect("well-shaped");
    let energy = partition.strain_energy(&state).expect("well-shaped");

    assert!(
        energy > 0.0,
        "the structure stored no energy, so the identity is satisfied trivially"
    );
    // The traction does work `integral(t . u)`; the structure stores half of
    // it as strain energy, the other half being the work of the linear
    // response. `W = 2 U` for a linear elastic static solve.
    let relative = (work - 2.0 * energy).abs() / work.abs();
    assert!(
        relative < 1e-9,
        "interface work {work:.6e} and strain energy {energy:.6e} violate W = 2U by \
         {relative:.3e} relative, far above the solver tolerance"
    );
}

#[test]
fn a_zero_traction_produces_exactly_zero_displacement() {
    // The coupling must add no loading of its own. A spurious constant would
    // be invisible against a real traction and dominate a small one.
    let fixture = Fixture::new(4, 2);
    let (_, state, _) = couple(&fixture, [0.0, 0.0]);
    for (dof, value) in state.iter().enumerate() {
        assert_eq!(*value, 0.0, "dof {dof} moved under zero traction");
    }
}

#[test]
fn the_exported_displacement_matches_the_interface_nodes() {
    let fixture = Fixture::new(4, 2);
    let (partition, state, _) = couple(&fixture, [2.0e6, 0.0]);

    let mut exported = vec![0.0_f64; partition.output_dimension()];
    partition
        .export(&state, &mut exported)
        .expect("well-shaped");

    for (position, node) in fixture.interface_nodes.iter().enumerate() {
        for component in 0..2 {
            assert_eq!(
                exported[position * 2 + component],
                state[node * 2 + component],
                "exchange position {position} does not carry node {node}"
            );
        }
    }
    assert!(
        exported.iter().any(|v| v.abs() > 0.0),
        "the interface did not move, so the export check is vacuous"
    );
}

#[test]
fn a_non_conforming_interface_is_rejected() {
    // Phase 0 requires a conforming interface (ADR 0059). A facet naming a
    // node the interface does not carry is the non-conforming case, and it
    // must fail with a typed error rather than transfer silently.
    let fixture = Fixture::new(4, 2);
    let mesh = SimplexMesh::try_new(&fixture.nodes, &fixture.cells).expect("valid mesh");

    let stray = vec![[fixture.interface_nodes[0], 0]];
    assert_eq!(
        StructuralInterface::try_new(&fixture.interface_nodes, &stray, &mesh)
            .expect_err("node 0 is not on the interface"),
        InvalidInterface::FacetNodeNotOnInterface {
            facet: 0,
            position: 1,
            node: 0,
        }
    );

    let duplicated = vec![fixture.interface_nodes[0], fixture.interface_nodes[0]];
    assert_eq!(
        StructuralInterface::try_new(&duplicated, &fixture.interface_facets, &mesh)
            .expect_err("the exchange ordering repeats a node"),
        InvalidInterface::DuplicateNode {
            position: 1,
            node: fixture.interface_nodes[0],
        }
    );
}

#[test]
fn a_misshaped_exchange_is_rejected() {
    let fixture = Fixture::new(4, 2);
    let mesh = SimplexMesh::try_new(&fixture.nodes, &fixture.cells).expect("valid mesh");
    let interface =
        StructuralInterface::try_new(&fixture.interface_nodes, &fixture.interface_facets, &mesh)
            .expect("conforming");
    let mut facets = vec![TractionFacet::new([0_usize; 2], [0.0_f64; 2]); 2];
    assert!(interface.read_traction(&[0.0; 3], &mut facets).is_err());
}

#[test]
fn a_uniform_normal_traction_extends_the_strip() {
    // Direction and magnitude together: pulling the free face outward must
    // move it outward, and the axial extension must match the closed form for
    // a bar in plane strain to within the discretisation.
    let fixture = Fixture::new(8, 4);
    let pressure = 1.0e6_f64;
    let (partition, state, _) = couple(&fixture, [pressure, 0.0]);
    let mut exported = vec![0.0_f64; partition.output_dimension()];
    partition
        .export(&state, &mut exported)
        .expect("well-shaped");

    let mean_extension: f64 = exported
        .as_chunks::<2>()
        .0
        .iter()
        .map(|u| u[0])
        .sum::<f64>()
        / (fixture.interface_nodes.len() as f64);
    assert!(
        mean_extension > 0.0,
        "an outward traction pulled the face inward: {mean_extension:.6e}"
    );

    // Uniaxial extension in plane strain: eps_xx = sigma (1 - nu^2) / E for a
    // laterally free bar, so delta = sigma L (1 - nu^2) / E. The strip is
    // laterally constrained only at the clamp, so this is approached rather
    // than matched; the assertion is order of magnitude, with the exact
    // comparison left to the Ares-side oracles that own it.
    let (young, poisson, length) = (200e9_f64, 0.3_f64, 2.0_f64);
    let closed_form = pressure * length * (1.0 - poisson * poisson) / young;
    assert!(
        mean_extension > closed_form * 0.5 && mean_extension < closed_form * 1.5,
        "extension {mean_extension:.6e} is far from the bar closed form {closed_form:.6e}"
    );
}
