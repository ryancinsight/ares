//! A validated traction boundary retains its immutable coordinate domain.

use ares::{InvalidBoundary, MisshapedLoad, TractionBoundary, TractionFacet};
use eunomia::{NumericElement, RealField};

fn verify<T: RealField>() {
    let zero = <T as NumericElement>::ZERO;
    let one = <T as NumericElement>::ONE;
    let two = one + one;
    let three = two + one;
    let four = two + two;
    let five = four + one;
    let six = three + three;
    let seven = six + one;
    let nine = seven + two;

    // One unit edge contributes traction / 2 to each endpoint. A foreign
    // coordinate slice is unrepresentable in this API (Rustdoc compile-fail).
    let nodes = [[zero, zero], [one, zero]];
    let facets = [TractionFacet::new([0, 1], [two, four])];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("unit edge");
    for found in [0, 2, 3, 5, 6] {
        let mut loads = vec![one; found];
        assert_eq!(
            boundary.add_consistent_loads(&mut loads),
            Err(MisshapedLoad { expected: 4, found })
        );
        assert_eq!(loads, vec![one; found], "shape rejection writes nothing");
    }
    let mut loads = [one; 4];
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(loads, [two, three, two, three]);

    // Geometry changes require a fresh validation, which rejects this collapse.
    let collapsed = [[zero, zero]; 2];
    assert_eq!(
        TractionBoundary::try_new(&facets, &collapsed).expect_err("collapsed edge"),
        InvalidBoundary::DegenerateFacet { facet: 0 }
    );
    assert_eq!(
        TractionBoundary::try_new(&facets, &nodes[..1]).expect_err("missing endpoint"),
        InvalidBoundary::NodeOutOfRange {
            facet: 0,
            position: 1,
            node: 1,
            nodes: 1
        }
    );

    // Two adjacent unit edges share a node. Both their contributions and the
    // existing load survive; repeating assembly adds the same physical load.
    let joined = [[zero, zero], [one, zero], [two, zero]];
    let facets = [
        TractionFacet::new([0, 1], [two, four]),
        TractionFacet::new([1, 2], [two, four]),
    ];
    let boundary = TractionBoundary::try_new(&facets, &joined).expect("joined edges");
    let mut loads = [one; 6];
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(loads, [two, three, three, five, two, three]);
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(loads, [three, five, five, nine, three, five]);

    // A 3-by-4 right triangle has area 6, so each of its three nodes gets
    // twice the traction. A fourth, interior node receives no surface load.
    let nodes = [
        [zero, zero, zero],
        [three, zero, zero],
        [zero, four, zero],
        [one, one, one],
    ];
    let facets = [TractionFacet::new([0, 1, 2], [one, two, three])];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("triangle area six");
    let mut loads = [one; 12];
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(
        loads,
        [
            three, five, seven, three, five, seven, three, five, seven, one, one, one
        ]
    );

    positive_dimensions::<T>();
}

fn positive_dimensions<T: RealField>() {
    let zero = <T as NumericElement>::ZERO;
    let one = <T as NumericElement>::ONE;
    let two = one + one;
    let three = two + one;
    let four = two + two;
    let five = four + one;

    // D = 1 is the positive-dimension endpoint: a point has unit 0-measure.
    let nodes = [[zero], [one]];
    let facets = [TractionFacet::new([1], [two])];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("one point");
    let mut loads = [one; 2];
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(loads, [one, three]);

    // In four dimensions this tetrahedral facet has 3-measure 2*3*4/6=4,
    // hence each of its four nodes receives exactly the applied traction.
    let nodes = [
        [zero; 4],
        [two, zero, zero, zero],
        [zero, three, zero, zero],
        [zero, zero, four, zero],
    ];
    let facets = [TractionFacet::new([0, 1, 2, 3], [one, two, three, four])];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("4-D facet");
    let mut loads = [one; 16];
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(
        loads,
        [
            two, three, four, five, two, three, four, five, two, three, four, five, two, three,
            four, five
        ]
    );

    let facets: [TractionFacet<T, 2>; 0] = [];
    let nodes = [[zero, zero]];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("no applied traction");
    let mut loads = [one; 2];
    boundary
        .add_consistent_loads(&mut loads)
        .expect("bound shape");
    assert_eq!(loads, [one; 2]);
    let boundary = TractionBoundary::try_new(&facets, &[]).expect("empty coordinate domain");
    boundary.add_consistent_loads(&mut []).expect("empty shape");
}

#[test]
fn bound_geometry_binary32() {
    verify::<f32>();
}
#[test]
fn bound_geometry_binary64() {
    verify::<f64>();
}
#[test]
fn bound_geometry_wrapped_binary32() {
    verify::<eunomia::F32>();
}
#[test]
fn bound_geometry_wrapped_binary64() {
    verify::<eunomia::F64>();
}
