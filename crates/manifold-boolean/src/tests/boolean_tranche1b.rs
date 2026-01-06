use glam::DVec3;
use manifold::Manifold;
use manifold_types::OpType;

use crate::tests::util::{assert_deterministic_runs, counts};

fn translate(m: &Manifold, v: DVec3) -> Manifold {
    m.translate(v)
}

fn boolean(m1: &Manifold, m2: &Manifold, op: OpType) -> Manifold {
    let op_fn = match op {
        OpType::Add => Manifold::add,
        OpType::Subtract => Manifold::sub,
        OpType::Intersect => Manifold::bitxor,
    };
    op_fn(m1.clone(), m2.clone())
}

#[test]
fn boolean_cubes_basic() {
    let mut result =
        Manifold::cube(DVec3::new(1.2, 1.0, 1.0), true).translate(DVec3::new(0.0, -0.5, 0.5));
    result = result
        + Manifold::cube(DVec3::new(1.0, 0.8, 0.5), false).translate(DVec3::new(-0.5, 0.0, 0.5));
    result = result
        + Manifold::cube(DVec3::new(1.2, 0.1, 0.5), false).translate(DVec3::new(-0.6, -0.1, 0.0));

    let (v, t) = counts(&result);
    assert!(v > 0);
    assert!(t > 0);
    assert_relative_eq!(result.volume(), 1.6, 0.001);
    assert_relative_eq!(result.surface_area(), 9.2, 0.01);

    assert_deterministic_runs(3, || result.clone());
}

#[test]
fn boolean_simple_cube_regression_like() {
    let a = Manifold::cube(DVec3::splat(1.0), false);
    let b = translate(&a, DVec3::new(1.0, 1.0, 0.0));
    let result = a + b;

    let (v, t) = counts(&result);
    assert!(v > 0);
    assert!(t > 0);

    assert_deterministic_runs(3, || result.clone());
}

#[test]
fn boolean_union_difference_basic() {
    let block = Manifold::cube(DVec3::splat(1.0), true) - Manifold::sphere(1.0, 0);
    let result = block.clone() + translate(&block, DVec3::new(0.0, 0.0, 1.0));

    let (v, t) = counts(&result);
    assert!(v > 0);
    assert!(t > 0);
    assert_relative_eq!(result.volume(), block.volume() * 2.0, 0.0001);

    assert_deterministic_runs(3, || result.clone());
}

#[test]
fn boolean_batch_boolean_like() {
    let base = Manifold::cube(DVec3::splat(1.0), false);
    let b1 = translate(&base, DVec3::new(0.5, 0.0, 0.0));
    let b2 = translate(&base, DVec3::new(0.0, 0.5, 0.0));
    let b3 = translate(&base, DVec3::new(0.0, 0.0, 0.5));

    let result = base + b1 + b2 + b3;
    let (v, t) = counts(&result);
    assert_eq!(v, 152);
    assert_eq!(t, 300);

    assert_deterministic_runs(3, || result.clone());
}
