use glam::DVec3;
use glam::DVec3;
use manifold::Manifold;
use manifold_types::OpType;

use crate::tests::util::{assert_deterministic_runs, canonical_signature};

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
fn boolean_tetra_basic() {
    let tetra = Manifold::tetrahedron();
    let tetra2 = translate(&tetra, DVec3::splat(0.5));
    let result = tetra2 - tetra;

    let (verts, tris) = canonical_signature(&result);
    assert_eq!(verts.len(), 8);
    assert_eq!(tris.len(), 12);

    assert_deterministic_runs(3, || result.clone());
}

#[test]
fn boolean_empty_original() {
    let cube = Manifold::cube(DVec3::splat(1.0), false);
    let tet = Manifold::tetrahedron();
    let far = translate(&cube, DVec3::new(3.0, 4.0, 5.0));
    let result = tet - far;

    let (verts, tris) = canonical_signature(&result);
    assert!(verts.is_empty());
    assert!(tris.is_empty());
    assert_eq!(result.volume(), 0.0);
    assert_eq!(result.surface_area(), 0.0);

    assert_deterministic_runs(3, || result.clone());
}

#[test]
fn boolean_non_intersecting() {
    let a = Manifold::cube(DVec3::splat(1.0), false);
    let b = translate(&a, DVec3::new(2.0, 2.0, 2.0));
    let result = a ^ b;

    let (verts, tris) = canonical_signature(&result);
    assert!(verts.is_empty());
    assert!(tris.is_empty());
    assert_eq!(result.volume(), 0.0);
    assert_eq!(result.surface_area(), 0.0);

    assert_deterministic_runs(3, || result.clone());
}

#[test]
fn boolean_self_subtract_empty() {
    let cube = Manifold::cube(DVec3::splat(1.0), false);
    let result = cube.clone() - cube;

    let (verts, tris) = canonical_signature(&result);
    assert!(verts.is_empty());
    assert!(tris.is_empty());
    assert_eq!(result.volume(), 0.0);
    assert_eq!(result.surface_area(), 0.0);

    assert_deterministic_runs(3, || result.clone());
}
