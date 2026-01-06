use manifold::Manifold;

#[test]
fn boolean_tetra_difference_is_nonempty() {
    let tetra = Manifold::tetrahedron();
    let tetra2 = tetra.translate(glam::DVec3::splat(0.5));

    let result = tetra2 - tetra;

    assert!(result.num_vert() > 0);
    assert!(result.num_tri() > 0);
}
