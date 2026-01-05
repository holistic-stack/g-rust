use manifold_boolean::kernel::boolean3::*;

#[test]
fn test_boolean_kernel_compilation() {
    let mesh_p = manifold_boolean::kernel::mod_rs::ManifoldImpl::new();
    let mesh_q = manifold_boolean::kernel::mod_rs::ManifoldImpl::new();

    // Test that Boolean3 constructor works
    let result = boolean3_new(&mesh_p, &mesh_q, manifold_types::OpType::Union);
    assert!(result.valid);
}
