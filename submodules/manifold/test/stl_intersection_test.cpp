// Test for csg_intersection_centered to verify expected triangle count
#include "manifold.h"
#include "test.h"
#include <iostream>

TEST(STLDebug, CsgIntersectionCentered) {
    // Exact parameters from generate_reference.cpp line 131-135
    Manifold cube = Manifold::Cube(vec3(10.0), true);
    Manifold sphere = Manifold::Sphere(7.0, 32);
    Manifold result = cube ^ sphere;  // Intersection operator

    auto mesh = result.GetMesh();
    int num_tri = mesh.NumTri();

    std::cout << "=== CSG Intersection Centered ===" << std::endl;
    std::cout << "Cube: vec3(10.0) with forceBounds=true" << std::endl;
    std::cout << "Sphere: radius=7.0, segments=32" << std::endl;
    std::cout << "Result triangles: " << num_tri << std::endl;
    std::cout << "Expected (from reference): 500" << std::endl;
    std::cout << "Match: " << (num_tri == 500 ? "YES" : "NO") << std::endl;

    // The reference produces 500 triangles
    EXPECT_EQ(num_tri, 500);
}

TEST(STLDebug, CsgIntersectionPartial) {
    // Exact parameters from generate_reference.cpp line 150-154
    Manifold cube = Manifold::Cube(vec3(8.0), true);
    Manifold sphere = Manifold::Sphere(6.0, 32).Translate(vec3(3.0, 0.0, 0.0));
    Manifold result = cube ^ sphere;

    auto mesh = result.GetMesh();
    int num_tri = mesh.NumTri();

    std::cout << "=== CSG Intersection Partial ===" << std::endl;
    std::cout << "Cube: vec3(8.0) with forceBounds=true" << std::endl;
    std::cout << "Sphere: radius=6.0, segments=32, translated by vec3(3.0, 0.0, 0.0)" << std::endl;
    std::cout << "Result triangles: " << num_tri << std::endl;
    std::cout << "Expected (from reference): 214" << std::endl;
    std::cout << "Match: " << (num_tri == 214 ? "YES" : "NO") << std::endl;

    EXPECT_EQ(num_tri, 214);
}
