// Copyright 2025 Debug Test
//
// Debug test for csg_difference_hole CSG operation
// Creates a cube and cylinder, then performs difference
// Outputs intermediate pipeline states for comparison with Rust
//

#include "manifold/manifold.h"
#include <iostream>
#include <iomanip>

using namespace manifold;

int main() {
    std::cout << "\n=== C++ Test: csg_difference_hole ===" << std::endl;

    // Create cube (size 10, centered at origin)
    // Matches Rust: Manifold::cube(DVec3::splat(10.0), true)
    Manifold cube = Manifold::Cube(vec3(10.0, 10.0, 10.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;

    // Create cylinder (radius=20, height=2, segments=32)
    // Matches Rust: Manifold::cylinder(20.0, 2.0, 2.0, 32, true)
    Manifold cylinder = Manifold::Cylinder(20.0, 2.0, 32);
    std::cout << "\nCylinder:" << std::endl;
    std::cout << "  NumVert: " << cylinder.NumVert() << std::endl;
    std::cout << "  NumTri: " << cylinder.NumTri() << std::endl;

    // Perform difference operation
    // Matches Rust: cube.difference(&cylinder)
    std::cout << "\nPerforming Difference: cube - cylinder" << std::endl;
    Manifold result = cube - cylinder;

    // Output final results
    std::cout << "\n=== Final Result ===" << std::endl;
    std::cout << "Result:" << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;
    std::cout << "  Status: " << (result.Status() == Manifold::Error::NoError ? "Valid" : "Invalid") << std::endl;

    // Expected: 144 triangles
    std::cout << "\nExpected: 144 triangles" << std::endl;

    return 0;
}
