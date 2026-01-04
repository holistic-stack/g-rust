// C++ test to debug boolean operations
#include <iostream>
#include <iomanip>
#include "manifold/manifold.h"

using namespace manifold;

void print_mesh_info(const std::string& name, const Manifold& m) {
    std::cout << name << ":" << std::endl;
    std::cout << "  NumTri: " << m.NumTri() << std::endl;
    std::cout << "  NumVert: " << m.NumVert() << std::endl;
    std::cout << "  Volume: " << std::fixed << std::setprecision(3) << m.Volume() << std::endl;
    std::cout << "  IsEmpty: " << (m.IsEmpty() ? "true" : "false") << std::endl;
    std::cout << "  NumDegenerateTris: " << m.NumDegenerateTris() << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << "=== Boolean Operations Debug Test ===" << std::endl << std::endl;

    // Test 1: Cube - Sphere Difference
    {
        std::cout << "Test 1: Cube - Sphere Difference" << std::endl;

        // Create cube centered at origin with size 2.0
        Manifold cube = Manifold::Cube({2.0, 2.0, 2.0}, true);

        // Create sphere with radius 0.8
        Manifold sphere = Manifold::Sphere(0.8, 32);

        print_mesh_info("Cube", cube);
        print_mesh_info("Sphere", sphere);

        // Perform difference operation (cube - sphere)
        Manifold result = cube - sphere;

        print_mesh_info("Difference", result);

        // Expected volume: 8.0 - 2.14 ≈ 5.86
        double expected_volume = 5.86;
        double result_volume = result.Volume();
        std::cout << "Expected volume: " << expected_volume << std::endl;
        std::cout << "Result volume: " << result_volume << std::endl;

        if (result_volume > 5.0 && result_volume < 7.0) {
            std::cout << "PASS: Volume in expected range" << std::endl;
        } else {
            std::cout << "FAIL: Volume out of expected range" << std::endl;
        }
        std::cout << std::endl;
    }

    // Test 2: Cube ^ Sphere Intersection
    {
        std::cout << "Test 2: Cube ^ Sphere Intersection" << std::endl;

        // Create cube centered at origin with size 2.0
        Manifold cube = Manifold::Cube({2.0, 2.0, 2.0}, true);

        // Create sphere with radius 0.8
        Manifold sphere = Manifold::Sphere(0.8, 32);

        // Perform intersection operation (cube ^ sphere)
        Manifold result = cube ^ sphere;

        print_mesh_info("Intersection", result);

        // Expected: Not empty, sphere is inside cube
        if (result.NumTri() > 0 && result.NumVert() > 0) {
            std::cout << "PASS: Intersection is not empty" << std::endl;
        } else {
            std::cout << "FAIL: Intersection is empty" << std::endl;
        }
        std::cout << std::endl;
    }

    // Test 3: Cube + Sphere Union
    {
        std::cout << "Test 3: Cube + Sphere Union" << std::endl;

        // Create cube centered at origin with size 2.0
        Manifold cube = Manifold::Cube({2.0, 2.0, 2.0}, true);

        // Create sphere with radius 0.8
        Manifold sphere = Manifold::Sphere(0.8, 32);

        // Perform union operation (cube + sphere)
        Manifold result = cube + sphere;

        print_mesh_info("Union", result);

        // Expected: Union should preserve both volumes (approximately)
        double cube_volume = cube.Volume();
        double sphere_volume = sphere.Volume();
        double union_volume = result.Volume();
        std::cout << "Cube volume: " << cube_volume << std::endl;
        std::cout << "Sphere volume: " << sphere_volume << std::endl;
        std::cout << "Union volume: " << union_volume << std::endl;

        if (union_volume > cube_volume && union_volume < cube_volume + sphere_volume + 0.1) {
            std::cout << "PASS: Union volume reasonable" << std::endl;
        } else {
            std::cout << "FAIL: Union volume unreasonable" << std::endl;
        }
        std::cout << std::endl;
    }

    // Test 4: Coplanar cubes (like test_boolean_operators)
    {
        std::cout << "Test 4: Coplanar Cubes Union/Diff/Intersect" << std::endl;

        Manifold cube1 = Manifold::Cube({1.0, 1.0, 1.0}, false).Translate({0.0, 0.0, 0.0});
        Manifold cube2 = Manifold::Cube({1.0, 1.0, 1.0}, false).Translate({0.5, 0.5, 0.5});

        print_mesh_info("Cube1", cube1);
        print_mesh_info("Cube2", cube2);

        Manifold union_result = cube1 + cube2;
        Manifold diff_result = cube1 - cube2;
        Manifold inter_result = cube1 ^ cube2;

        print_mesh_info("Union", union_result);
        print_mesh_info("Difference", diff_result);
        print_mesh_info("Intersection", inter_result);

        std::cout << std::endl;
    }

    // Test 5: Tangent face 50% overlap along X
    {
        std::cout << "Test 5: Tangent Face 50% Overlap (x=0.5)" << std::endl;

        Manifold cube1 = Manifold::Cube({1.0, 1.0, 1.0}, false);
        Manifold cube2 = Manifold::Cube({1.0, 1.0, 1.0}, false).Translate({0.5, 0.0, 0.0});

        Manifold union_result = cube1 + cube2;
        Manifold inter_result = cube1 ^ cube2;

        print_mesh_info("Cube1", cube1);
        print_mesh_info("Cube2", cube2);
        print_mesh_info("Union", union_result);
        print_mesh_info("Intersection", inter_result);

        const double union_vol = union_result.Volume();
        const double inter_vol = inter_result.Volume();

        std::cout << std::setprecision(6)
                  << "Expected union ~1.5, got " << union_vol << std::endl
                  << "Expected intersection ~0.5, got " << inter_vol << std::endl
                  << std::endl;
    }

    std::cout << "=== All Tests Complete ===" << std::endl;
    return 0;
}
