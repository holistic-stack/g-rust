// Copyright 2025 Debug Test
//
// Debug test for CSG boolean operations to compare C++ vs Rust
//

#include "manifold/manifold.h"
#include <iostream>
#include <iomanip>

using namespace manifold;

// Helper to print bbox
void printBBox(const std::string& name, const Box& box) {
    std::cout << "  " << name << ": ("
              << std::fixed << std::setprecision(3) << box.min.x << ", "
              << box.min.y << ", "
              << box.min.z << ") - ("
              << box.max.x << ", "
              << box.max.y << ", "
              << box.max.z << ")" << std::endl;
}

// Helper to export sphere vertex positions
void exportSphereVertices(const Manifold& sphere) {
    std::cout << "\n=== Exporting Sphere Vertices ===\n" << std::endl;

    auto mesh = sphere.GetMeshGL();
    std::cout << "NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "NumTri: " << sphere.NumTri() << std::endl;

    std::cout << "\nVertex positions:" << std::endl;

    double min_x = 1e100;
    double max_x = -1e100;
    double min_y = 1e100;
    double max_y = -1e100;
    double min_z = 1e100;
    double max_z = -1e100;

    for (size_t i = 0; i < mesh.vertPos.size() / 3; i++) {
        double x = mesh.vertPos[3 * i + 0];
        double y = mesh.vertPos[3 * i + 1];
        double z = mesh.vertPos[3 * i + 2];

        min_x = fmin(min_x, x);
        max_x = fmax(max_x, x);
        min_y = fmin(min_y, y);
        max_y = fmax(max_y, y);
        min_z = fmin(min_z, z);
        max_z = fmax(max_z, z);
    }

    std::cout << "\nVertex statistics:" << std::endl;
    std::cout << "  X: [" << min_x << ", " << max_x << "]" << std::endl;
    std::cout << "  Y: [" << min_y << ", " << max_y << "]" << std::endl;
    std::cout << "  Z: [" << min_z << ", " << max_z << "]" << std::endl;

    int pos_z_count = 0;
    int neg_z_count = 0;
    int mid_z_count = 0;

    for (size_t i = 0; i < mesh.vertPos.size() / 3; i++) {
        double z = mesh.vertPos[3 * i + 2];

        if (z > 0.1) {
            pos_z_count++;
        } else if (z < -0.1) {
            neg_z_count++;
        } else {
            mid_z_count++;
        }
    }

    std::cout << "\nZ distribution:" << std::endl;
    std::cout << "  Positive Z (> 0.1): " << pos_z_count << std::endl;
    std::cout << "  Negative Z (< -0.1): " << neg_z_count << std::endl;
    std::cout << "  Mid Z (-0.1 to 0.1): " << mid_z_count << std::endl;

    std::cout << "\nSample vertices:" << std::endl;
    int sampled = 0;
    for (size_t i = 0; i < mesh.vertPos.size() / 3 && sampled < 15; i++) {
        double x = mesh.vertPos[3 * i + 0];
        double y = mesh.vertPos[3 * i + 1];
        double z = mesh.vertPos[3 * i + 2];

        std::cout << "  Vert " << std::setw(3) << i
                  << ": (" << std::fixed << std::setprecision(6)
                  << x << ", " << y << ", " << z << ")" << std::endl;
        sampled++;
    }
}

// Test cube-sphere intersection with detailed output
void testCubeSphereIntersection() {
    std::cout << "\n=== C++ Test: Cube-Sphere Intersection ===\n" << std::endl;

    // Create cube (size 2.0, centered at origin)
    Manifold cube = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;
    printBBox("BBox", cube.BoundingBox());

    // Create sphere (radius 0.8, 32 segments)
    Manifold sphere = Manifold::Sphere(0.8, 32);
    std::cout << "\nSphere:" << std::endl;
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;
    printBBox("BBox", sphere.BoundingBox());

    // Export sphere vertices for comparison
    exportSphereVertices(sphere);

    // Check overlap
    bool overlap = cube.BoundingBox().DoesOverlap(sphere.BoundingBox());
    std::cout << "\nBounding boxes overlap: " << (overlap ? "YES" : "NO") << std::endl;

    // Perform intersection
    std::cout << "\nPerforming intersection..." << std::endl;
    Manifold result = cube ^ sphere;

    std::cout << "\nIntersection result:" << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.Status()) << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;

    if (!result.IsEmpty()) {
        printBBox("BBox", result.BoundingBox());
        double volume = result.Volume();
        std::cout << "  Volume: " << std::fixed << std::setprecision(6) << volume << std::endl;
    } else {
        std::cout << "  Result is EMPTY!" << std::endl;
    }
}

// Test cube-sphere union with detailed output
void testCubeSphereUnion() {
    std::cout << "\n=== C++ Test: Cube-Sphere Union ===\n" << std::endl;

    // Create cube (size 2.0, centered at origin)
    Manifold cube = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;
    printBBox("BBox", cube.BoundingBox());

    // Create sphere (radius 0.8, 32 segments)
    Manifold sphere = Manifold::Sphere(0.8, 32);
    std::cout << "\nSphere:" << std::endl;
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;
    printBBox("BBox", sphere.BoundingBox());

    // Export sphere vertices for comparison
    exportSphereVertices(sphere);

    // Perform union
    std::cout << "\nPerforming union..." << std::endl;
    Manifold result = cube + sphere;

    std::cout << "\nUnion result:" << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.Status()) << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;

    if (!result.IsEmpty()) {
        printBBox("BBox", result.BoundingBox());
        double volume = result.Volume();
        std::cout << "  Volume: " << std::fixed << std::setprecision(6) << volume << std::endl;
    } else {
        std::cout << "  Result is EMPTY!" << std::endl;
    }
}

// Test cube-sphere difference with detailed output
void testCubeSphereDifference() {
    std::cout << "\n=== C++ Test: Cube-Sphere Difference ===\n" << std::endl;

    // Create cube (size 2.0, centered at origin)
    Manifold cube = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;
    printBBox("BBox", cube.BoundingBox());

    // Create sphere (radius 0.8, 32 segments)
    Manifold sphere = Manifold::Sphere(0.8, 32);
    std::cout << "\nSphere:" << std::endl;
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;
    printBBox("BBox", sphere.BoundingBox());

    // Export sphere vertices for comparison
    exportSphereVertices(sphere);

    // Perform difference
    std::cout << "\nPerforming difference..." << std::endl;
    Manifold result = cube - sphere;

    std::cout << "\nDifference result:" << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.Status()) << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;

    if (!result.IsEmpty()) {
        printBBox("BBox", result.BoundingBox());
        double volume = result.Volume();
        std::cout << "  Volume: " << std::fixed << std::setprecision(6) << volume << std::endl;
    } else {
        std::cout << "  Result is EMPTY!" << std::endl;
    }
}

int main() {
    testCubeSphereIntersection();
    testCubeSphereUnion();
    testCubeSphereDifference();

    return 0;
}

// Test cube-sphere intersection with detailed output
void testCubeSphereIntersection() {
    std::cout << "\n=== C++ Test: Cube-Sphere Intersection ===\n" << std::endl;

    // Create cube (size 2.0, centered at origin)
    Manifold cube = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;
    printBBox("BBox", cube.BoundingBox());

    // Create sphere (radius 0.8, 32 segments)
    Manifold sphere = Manifold::Sphere(0.8, 32);
    std::cout << "\nSphere:" << std::endl;
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;
    printBBox("BBox", sphere.BoundingBox());

    // Check overlap
    bool overlap = cube.BoundingBox().DoesOverlap(sphere.BoundingBox());
    std::cout << "\nBounding boxes overlap: " << (overlap ? "YES" : "NO") << std::endl;

    // Perform intersection
    std::cout << "\nPerforming intersection..." << std::endl;
    Manifold result = cube ^ sphere;  // ^ operator is intersection

    std::cout << "\nIntersection result:" << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.Status()) << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;

    if (!result.IsEmpty()) {
        printBBox("BBox", result.BoundingBox());
        double volume = result.Volume();
        std::cout << "  Volume: " << std::fixed << std::setprecision(6) << volume << std::endl;
    } else {
        std::cout << "  Result is EMPTY!" << std::endl;
    }
}

// Test cube-sphere union with detailed output
void testCubeSphereUnion() {
    std::cout << "\n=== C++ Test: Cube-Sphere Union ===\n" << std::endl;

    // Create cube (size 2.0, centered at origin)
    Manifold cube = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;
    printBBox("BBox", cube.BoundingBox());

    // Create sphere (radius 0.8, 32 segments)
    Manifold sphere = Manifold::Sphere(0.8, 32);
    std::cout << "\nSphere:" << std::endl;
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;
    printBBox("BBox", sphere.BoundingBox());

    // Perform union
    std::cout << "\nPerforming union..." << std::endl;
    Manifold result = cube + sphere;

    std::cout << "\nUnion result:" << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.Status()) << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;

    if (!result.IsEmpty()) {
        printBBox("BBox", result.BoundingBox());
        double volume = result.Volume();
        std::cout << "  Volume: " << std::fixed << std::setprecision(6) << volume << std::endl;
    } else {
        std::cout << "  Result is EMPTY!" << std::endl;
    }
}

// Test cube-sphere difference with detailed output
void testCubeSphereDifference() {
    std::cout << "\n=== C++ Test: Cube-Sphere Difference ===\n" << std::endl;

    // Create cube (size 2.0, centered at origin)
    Manifold cube = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
    std::cout << "Cube:" << std::endl;
    std::cout << "  NumVert: " << cube.NumVert() << std::endl;
    std::cout << "  NumTri: " << cube.NumTri() << std::endl;
    printBBox("BBox", cube.BoundingBox());

    // Create sphere (radius 0.8, 32 segments)
    Manifold sphere = Manifold::Sphere(0.8, 32);
    std::cout << "\nSphere:" << std::endl;
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;
    printBBox("BBox", sphere.BoundingBox());

    // Perform difference
    std::cout << "\nPerforming difference..." << std::endl;
    Manifold result = cube - sphere;

    std::cout << "\nDifference result:" << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.Status()) << std::endl;
    std::cout << "  NumVert: " << result.NumVert() << std::endl;
    std::cout << "  NumTri: " << result.NumTri() << std::endl;

    if (!result.IsEmpty()) {
        printBBox("BBox", result.BoundingBox());
        double volume = result.Volume();
        std::cout << "  Volume: " << std::fixed << std::setprecision(6) << volume << std::endl;
    } else {
        std::cout << "  Result is EMPTY!" << std::endl;
    }
}

int main() {
    std::cout << std::setprecision(10);

    testCubeSphereIntersection();
    testCubeSphereUnion();
    testCubeSphereDifference();

    return 0;
}
