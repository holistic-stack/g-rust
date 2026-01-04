// Test sphere subdivision levels
#include "../src/utils.h"
#include "manifold/manifold.h"
#include <iostream>
#include <iomanip>

using namespace manifold;

void testSphere(int segments, const char* name) {
    std::cout << "\n=== Sphere: " << name << " (" << segments << " segments) ===" << std::endl;

    Manifold sphere = Manifold::Sphere(0.8, segments);
    std::cout << "  NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "  NumTri: " << sphere.NumTri() << std::endl;

    auto bbox = sphere.BoundingBox();
    std::cout << "  BBox: ("
              << std::fixed << std::setprecision(3) << bbox.min.x << ", "
              << bbox.min.y << ", "
              << bbox.min.z << ") - ("
              << bbox.max.x << ", "
              << bbox.max.y << ", "
              << bbox.max.z << ")" << std::endl;
}

int main() {
    int n;

    // Test various circular segment values
    testSphere(32, "32 segments");

    // Calculate n = (circularSegments + 3) / 4
    n = (32 + 3) / 4;
    std::cout << "\nFor 32 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    testSphere(28, "28 segments");
    n = (28 + 3) / 4;
    std::cout << "\nFor 28 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    testSphere(24, "24 segments");
    n = (24 + 3) / 4;
    std::cout << "\nFor 24 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    testSphere(16, "16 segments");
    n = (16 + 3) / 4;
    std::cout << "\nFor 16 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    testSphere(12, "12 segments");
    n = (12 + 3) / 4;
    std::cout << "\nFor 12 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    testSphere(8, "8 segments");
    n = (8 + 3) / 4;
    std::cout << "\nFor 8 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    testSphere(4, "4 segments");
    n = (4 + 3) / 4;
    std::cout << "\nFor 4 segments: n = " << n << ", subdivisions = n - 1 = " << n - 1 << std::endl;

    return 0;
}
