// Copyright 2025 Debug Test
//
// Export sphere vertex positions for comparison with Rust
//

#include "manifold/manifold.h"
#include <iostream>
#include <cmath>

using namespace manifold;

int main() {
    std::cout << std::setprecision(10) << std::fixed;

    std::cout << "\n=== C++ Sphere Export (32 segments) ===" << std::endl;

    Manifold sphere = Manifold::Sphere(0.8, 32);

    std::cout << "NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "NumTri: " << sphere.NumTri() << std::endl;

    auto mesh = sphere.GetMeshGL();

    std::cout << "\nVertex positions:" << std::endl;

    double min_x = 1e100;
    double max_x = -1e100;
    double min_y = 1e100;
    double max_y = -1e100;
    double min_z = 1e100;
    double max_z = -1e100;

    for (size_t i = 0; i < mesh.vertProperties.size() / 3; i++) {
        double x = mesh.vertProperties[i * 3 + 0];
        double y = mesh.vertProperties[i * 3 + 1];
        double z = mesh.vertProperties[i * 3 + 2];

        min_x = std::fmin(min_x, x);
        max_x = std::fmax(max_x, x);
        min_y = std::fmin(min_y, y);
        max_y = std::fmax(max_y, y);
        min_z = std::fmin(min_z, z);
        max_z = std::fmax(max_z, z);
    }

    std::cout << "\nVertex statistics:" << std::endl;
    std::cout << "  X: [" << min_x << ", " << max_x << "]" << std::endl;
    std::cout << "  Y: [" << min_y << ", " << max_y << "]" << std::endl;
    std::cout << "  Z: [" << min_z << ", " << max_z << "]" << std::endl;

    int pos_z_count = 0;
    int neg_z_count = 0;
    int mid_z_count = 0;

    for (size_t i = 0; i < mesh.vertProperties.size() / 3; i++) {
        double z = mesh.vertProperties[i * 3 + 2];

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
    for (size_t i = 0; i < mesh.vertProperties.size() / 3 && sampled < 15; i++) {
        double x = mesh.vertProperties[i * 3 + 0];
        double y = mesh.vertProperties[i * 3 + 1];
        double z = mesh.vertProperties[i * 3 + 2];

        std::cout << "  Vert " << std::setw(3) << i
                  << ": (" << std::fixed << std::setprecision(6)
                  << x << ", " << y << ", " << z << ")" << std::endl;
        sampled++;
    }

    return 0;
}
