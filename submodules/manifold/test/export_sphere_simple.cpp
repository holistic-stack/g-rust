// Simple C++ program to export sphere vertex positions
#include "manifold/manifold.h"
#include <iostream>
#include <iomanip>

using namespace manifold;

int main() {
    std::cout << std::setprecision(10) << std::fixed;

    // Create sphere (matching Rust test parameters)
    std::cout << "\n=== C++ Sphere Export (32 segments, radius 0.8) ===" << std::endl;

    Manifold sphere = Manifold::Sphere(0.8, 32);

    std::cout << "NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "NumTri: " << sphere.NumTri() << std::endl;

    auto mesh = sphere.GetMeshGL();

    // Access vertex positions directly from vertProperties array
    // MeshGLP structure: prop1, prop2, prop3, prop4 for each vertex
    std::cout << "\nVertex positions (from vertProperties):" << std::endl;

    double min_x = 1e100;
    double max_x = -1e100;
    double min_y = 1e100;
    double max_y = -1e100;
    double min_z = 1e100;
    double max_z = -1e100;

    for (size_t i = 0; i < mesh.vertProperties.size() / 4; i++) {
        double x = mesh.vertProperties[i * 4 + 0];
        double y = mesh.vertProperties[i * 4 + 1];
        double z = mesh.vertProperties[i * 4 + 2];

        min_x = fmin(min_x, x);
        max_x = fmax(max_x, x);
        min_y = fmin(min_y, y);
        max_y = fmax(max_y, y);
        min_z = fmin(min_z, z);
        max_z = fmax(max_z, z);

        if (i < 15) {  // Print first 15 vertices
            std::cout << i
                      << ": (" << std::setprecision(6) << x << ", "
                      << y << ", " << z << ")" << std::endl;
        }
    }

    std::cout << "\nVertex statistics:" << std::endl;
    std::cout << "  X: [" << min_x << ", " << max_x << "]" << std::endl;
    std::cout << "  Y: [" << min_y << ", " << max_y << "]" << std::endl;
    std::cout << "  Z: [" << min_z << ", " << max_z << "]" << std::endl;

    int pos_z_count = 0;
    int neg_z_count = 0;
    int mid_z_count = 0;

    for (size_t i = 0; i < mesh.vertProperties.size() / 4; i++) {
        double z = mesh.vertProperties[i * 4 + 2];
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

    return 0;
}
