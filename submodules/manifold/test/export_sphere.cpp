// Export sphere vertex positions for comparison with Rust
#include <iostream>
#include <iomanip>
#include "manifold/manifold.h"

using namespace manifold;

int main() {
    std::cout << std::setprecision(10) << std::fixed;

    std::cout << "\n=== C++ Sphere Export (32 segments) ===" << std::endl;

    Manifold sphere = Manifold::Sphere(0.8, 32);

    std::cout << "NumVert: " << sphere.NumVert() << std::endl;
    std::cout << "NumTri: " << sphere.NumTri() << std::endl;

    auto mesh = sphere.GetMeshGL();

    std::cout << "\nVertex positions:" << std::endl;

    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    double min_z = std::numeric_limits<double>::max();
    double max_z = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < mesh.vertPos.size() / 3; i++) {
        double x = mesh.vertPos[3 * i + 0];
        double y = mesh.vertPos[3 * i + 1];
        double z = mesh.vertPos[3 * i + 2];

        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
        min_z = std::min(min_z, z);
        max_z = std::max(max_z, z);
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
                  << ": (" << std::setprecision(6)
                  << std::setw(12) << x << ", "
                  << std::setw(12) << y << ", "
                  << std::setw(12) << z << ")" << std::endl;
        sampled++;
    }

    return 0;
}
