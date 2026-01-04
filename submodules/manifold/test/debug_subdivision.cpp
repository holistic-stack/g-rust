// Debug test for subdivision algorithm
// Compares C++ vs Rust implementation

#include "manifold/manifold.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace manifold;

void testSphere() {
    std::cout << "\n=== Test 2: Sphere Subdivision (32 segments, n=8) ===" << std::endl;
    
    // Radius 0.8, 32 segments -> n=8
    Manifold sphere = Manifold::Sphere(0.8, 32);
    
    auto mesh = sphere.GetMeshGL();
    std::cout << "NumVert: " << mesh.NumVert() << std::endl;
    std::cout << "NumTri: " << mesh.NumTri() << std::endl;
    
    // Print bounding box
    auto bbox = sphere.BoundingBox();
    std::cout << "BBox: ("
              << std::fixed << std::setprecision(3) 
              << bbox.min.x << ", " << bbox.min.y << ", " << bbox.min.z << ") - ("
              << bbox.max.x << ", " << bbox.max.y << ", " << bbox.max.z << ")" << std::endl;
    
    // Print vertex positions to see Z distribution
    auto& vertPos = mesh.vertProperties;
    int pos_z = 0, neg_z = 0, mid_z = 0;
    std::cout << "  Z distribution (first 30 vertices):" << std::endl;
    for (int i = 0; i < 30 && i < vertPos.size() / 3; i++) {
        double z = vertPos[i * 3 + 2];
        if (z > 0.1) pos_z++;
        else if (z < -0.1) neg_z++;
        else mid_z++;
        
        if (i < 15) {
            std::cout << "    Vert " << std::setw(2) << i << ": z=" 
                      << std::fixed << std::setprecision(3) << z << std::endl;
        }
    }
    std::cout << "  Summary: +" << pos_z << ", -" << neg_z 
              << ", mid " << mid_z << std::endl;
    
    // Print first 15 triangles and their vertex indices
    auto& triVerts = mesh.triVerts;
    std::cout << "  First 10 triangles:" << std::endl;
    for (int i = 0; i < 10 && i < triVerts.size() / 3; i++) {
        std::cout << "    Tri " << std::setw(2) << i << ": (" 
                  << std::setw(4) << triVerts[3*i] << ", "
                  << std::setw(4) << triVerts[3*i+1] << ", "
                  << std::setw(4) << triVerts[3*i+2] << ")" << std::endl;
    }
}

void testSphere16Segments() {
    std::cout << "\n=== Test 3: Sphere Subdivision (16 segments, n=4) ===" << std::endl;
    
    Manifold sphere = Manifold::Sphere(0.8, 16);
    
    auto mesh = sphere.GetMeshGL();
    std::cout << "NumVert: " << mesh.NumVert() << std::endl;
    std::cout << "NumTri: " << mesh.NumTri() << std::endl;

    // Print first 10 triangles to compare with Rust
    auto& triVerts = mesh.triVerts;
    std::cout << "  First 10 triangles:" << std::endl;
    for (int i = 0; i < 10 && i < triVerts.size() / 3; i++) {
        std::cout << "    Tri " << std::setw(2) << i << ": (" 
                  << std::setw(4) << triVerts[3*i] << ", "
                  << std::setw(4) << triVerts[3*i+1] << ", "
                  << std::setw(4) << triVerts[3*i+2] << ")" << std::endl;
    }
}

int main() {
    std::cout << std::setprecision(10) << std::fixed;
    
    testSphere();
    testSphere16Segments();
    
    return 0;
}
