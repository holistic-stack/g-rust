// Copyright 2024 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Generate reference STL files for comparison with Rust port.
// This program generates test geometries from simple to complex
// and exports them as ASCII STL files using Assimp via meshIO.

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "manifold/manifold.h"
#include "manifold/meshIO.h"

using namespace manifold;
namespace fs = std::filesystem;

// Forward declaration for MengerSponge (from samples)
namespace {
void Fractal(std::vector<Manifold>& holes, Manifold& hole, double w,
             vec2 position, int depth, int maxDepth) {
  w /= 3;
  holes.push_back(hole.Scale({w, w, 1.0}).Translate(vec3(position, 0.0)));
  if (depth == maxDepth) return;

  vec2 offsets[8] = {{-w, -w}, {-w, 0.0}, {-w, w}, {0.0, w},
                     {w, w},   {w, 0.0},  {w, -w}, {0.0, -w}};
  for (int i = 0; i < 8; ++i) {
    Fractal(holes, hole, w, position + offsets[i], depth + 1, maxDepth);
  }
}

Manifold MengerSponge(int n) {
  Manifold result = Manifold::Cube(vec3(1.0), true);
  std::vector<Manifold> holes;
  Fractal(holes, result, 1.0, {0.0, 0.0}, 1, n);
  Manifold hole = Manifold::Compose(holes);

  result -= hole;
  hole = hole.Rotate(90);
  result -= hole;
  hole = hole.Rotate(0, 0, 90);
  result -= hole;

  return result;
}
}  // namespace

// Export helper - uses Assimp for ASCII STL
void ExportSTL(const Manifold& m, const std::string& output_dir,
               const std::string& case_id) {
  std::string path = output_dir + "/" + case_id + ".stl";
  ExportMesh(path, m.GetMeshGL(), {});
  std::cout << "  Exported: " << path << " (" << m.NumTri() << " triangles)"
            << std::endl;
}

// Test case generator
struct TestCase {
  std::string id;
  std::function<Manifold()> generate;
};

std::vector<TestCase> GetTestCases() {
  std::vector<TestCase> cases;

  // ============== TIER 1: Basic Primitives ==============
  cases.push_back({"primitive_tetrahedron", []() {
                     return Manifold::Tetrahedron();
                   }});

  cases.push_back({"primitive_cube", []() {
                     return Manifold::Cube(vec3(10.0), true);
                   }});

  cases.push_back({"primitive_sphere_lo", []() {
                     return Manifold::Sphere(5.0, 16);
                   }});

  cases.push_back({"primitive_cylinder", []() {
                     return Manifold::Cylinder(10.0, 3.0, 3.0, 16, true);
                   }});

  // ============== TIER 2: High-Res & Transforms ==============
  cases.push_back({"primitive_sphere_hi", []() {
                     return Manifold::Sphere(5.0, 64);
                   }});

  cases.push_back({"transform_translate", []() {
                     return Manifold::Cube(vec3(10.0), true)
                         .Translate(vec3(5.0, 5.0, 5.0));
                   }});

  cases.push_back({"transform_rotate", []() {
                     return Manifold::Cube(vec3(10.0), true).Rotate(45.0, 0.0, 0.0);
                   }});

  cases.push_back({"transform_scale", []() {
                     return Manifold::Sphere(5.0, 32).Scale(vec3(2.0, 1.0, 1.0));
                   }});

  cases.push_back({"transform_mirror", []() {
                     return Manifold::Cube(vec3(10.0), true).Mirror(vec3(1.0, 0.0, 0.0));
                   }});

  // ============== TIER 3: Simple CSG Operations ==============
  cases.push_back({"csg_union_centered", []() {
                     auto cube = Manifold::Cube(vec3(8.0), true);
                     auto sphere = Manifold::Sphere(5.0, 32);
                     return cube + sphere;
                   }});

  cases.push_back({"csg_difference_centered", []() {
                     auto cube = Manifold::Cube(vec3(10.0), true);
                     auto sphere = Manifold::Sphere(6.0, 32);
                     return cube - sphere;
                   }});

  cases.push_back({"csg_intersection_centered", []() {
                     auto cube = Manifold::Cube(vec3(10.0), true);
                     auto sphere = Manifold::Sphere(7.0, 32);
                     return cube ^ sphere;
                   }});

  // ============== TIER 4: Offset CSG & Edge Cases ==============
  cases.push_back({"csg_union_offset", []() {
                     auto sphere = Manifold::Sphere(5.0, 32);
                     auto cube = Manifold::Cube(vec3(8.0), true).Translate(vec3(5.0, 0.0, 0.0));
                     return sphere + cube;
                   }});

  cases.push_back({"csg_difference_hole", []() {
                     auto cube = Manifold::Cube(vec3(10.0), true);
                     auto cylinder = Manifold::Cylinder(20.0, 2.0, 2.0, 32, true);
                     return cube - cylinder;
                   }});

  cases.push_back({"csg_intersection_partial", []() {
                     auto cube = Manifold::Cube(vec3(8.0), true);
                     auto sphere = Manifold::Sphere(6.0, 32).Translate(vec3(3.0, 0.0, 0.0));
                     return cube ^ sphere;
                   }});

  // Edge case: externally tangent spheres
  cases.push_back({"edge_tangent_spheres", []() {
                     auto s1 = Manifold::Sphere(5.0, 32);
                     auto s2 = Manifold::Sphere(5.0, 32).Translate(vec3(10.0, 0.0, 0.0));
                     return s1 + s2;
                   }});

  // Edge case: coincident faces (shared face between cubes)
  cases.push_back({"edge_coincident_faces", []() {
                     auto c1 = Manifold::Cube(vec3(10.0), true);
                     auto c2 = Manifold::Cube(vec3(10.0), true).Translate(vec3(10.0, 0.0, 0.0));
                     return c1 + c2;
                   }});

  // Edge case: fully contained operand
  cases.push_back({"edge_contained", []() {
                     auto outer = Manifold::Cube(vec3(20.0), true);
                     auto inner = Manifold::Cube(vec3(5.0), true);
                     return outer - inner;
                   }});

  // ============== TIER 5: Chained CSG Operations ==============
  cases.push_back({"csg_chain_2ops", []() {
                     auto cube = Manifold::Cube(vec3(10.0), true);
                     auto sphere = Manifold::Sphere(5.0, 32);
                     auto cylinder = Manifold::Cylinder(20.0, 2.0, 2.0, 32, true);
                     return (cube + sphere) - cylinder;
                   }});

  cases.push_back({"csg_nested_3ops", []() {
                     auto cube = Manifold::Cube(vec3(10.0), true);
                     auto sphere = Manifold::Sphere(6.0, 32);
                     auto cylinder = Manifold::Cylinder(8.0, 3.0, 3.0, 32, true);
                     auto box = Manifold::Cube(vec3(15.0), true);
                     return ((cube - sphere) + cylinder) ^ box;
                   }});

  cases.push_back({"csg_symmetric_diff", []() {
                     auto a = Manifold::Cube(vec3(10.0), true);
                     auto b = Manifold::Sphere(6.0, 32).Translate(vec3(3.0, 0.0, 0.0));
                     // Symmetric difference: (A ∪ B) − (A ∩ B)
                     return (a + b) - (a ^ b);
                   }});

  // ============== TIER 6: Complex CSG (Fractals) ==============
  cases.push_back({"complex_menger_d2", []() {
                     return MengerSponge(2);
                   }});

  return cases;
}

int main(int argc, char** argv) {
  std::string output_dir = "reference_stl";

  if (argc > 1) {
    output_dir = argv[1];
  }

  // Create output directory
  fs::create_directories(output_dir);

  std::cout << "=== Manifold C++ Reference STL Generator ===" << std::endl;
  std::cout << "Output directory: " << output_dir << std::endl;
  std::cout << std::endl;

  auto cases = GetTestCases();
  std::cout << "Generating " << cases.size() << " test cases..." << std::endl;

  int success = 0;
  int failed = 0;

  for (const auto& tc : cases) {
    std::cout << "[" << tc.id << "]" << std::endl;
    try {
      Manifold m = tc.generate();
      if (m.Status() == Manifold::Error::NoError && m.NumTri() > 0) {
        ExportSTL(m, output_dir, tc.id);
        success++;
      } else {
        std::cerr << "  ERROR: Invalid manifold or empty mesh" << std::endl;
        failed++;
      }
    } catch (const std::exception& e) {
      std::cerr << "  EXCEPTION: " << e.what() << std::endl;
      failed++;
    }
  }

  std::cout << std::endl;
  std::cout << "=== Summary ===" << std::endl;
  std::cout << "Success: " << success << "/" << cases.size() << std::endl;
  std::cout << "Failed:  " << failed << "/" << cases.size() << std::endl;

  return failed > 0 ? 1 : 0;
}
