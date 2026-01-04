# Detailed File Migration Mapping

This document maps every file in `submodules/manifold` to its intended location in the Rust project.

> **CRITICAL**: Only use the C++ code in `submodules/manifold/` as reference. No external algorithms permitted.

## Core Library - Source Files

| C++ File (Source/Header) | Rust Destination | Description | Lines | Priority |
| :--- | :--- | :--- | :--- | :--- |
| `src/boolean3.cpp` | `crates/manifold-boolean/src/kernel/mod.rs` | Core Boolean logic (Intersect, Interpolate, Shadows) | 553 | P0 |
| `src/boolean3.h` | `crates/manifold-boolean/src/kernel/mod.rs` | Boolean3 struct, Intersections | 66 | P0 |
| `src/boolean_result.cpp` | `crates/manifold-boolean/src/kernel/assembly.rs` | Boolean output assembly | ~600 | P0 |
| `src/impl.cpp` | `crates/manifold-boolean/src/impl/mod.rs` | Manifold::Impl implementation | ~1200 | P0 |
| `src/impl.h` | `crates/manifold-boolean/src/impl/mod.rs` | Impl struct, MeshRelationD | 405 | P0 |
| `src/collider.h` | `crates/manifold-collider/src/lib.rs` | BVH implementation | 371 | P0 |
| `src/csg_tree.cpp` | `crates/manifold-boolean/src/csg.rs` | CSG tree lazy evaluation | ~800 | P1 |
| `src/csg_tree.h` | `crates/manifold-boolean/src/csg.rs` | CsgNode, CsgLeafNode | ~150 | P1 |
| `src/constructors.cpp` | `crates/manifold/src/constructors.rs` | Mesh primitives (Cube, Sphere, etc.) | 525 | P2 |
| `src/manifold.cpp` | `crates/manifold/src/lib.rs` | Public Manifold API | ~600 | P1 |
| `src/quickhull.cpp` | `crates/manifold-boolean/src/quickhull.rs` | Convex Hull algorithm | ~600 | P2 |
| `src/quickhull.h` | `crates/manifold-boolean/src/quickhull.rs` | QuickHull structs | 289 | P2 |
| `src/polygon.cpp` | `crates/manifold-polygon/src/lib.rs` | Triangulation logic | ~800 | P1 |
| `src/tree2d.cpp` | `crates/manifold-polygon/src/tree2d.rs` | 2D spatial tree | ~300 | P1 |
| `src/tree2d.h` | `crates/manifold-polygon/src/tree2d.rs` | Tree2D struct | ~100 | P1 |
| `src/smoothing.cpp` | `crates/manifold-boolean/src/smoothing.rs` | Mesh smoothing | ~600 | P3 |
| `src/subdivision.cpp` | `crates/manifold-boolean/src/subdivision.rs` | Mesh subdivision | ~500 | P3 |
| `src/properties.cpp` | `crates/manifold-boolean/src/properties.rs` | Property interpolation | ~400 | P3 |
| `src/sort.cpp` | `crates/manifold-boolean/src/sort.rs` | Morton sorting | ~200 | P1 |
| `src/edge_op.cpp` | `crates/manifold-boolean/src/impl/edge_op.rs` | Edge operations | ~400 | P1 |
| `src/face_op.cpp` | `crates/manifold-boolean/src/impl/face_op.rs` | Face operations | ~300 | P1 |
| `src/sdf.cpp` | `crates/manifold/src/sdf.rs` | Signed Distance Field | ~300 | P3 |

## Core Library - Header-Only Files

| C++ File | Rust Destination | Description | Lines | Priority |
| :--- | :--- | :--- | :--- | :--- |
| `src/parallel.h` | `crates/manifold-parallel/src/lib.rs` | **ALL** parallel primitives | 1162 | P0 |
| `src/shared.h` | `crates/manifold-types/src/lib.rs` | Halfedge, TriRef, Barycentric, utilities | 212 | P0 |
| `src/disjoint_sets.h` | `crates/manifold-types/src/disjoint_sets.rs` | Thread-safe Union-Find | 122 | P0 |
| `src/hashtable.h` | `crates/manifold-parallel/src/hashtable.rs` | Concurrent hash table | ~200 | P1 |
| `src/vec.h` | Native Rust `Vec<T>` | Custom vector (use std) | ~100 | - |
| `src/utils.h` | `crates/manifold-math/src/utils.rs` | General utilities | ~150 | P1 |
| `src/iters.h` | `crates/manifold-parallel/src/iters.rs` | Iterator utilities | ~100 | P1 |
| `src/svd.h` | `crates/manifold-math/src/svd.rs` | SVD decomposition | ~200 | P2 |
| `src/tri_dist.h` | `crates/manifold-math/src/tri_dist.rs` | Triangle distance math | ~150 | P2 |
| `src/mesh_fixes.h` | `crates/manifold-boolean/src/mesh_fixes.rs` | Mesh repair utilities | ~100 | P2 |

## Public Headers

| C++ File | Rust Destination | Description |
| :--- | :--- | :--- |
| `include/manifold/common.h` | `crates/manifold-math/src/common.rs` + `crates/manifold-types/src/lib.rs` | Constants, Enums, Box, types |
| `include/manifold/linalg.h` | Uses `glam` crate | Math type aliases |
| `include/manifold/manifold.h` | `crates/manifold/src/lib.rs` | Public API definitions |
| `include/manifold/polygon.h` | `crates/manifold-polygon/src/lib.rs` | Polygon triangulation API |
| `include/manifold/cross_section.h` | `crates/manifold-polygon/src/cross_section.rs` | 2D cross-section API |
| `include/manifold/meshIO.h` | `crates/manifold-meshio/src/lib.rs` | Mesh I/O API |
| `include/manifold/vec_view.h` | Native Rust `&[T]` / `&mut [T]` | View types |
| `include/manifold/optional_assert.h` | Rust `debug_assert!` | Debug assertions |

## Optional Modules

| C++ File | Rust Destination | Description |
| :--- | :--- | :--- |
| `src/cross_section/cross_section.cpp` | `crates/manifold-polygon/src/cross_section.rs` | 2D cross sections (Clipper2 dependency) |
| `src/meshIO/meshIO.cpp` | `crates/manifold-meshio/src/lib.rs` | Mesh I/O (Assimp-like functionality) |

## Test Files (ALL MUST BE PORTED)

| C++ Test File | Rust Test File | Test Count |
| :--- | :--- | :--- |
| `test/boolean_test.cpp` | `crates/manifold-boolean/tests/boolean.rs` | 35+ |
| `test/boolean_complex_test.cpp` | `crates/manifold-boolean/tests/boolean_complex.rs` | 40+ |
| `test/manifold_test.cpp` | `crates/manifold/tests/manifold.rs` | 40 |
| `test/polygon_test.cpp` | `crates/manifold-polygon/tests/polygon.rs` | 6 |
| `test/hull_test.cpp` | `crates/manifold-boolean/tests/hull.rs` | 9 |
| `test/smooth_test.cpp` | `crates/manifold-boolean/tests/smooth.rs` | 14 |
| `test/sdf_test.cpp` | `crates/manifold/tests/sdf.rs` | 9 |
| `test/cross_section_test.cpp` | `crates/manifold-polygon/tests/cross_section.rs` | 14 |
| `test/properties_test.cpp` | `crates/manifold-boolean/tests/properties.rs` | 21 |
| `test/samples_test.cpp` | `crates/manifold/tests/samples.rs` | 12 |
| `test/manifoldc_test.cpp` | `crates/manifold/tests/cbind.rs` | 8 (C bindings) |
| `test/sphere_segments_test.cpp` | `crates/manifold/tests/sphere_segments.rs` | 2 |
| `test/stl_intersection_test.cpp` | `crates/manifold-boolean/tests/stl_intersection.rs` | 2 |
| `test/manifold_fuzz.cpp` | `fuzz/fuzz_targets/manifold_fuzz.rs` | Fuzz target |
| `test/polygon_fuzz.cpp` | `fuzz/fuzz_targets/polygon_fuzz.rs` | Fuzz target |

**Total: 197 C++ tests to port**

## Test Support Files

| C++ File | Rust File | Description |
| :--- | :--- | :--- |
| `test/test.h` | `crates/manifold-test-utils/src/lib.rs` | Test utilities, fixtures |
| `test/test_main.cpp` | N/A (Rust test harness) | Test runner |
| `test/samples.h` | `crates/manifold-test-utils/src/samples.rs` | Sample model generators |
| `test/polygons/*.txt` | `tests/data/polygons/*.txt` | Polygon test corpus |
| `test/models/*` | `tests/data/models/*` | Test model files |

## Debug Test Files (Optional - for development)

| C++ File | Description | Port Status |
| :--- | :--- | :--- |
| `test/debug_boolean_test.cpp` | Boolean debugging utilities | Optional |
| `test/debug_csg_test.cpp` | CSG tree debugging | Optional |
| `test/debug_difference_hole.cpp` | Hole difference debugging | Optional |
| `test/debug_subdivision.cpp` | Subdivision debugging | Optional |
| `test/export_sphere.cpp` | Sphere export utility | Optional |
| `test/export_sphere_simple.cpp` | Simple sphere export | Optional |
| `test/export_sphere_vertices.cpp` | Sphere vertices export | Optional |

## Critical Functions to Port Exactly

These functions in `boolean3.cpp` (553 lines) are the heart of the Boolean algorithm:

| Function | Lines | Template Params | Priority |
| :--- | :--- | :--- | :--- |
| `withSign()` | 38 | - | P0 |
| `Interpolate()` | 40-53 | - | P0 |
| `Intersect()` | 55-72 | - | P0 |
| `Shadows()` | 74-76 | - | P0 |
| `Shadow01<>` | 78-122 | `<expandP, forward>` | P0 |
| `Kernel11<>` | 124-180 | `<expandP>` | P0 |
| `Kernel02<>` | 182-248 | `<expandP, forward>` | P0 |
| `Kernel12<>` | 250-320 | `<expandP, forward>` | P0 |
| `Kernel12Recorder<>` | 322-380 | `<expandP, forward>` | P0 |
| `Intersect12_<>` | 382-430 | `<expandP, forward>` | P0 |
| `Intersect12<>` | 432-438 | `<forward>` | P0 |
| `Winding03_<>` | 440-520 | `<expandP, forward>` | P0 |
| `Winding03<>` | 522-528 | `<forward>` | P0 |
| `Boolean3::Boolean3` | 530-553 | - | P0 |

These MUST be ported character-by-character from the C++ reference.

---

## Source Line Count Summary

| File | Lines | Priority |
| :--- | :--- | :--- |
| `parallel.h` | 1162 | P0 |
| `impl.cpp` | ~1200 | P0 |
| `boolean3.cpp` | 553 | P0 |
| `boolean_result.cpp` | ~600 | P0 |
| `constructors.cpp` | 525 | P2 |
| `polygon.cpp` | ~800 | P1 |
| `collider.h` | 371 | P0 |
| `impl.h` | 405 | P0 |
| `quickhull.cpp` | ~600 | P2 |
| `csg_tree.cpp` | ~800 | P1 |
| **Total Core** | **~7000+** | |
