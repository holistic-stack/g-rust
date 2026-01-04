# Manifold Test Specification

## Overview

This document provides a complete inventory of all tests that must be ported from C++ to Rust, plus additional tests required for comprehensive coverage.

> **CRITICAL**: All 197 C++ tests must be ported exactly. This document catalogs every test.

---

## Test Summary

| Test File | C++ Test Count | Status |
|-----------|----------------|--------|
| `boolean_test.cpp` | 39 | Must port |
| `boolean_complex_test.cpp` | 19 | Must port |
| `manifold_test.cpp` | 49 | Must port |
| `hull_test.cpp` | 9 | Must port |
| `polygon_test.cpp` | ~30 (file-based) | Must port |
| `smooth_test.cpp` | 13 | Must port |
| `sdf_test.cpp` | 9 | Must port |
| `cross_section_test.cpp` | 14 | Must port |
| `properties_test.cpp` | 22 | Must port |
| `samples_test.cpp` | 12 | Must port |
| `manifoldc_test.cpp` | 9 | Must port (C bindings) |
| `stl_intersection_test.cpp` | 2 | Must port |
| **Total** | **197+** | |

> **Note**: `polygon_test.cpp` uses file-based tests loaded from `test/polygons/` directory. These must also be ported.

---

## 1. C++ Test File Inventory

### 1.1 boolean_test.cpp (39 tests)

| Test Name | Line | Description | Complexity |
|-----------|------|-------------|------------|
| `Boolean::Tetra` | - | Simplest boolean: tetrahedron difference | Basic |
| `Boolean::MeshGLRoundTrip` | - | MeshGL serialization round-trip | Basic |
| `Boolean::Normals` | - | Normal preservation through boolean ops | Medium |
| `Boolean::EmptyOriginal` | - | Boolean with disjoint shapes | Basic |
| `Boolean::Mirrored` | - | Boolean with mirrored (negative scale) shapes | Medium |
| `Boolean::Cubes` | - | Union of multiple cubes | Medium |
| `Boolean::Simplify` | - | Mesh simplification after boolean | Medium |
| `Boolean::DISABLED_SimplifyCracks` | - | Crack detection in simplification | Complex |
| `Boolean::NoRetainedVerts` | - | Boolean without retained vertices | Medium |
| `Boolean::PropertiesNoIntersection` | - | Properties with no geometric intersection | Basic |
| `Boolean::MixedProperties` | - | Boolean with mixed property counts | Medium |
| `Boolean::MixedNumProp` | - | Different numProp between operands | Medium |
| `Boolean::PropsMismatch` | - | Property count mismatch handling | Medium |
| `Boolean::UnionDifference` | - | Union followed by difference | Medium |
| `Boolean::TreeTransforms` | - | CSG tree with transforms | Medium |
| `Boolean::CreatePropertiesSlow` | - | Performance: large property creation | Stress |
| `Boolean::SelfSubtract` | - | Self-subtraction yields empty | Basic |
| `Boolean::Perturb` | - | Symbolic perturbation basic | Critical |
| `Boolean::Perturb1` | - | Symbolic perturbation complex | Critical |
| `Boolean::Perturb2` | - | Symbolic perturbation stress | Critical |
| `Boolean::Coplanar` | - | Coplanar face handling | Critical |
| `Boolean::MultiCoplanar` | - | Multiple coplanar faces | Critical |
| `Boolean::AlmostCoplanar` | - | Near-coplanar perturbation | Critical |
| `Boolean::FaceUnion` | - | Union at shared face | Medium |
| `Boolean::EdgeUnion` | - | Union at shared edge | Medium |
| `Boolean::EdgeUnion2` | - | Union at shared edge variant | Medium |
| `Boolean::CornerUnion` | - | Union at shared corner | Medium |
| `Boolean::DebugTangentFaces` | - | Tangent face debugging | Complex |
| `Boolean::Split` | - | Split operation | Medium |
| `Boolean::SplitByPlane` | - | Split by plane | Medium |
| `Boolean::SplitByPlane60` | - | Split by angled plane | Medium |
| `Boolean::Vug` | - | Vug cavity creation | Complex |
| `Boolean::Empty` | - | Empty result handling | Basic |
| `Boolean::Winding` | - | Winding number tests | Medium |
| `Boolean::NonIntersecting` | - | Disjoint shapes boolean | Basic |
| `Boolean::Precision` | - | Numerical precision test | Critical |
| `Boolean::Precision2` | - | Numerical precision variant | Critical |
| `Boolean::SimpleCubeRegression` | - | Simple cube regression test | Basic |
| `Boolean::BatchBoolean` | - | Batch boolean operations | Medium |

### 1.2 boolean_complex_test.cpp (19 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `BooleanComplex::Sphere` | Sphere boolean operations | Medium |
| `BooleanComplex::MeshRelation` | Mesh relation tracking through booleans | Complex |
| `BooleanComplex::Cylinders` | Multiple cylinder intersections | Complex |
| `BooleanComplex::Subtract` | Complex subtraction | Complex |
| `BooleanComplex::Close` | Nearly touching geometry | Critical |
| `BooleanComplex::BooleanVolumes` | Volume preservation tests | Medium |
| `BooleanComplex::Spiral` | Spiral shape operations | Complex |
| `BooleanComplex::Sweep` | Sweep operations | Complex |
| `BooleanComplex::InterpolatedNormals` | Normal interpolation | Medium |
| `BooleanComplex::Ring` | Ring/torus operations | Complex |
| `BooleanComplex::SelfIntersect` | Self-intersection handling | Critical |
| `BooleanComplex::GenericTwinBooleanTest7081` | Regression test #7081 | Critical |
| `BooleanComplex::GenericTwinBooleanTest7863` | Regression test #7863 | Critical |
| `BooleanComplex::Havocglass8Bool` | Havocglass regression | Complex |
| `BooleanComplex::CraycloudBool` | Complex cloud boolean | Complex |
| `BooleanComplex::HullMask` | Hull with mask | Medium |
| `BooleanComplex::SimpleOffset` | Simple offset operations | Medium |
| `BooleanComplex::DISABLED_OffsetTriangulationFailure` | Offset triangulation failure (disabled) | Complex |
| `BooleanComplex::DISABLED_OffsetSelfIntersect` | Offset self-intersection (disabled) | Complex |

### 1.3 manifold_test.cpp (49 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `Manifold::GetMeshGL` | MeshGL export | Basic |
| `Manifold::MeshDeterminism` | Deterministic mesh output | Critical |
| `Manifold::Empty` | Empty manifold handling | Basic |
| `Manifold::ValidInput` | Valid mesh construction | Basic |
| `Manifold::ValidInputOneRunIndex` | Single run index | Basic |
| `Manifold::InvalidInput1` | NaN vertex rejection | Basic |
| `Manifold::InvalidInput2` | Non-manifold rejection | Basic |
| `Manifold::InvalidInput3` | Negative vertex index | Basic |
| `Manifold::InvalidInput4` | Out of bounds vertex | Basic |
| `Manifold::InvalidInput5` | Merge index out of bounds | Basic |
| `Manifold::InvalidInput6` | Triangle index out of bounds | Basic |
| `Manifold::InvalidInput7` | Run index wrong length | Basic |
| `Manifold::OppositeFace` | Opposite face handling | Medium |
| `Manifold::Decompose` | Mesh decomposition | Medium |
| `Manifold::DecomposeProps` | Decomposition with properties | Medium |
| `Manifold::Sphere` | Sphere construction | Basic |
| `Manifold::Cylinder` | Cylinder construction | Basic |
| `Manifold::Extrude` | Extrusion | Medium |
| `Manifold::ExtrudeCone` | Cone extrusion | Medium |
| `Manifold::Revolve` | Revolve operation | Medium |
| `Manifold::Revolve2` | Revolve variant 2 | Medium |
| `Manifold::Revolve3` | Revolve variant 3 | Medium |
| `Manifold::RevolveClip` | Revolve with clipping | Medium |
| `Manifold::PartialRevolveOnYAxis` | Partial Y-axis revolve | Medium |
| `Manifold::PartialRevolveOffset` | Partial revolve with offset | Medium |
| `Manifold::Warp` | Warp function application | Medium |
| `Manifold::Warp2` | Warp variant 2 | Medium |
| `Manifold::WarpBatch` | Batch warp | Medium |
| `Manifold::Project` | Projection operation | Medium |
| `Manifold::Transform` | Transform application | Basic |
| `Manifold::Slice` | Slice operation | Medium |
| `Manifold::SliceEmptyObject` | Slice empty object | Basic |
| `Manifold::Simplify` | Mesh simplification | Medium |
| `Manifold::MeshID` | Mesh ID tracking | Medium |
| `Manifold::MeshRelation` | Mesh relation | Medium |
| `Manifold::MeshRelationTransform` | Mesh relation with transform | Medium |
| `Manifold::MeshRelationRefine` | Mesh relation with refinement | Medium |
| `Manifold::MeshRelationRefinePrecision` | Refinement precision | Medium |
| `Manifold::MeshGLRoundTrip` | MeshGL round-trip | Basic |
| `Manifold::Merge` | Mesh merging | Medium |
| `Manifold::MergeEmpty` | Merge empty | Basic |
| `Manifold::PinchedVert` | Pinched vertex handling | Medium |
| `Manifold::FaceIDRoundTrip` | Face ID round-trip | Medium |
| `Manifold::MirrorUnion` | Mirror union | Medium |
| `Manifold::MirrorUnion2` | Mirror union variant | Medium |
| `Manifold::Invalid` | General invalid input | Basic |
| `Manifold::MultiCompose` | Multiple composition | Medium |
| `Manifold::MergeDegenerates` | Merge degenerate triangles | Medium |
| `Manifold::MergeRefine` | Merge with refinement | Medium |
| `Manifold::MergeDegenerates` | Merge degenerates | Medium |
| `Manifold::MeshRelation` | Mesh relation tracking | Medium |
| `Manifold::MeshRelationRefine` | Mesh relation with refine | Medium |
| `Manifold::MeshRelationRefinePrecision` | Precision refine | Critical |
| `Manifold::MeshRelationTransform` | Relation with transform | Medium |
| `Manifold::MeshID` | Mesh ID handling | Basic |
| `Manifold::FaceIDRoundTrip` | Face ID round-trip | Basic |
| `Manifold::PinchedVert` | Pinched vertex handling | Medium |

### 1.4 hull_test.cpp (9 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `Hull::Tictac` | Tic-tac shape hull | Medium |
| `Hull::Hollow` | Hollow shape hull | Medium |
| `Hull::Cube` | Cube points hull | Basic |
| `Hull::Empty` | Empty/coplanar rejection | Basic |
| `Hull::MengerSponge` | Menger sponge hull | Medium |
| `Hull::Sphere` | Sphere hull | Basic |
| `Hull::FailingTest1` | Regression case 1 | Critical |
| `Hull::FailingTest2` | Regression case 2 | Critical |
| `Hull::DisabledFaceTest` | Disabled face handling | Medium |

### 1.5 polygon_test.cpp (6 tests)

Uses polygon corpus files for data-driven testing.

| Test Name | Description | Source File |
|-----------|-------------|-------------|
| Polygon corpus | Various polygon shapes | `polygons/polygon_corpus.txt` |
| Sponge | Sponge-like polygons | `polygons/sponge.txt` |
| Zebra | Zebra pattern polygons | `polygons/zebra.txt` |
| Zebra3 | Complex zebra patterns | `polygons/zebra3.txt` |

**Test Properties:**
- Turn 180° invariance
- Duplication correctness
- Epsilon tolerance

### 1.6 smooth_test.cpp (14 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `Smooth::Tetrahedron` | Smooth tetrahedron | Basic |
| `Smooth::RefineQuads` | Quad refinement | Medium |
| `Smooth::TruncatedCone` | Truncated cone smoothing | Medium |
| `Smooth::ToLength` | Smooth to length | Medium |
| `Smooth::Sphere` | Sphere smoothing | Medium |
| `Smooth::Precision` | Smoothing precision | Critical |
| `Smooth::Normals` | Normal-based smoothing | Medium |
| `Smooth::Manual` | Manual tangent control | Complex |
| `Smooth::Mirrored` | Mirrored smooth shapes | Medium |
| `Smooth::Csaszar` | Csaszar polyhedron | Complex |
| `Smooth::Torus` | Torus smoothing | Complex |
| `Smooth::SDF` | SDF-based smoothing | Complex |
| `Smooth::SineSurface` | Sine surface smoothing | Complex |

### 1.7 sdf_test.cpp (9 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `SDF::Blobs` | Blob SDF | Medium |
| `SDF::Bounds` | SDF bounds | Basic |
| `SDF::Bounds2` | SDF bounds variant 2 | Basic |
| `SDF::Bounds3` | SDF bounds variant 3 | Basic |
| `SDF::CubeVoid` | Cube with void SDF | Medium |
| `SDF::Resize` | SDF resize | Medium |
| `SDF::SineSurface` | Sine surface SDF | Complex |
| `SDF::SphereShell` | Sphere shell SDF | Medium |
| `SDF::Void` | Void SDF | Basic |

### 1.8 cross_section_test.cpp (14 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `CrossSection::Square` | Square cross-section | Basic |
| `CrossSection::Rect` | Rectangle cross-section | Basic |
| `CrossSection::Empty` | Empty cross-section | Basic |
| `CrossSection::Transform` | Cross-section transform | Basic |
| `CrossSection::Warp` | Cross-section warp | Medium |
| `CrossSection::Hull` | Cross-section convex hull | Medium |
| `CrossSection::HullError` | Hull error handling | Basic |
| `CrossSection::Decompose` | Cross-section decompose | Medium |
| `CrossSection::BatchBoolean` | Batch 2D boolean ops | Medium |
| `CrossSection::BevelOffset` | Bevel offset | Medium |
| `CrossSection::RoundOffset` | Round offset | Medium |
| `CrossSection::FillRule` | Fill rule handling | Medium |
| `CrossSection::MirrorCheckAxis` | Mirror axis check | Basic |
| `CrossSection::MirrorUnion` | Mirror union | Medium |

### 1.9 properties_test.cpp (21 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `Properties::Measurements` | Basic measurements | Basic |
| `Properties::Epsilon` | Epsilon handling | Basic |
| `Properties::Epsilon2` | Epsilon variant 2 | Basic |
| `Properties::Tolerance` | Tolerance testing | Medium |
| `Properties::ToleranceSphere` | Sphere tolerance | Medium |
| `Properties::CalculateCurvature` | Curvature calculation | Medium |
| `Properties::Coplanar` | Coplanar handling | Medium |
| `Properties::MinGapCubeCube` | Min gap cube-cube | Medium |
| `Properties::MinGapCubeCube2` | Min gap variant 2 | Medium |
| `Properties::MinGapCubeSphereOverlapping` | Overlapping min gap | Medium |
| `Properties::MinGapSphereSphere` | Sphere-sphere min gap | Medium |
| `Properties::MinGapSphereSphereOutOfBounds` | Out of bounds min gap | Medium |
| `Properties::MinGapClosestPointOnEdge` | Edge closest point | Medium |
| `Properties::MinGapClosestPointOnTriangleFace` | Face closest point | Medium |
| `Properties::MingapAfterTransformations` | Post-transform min gap | Medium |
| `Properties::MingapStretchyBracelet` | Stretchy bracelet | Complex |
| `Properties::MinGapAfterTransformationsOutOfBounds` | Transform OOB | Medium |
| `Properties::TriangleDistanceClosestPointsOnVertices` | Vertex distance | Medium |
| `Properties::TriangleDistanceClosestPointOnEdge` | Edge distance | Medium |
| `Properties::TriangleDistanceClosestPointOnEdge2` | Edge distance variant | Medium |
| `Properties::TriangleDistanceClosestPointOnFace` | Face distance | Medium |
| `Properties::TriangleDistanceOverlapping` | Overlapping triangles | Medium |

### 1.10 samples_test.cpp (12 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `Samples::Bracelet` | Bracelet sample | Medium |
| `Samples::CondensedMatter16` | Condensed matter 16 | Complex |
| `Samples::CondensedMatter64` | Condensed matter 64 | Complex |
| `Samples::Frame` | Frame sample | Medium |
| `Samples::FrameReduced` | Reduced frame | Medium |
| `Samples::GyroidModule` | Gyroid module | Complex |
| `Samples::Knot13` | Knot 1,3 | Medium |
| `Samples::Knot42` | Knot 4,2 | Medium |
| `Samples::Scallop` | Scallop sample | Medium |
| `Samples::Sponge1` | Sponge depth 1 | Medium |
| `Samples::Sponge4` | Sponge depth 4 | Stress |
| `Samples::TetPuzzle` | Tetrahedral puzzle | Complex |

### 1.11 manifoldc_test.cpp (8 tests - C Bindings)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `CBIND::sphere` | C binding sphere | Basic |
| `CBIND::extrude` | C binding extrude | Basic |
| `CBIND::level_set` | C binding level set | Medium |
| `CBIND::level_set_64` | C binding level set 64 | Medium |
| `CBIND::compose_decompose` | C binding compose/decompose | Medium |
| `CBIND::polygons` | C binding polygons | Basic |
| `CBIND::properties` | C binding properties | Medium |
| `CBIND::triangulation` | C binding triangulation | Basic |
| `CBIND::warp_translation` | C binding warp | Medium |

### 1.12 sphere_segments_test.cpp (2 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| Sphere segments basic | Basic segment test | Basic |
| Sphere segments high | High segment count | Medium |

### 1.13 stl_intersection_test.cpp (2 tests)

| Test Name | Description | Complexity |
|-----------|-------------|------------|
| `STLDebug::CsgIntersectionCentered` | Centered CSG intersection | Medium |
| `STLDebug::CsgIntersectionPartial` | Partial CSG intersection | Medium |

### 1.14 Debug Test Files (For Development)

These are not required for the Rust port but useful for debugging:

| File | Description |
|------|-------------|
| `debug_boolean_test.cpp` | Boolean debugging |
| `debug_csg_test.cpp` | CSG debugging |
| `debug_difference_hole.cpp` | Difference hole debugging |
| `debug_subdivision.cpp` | Subdivision debugging |

---

## 2. Additional Required Tests (Beyond C++)

### 2.1 Atomic Unit Tests (Not in C++ directly)

```rust
mod atomic_math_tests {
    // Interpolate function edge cases
    #[test] fn interpolate_at_left_boundary() { }
    #[test] fn interpolate_at_right_boundary() { }
    #[test] fn interpolate_exact_midpoint() { }
    #[test] fn interpolate_nan_left() { }
    #[test] fn interpolate_nan_right() { }
    #[test] fn interpolate_inf_values() { }
    #[test] fn interpolate_zero_delta() { }
    
    // Intersect function edge cases
    #[test] fn intersect_parallel_lines() { }
    #[test] fn intersect_coincident_lines() { }
    #[test] fn intersect_at_endpoint() { }
    #[test] fn intersect_nan_input() { }
    
    // Shadows predicate
    #[test] fn shadows_p_equals_q_negative_dir() { }
    #[test] fn shadows_p_equals_q_positive_dir() { }
    #[test] fn shadows_p_equals_q_zero_dir() { }
    #[test] fn shadows_p_less_q() { }
    #[test] fn shadows_p_greater_q() { }
    #[test] fn shadows_epsilon_difference() { }
    
    // Trigonometric exactness
    #[test] fn sind_exact_0() { }
    #[test] fn sind_exact_30() { }
    #[test] fn sind_exact_45() { }
    #[test] fn sind_exact_60() { }
    #[test] fn sind_exact_90() { }
    #[test] fn sind_exact_180() { }
    #[test] fn sind_exact_270() { }
    #[test] fn sind_exact_360() { }
    #[test] fn sind_exact_negative() { }
    #[test] fn sind_exact_large() { }
    #[test] fn cosd_exact_all_quadrants() { }
}
```

### 2.2 Box Tests

```rust
mod box_tests {
    #[test] fn box_default_contains_all() { }
    #[test] fn box_from_two_points() { }
    #[test] fn box_size() { }
    #[test] fn box_center() { }
    #[test] fn box_scale() { }
    #[test] fn box_contains_point_inside() { }
    #[test] fn box_contains_point_on_face() { }
    #[test] fn box_contains_point_on_edge() { }
    #[test] fn box_contains_point_on_corner() { }
    #[test] fn box_contains_point_outside() { }
    #[test] fn box_overlap_partial() { }
    #[test] fn box_overlap_contained() { }
    #[test] fn box_overlap_touching_face() { }
    #[test] fn box_overlap_touching_edge() { }
    #[test] fn box_overlap_touching_corner() { }
    #[test] fn box_no_overlap_separated() { }
    #[test] fn box_union() { }
    #[test] fn box_transform() { }
}
```

### 2.3 Halfedge Data Structure Tests

```rust
mod halfedge_tests {
    #[test] fn halfedge_is_forward_positive() { }
    #[test] fn halfedge_is_forward_negative() { }
    #[test] fn halfedge_is_forward_equal() { }
    #[test] fn halfedge_ordering() { }
    #[test] fn next_halfedge_within_triangle() { }
    #[test] fn next_halfedge_wrap_around() { }
    #[test] fn create_tmp_edges() { }
    #[test] fn tmp_edge_ordering() { }
}
```

### 2.4 Collider BVH Tests

```rust
mod collider_tests {
    #[test] fn collider_empty() { }
    #[test] fn collider_single_box() { }
    #[test] fn collider_two_boxes_overlap() { }
    #[test] fn collider_two_boxes_separate() { }
    #[test] fn collider_many_boxes_random() { }
    #[test] fn collider_morton_code_order() { }
    #[test] fn collider_tree_structure() { }
    #[test] fn collider_self_collision() { }
    #[test] fn collider_query_all() { }
    #[test] fn collider_query_partial() { }
}
```

### 2.5 Parallel Operations Tests

```rust
mod parallel_tests {
    #[test] fn auto_policy_below_threshold() { }
    #[test] fn auto_policy_at_threshold() { }
    #[test] fn auto_policy_above_threshold() { }
    
    #[test] fn for_each_seq_small() { }
    #[test] fn for_each_par_large() { }
    #[test] fn for_each_n_range() { }
    
    #[test] fn transform_seq() { }
    #[test] fn transform_par() { }
    
    #[test] fn reduce_sum() { }
    #[test] fn reduce_min() { }
    #[test] fn reduce_max() { }
    #[test] fn reduce_custom() { }
    
    #[test] fn exclusive_scan_sum() { }
    #[test] fn exclusive_scan_custom() { }
    #[test] fn inclusive_scan_sum() { }
    
    #[test] fn copy_if_all() { }
    #[test] fn copy_if_none() { }
    #[test] fn copy_if_partial() { }
    
    #[test] fn stable_sort_already_sorted() { }
    #[test] fn stable_sort_reversed() { }
    #[test] fn stable_sort_random() { }
    #[test] fn stable_sort_preserves_order() { }
}
```

### 2.6 DisjointSets Tests

```rust
mod disjoint_sets_tests {
    #[test] fn disjoint_sets_initial() { }
    #[test] fn disjoint_sets_find_self() { }
    #[test] fn disjoint_sets_unite_two() { }
    #[test] fn disjoint_sets_unite_chain() { }
    #[test] fn disjoint_sets_same_component() { }
    #[test] fn disjoint_sets_different_component() { }
    #[test] fn disjoint_sets_path_compression() { }
    #[test] fn disjoint_sets_rank_heuristic() { }
    #[test] fn disjoint_sets_connected_components() { }
    #[test] fn disjoint_sets_concurrent_unite() { }
    #[test] fn disjoint_sets_concurrent_find() { }
}
```

### 2.7 Numerical Edge Cases

```rust
mod numerical_edge_cases {
    // Epsilon handling
    #[test] fn epsilon_equal_vertices() { }
    #[test] fn epsilon_nearly_coplanar() { }
    #[test] fn epsilon_nearly_collinear() { }
    
    // Degenerate triangles
    #[test] fn degenerate_zero_area() { }
    #[test] fn degenerate_collinear_points() { }
    #[test] fn degenerate_duplicate_vertices() { }
    
    // Scale extremes
    #[test] fn scale_very_small() { }
    #[test] fn scale_very_large() { }
    #[test] fn scale_mixed() { }
    
    // Coordinate extremes
    #[test] fn coords_near_origin() { }
    #[test] fn coords_far_from_origin() { }
    #[test] fn coords_mixed_magnitude() { }
}
```

### 2.8 Stress Tests (Beyond C++)

```rust
mod stress_tests {
    // Many operations
    #[test] fn stress_100_unions() { }
    #[test] fn stress_100_differences() { }
    #[test] fn stress_100_intersections() { }
    #[test] fn stress_mixed_1000_ops() { }
    
    // Large meshes
    #[test] fn stress_sphere_10000_tris() { }
    #[test] fn stress_sphere_100000_tris() { }
    
    // Deep nesting
    #[test] fn stress_csg_depth_20() { }
    #[test] fn stress_csg_depth_50() { }
    #[test] fn stress_csg_depth_100() { }
    
    // Many small objects
    #[test] fn stress_1000_small_cubes() { }
    #[test] fn stress_10000_small_cubes() { }
    
    // Repeated operations
    #[test] fn stress_repeated_union_same() { }
    #[test] fn stress_repeated_difference_same() { }
}
```

### 2.9 Regression Tests

```rust
mod regression_tests {
    // From GitHub issues (port any failing cases)
    #[test] fn regression_issue_xxx() { }
    
    // From fuzzing discoveries
    #[test] fn fuzz_discovered_case_001() { }
}
```

---

## 3. Golden Test Data Generation

### 3.1 Golden Test Categories

| Category | Count | Description |
|----------|-------|-------------|
| Primitive shapes | 50 | Cubes, spheres, cylinders, cones |
| Basic booleans | 100 | Simple two-operand operations |
| Complex booleans | 200 | Multi-operand, nested CSG |
| Edge cases | 150 | Coplanar, touching, near-miss |
| Transformations | 100 | Scale, rotate, translate |
| Mesh operations | 100 | Refine, simplify, decompose |
| Hull operations | 50 | Various point cloud hulls |
| Smoothing | 50 | Various smooth shapes |
| Properties | 50 | Property interpolation |
| **Total** | **850+** | |

### 3.2 Golden Data Format

```json
{
  "test_name": "cube_sphere_difference",
  "operation": "difference",
  "operands": ["cube_1x1x1_centered", "sphere_r0.7"],
  "expected": {
    "num_vert": 42,
    "num_tri": 80,
    "num_edge": 120,
    "genus": 0,
    "volume": 0.5648,
    "surface_area": 5.234,
    "bounding_box": {
      "min": [-0.5, -0.5, -0.5],
      "max": [0.5, 0.5, 0.5]
    },
    "is_manifold": true,
    "matches_tri_normals": true
  },
  "mesh_hash": "sha256:abcd1234...",
  "tolerance": 1e-10
}
```

### 3.3 Golden Test Generation Script

```rust
// golden_gen/main.rs
fn main() {
    let tests = vec![
        ("tetra_diff", tetra_difference()),
        ("cube_union", cube_union()),
        // ... all test cases
    ];
    
    for (name, manifold) in tests {
        let golden = GoldenData {
            num_vert: manifold.num_vert(),
            num_tri: manifold.num_tri(),
            genus: manifold.genus(),
            volume: manifold.volume(),
            surface_area: manifold.surface_area(),
            // ...
        };
        save_golden(&format!("{}.json", name), &golden);
    }
}
```

---

## 4. Property-Based Test Specifications

### 4.1 Boolean Properties

```rust
proptest! {
    // Volume bounds
    #[test]
    fn prop_union_volume_bounds(a: TestManifold, b: TestManifold) {
        let union = a.clone() + b.clone();
        let intersect = a.clone() ^ b.clone();
        
        // Union volume ≤ sum of volumes
        assert!(union.volume() <= a.volume() + b.volume() + EPSILON);
        
        // Intersection volume ≤ min volume
        assert!(intersect.volume() <= a.volume().min(b.volume()) + EPSILON);
        
        // Inclusion-exclusion (approximately)
        let expected = a.volume() + b.volume() - intersect.volume();
        assert!((union.volume() - expected).abs() < EPSILON);
    }
    
    // Manifold preservation
    #[test]
    fn prop_boolean_preserves_manifold(a: TestManifold, b: TestManifold, op: BoolOp) {
        let result = apply_op(a, b, op);
        assert!(result.is_empty() || result.is_manifold());
    }
    
    // Self operations
    #[test]
    fn prop_self_subtract_empty(m: TestManifold) {
        let result = m.clone() - m;
        assert!(result.is_empty() || result.volume() < EPSILON);
    }
    
    #[test]
    fn prop_self_union_same(m: TestManifold) {
        let result = m.clone() + m.clone();
        assert!((result.volume() - m.volume()).abs() < EPSILON);
    }
    
    #[test]
    fn prop_self_intersect_same(m: TestManifold) {
        let result = m.clone() ^ m.clone();
        assert!((result.volume() - m.volume()).abs() < EPSILON);
    }
    
    // Commutativity
    #[test]
    fn prop_union_commutative(a: TestManifold, b: TestManifold) {
        let r1 = a.clone() + b.clone();
        let r2 = b + a;
        assert!((r1.volume() - r2.volume()).abs() < EPSILON);
    }
    
    #[test]
    fn prop_intersect_commutative(a: TestManifold, b: TestManifold) {
        let r1 = a.clone() ^ b.clone();
        let r2 = b ^ a;
        assert!((r1.volume() - r2.volume()).abs() < EPSILON);
    }
    
    // Associativity
    #[test]
    fn prop_union_associative(a: TestManifold, b: TestManifold, c: TestManifold) {
        let r1 = (a.clone() + b.clone()) + c.clone();
        let r2 = a + (b + c);
        assert!((r1.volume() - r2.volume()).abs() < EPSILON);
    }
}
```

### 4.2 Transform Properties

```rust
proptest! {
    // Scale volume
    #[test]
    fn prop_scale_volume(m: TestManifold, s: f64) {
        let scaled = m.scale(DVec3::splat(s.abs() + 0.001));
        let expected = m.volume() * (s.abs() + 0.001).powi(3);
        assert!((scaled.volume() - expected).abs() < EPSILON * expected);
    }
    
    // Translate preserves volume
    #[test]
    fn prop_translate_preserves_volume(m: TestManifold, t: DVec3) {
        let translated = m.translate(t);
        assert!((translated.volume() - m.volume()).abs() < EPSILON);
    }
    
    // Rotate preserves volume
    #[test]
    fn prop_rotate_preserves_volume(m: TestManifold, angles: (f64, f64, f64)) {
        let rotated = m.rotate(angles.0, angles.1, angles.2);
        assert!((rotated.volume() - m.volume()).abs() < EPSILON);
    }
}
```

---

## 5. Fuzz Test Specifications

### 5.1 Boolean Fuzzing

```rust
// fuzz/fuzz_targets/boolean_fuzz.rs
#![no_main]
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary)]
struct FuzzInput {
    seed1: u64,
    seed2: u64,
    op: u8,
    transform: [f64; 12],
}

fuzz_target!(|input: FuzzInput| {
    let m1 = make_manifold(input.seed1);
    let m2 = make_manifold(input.seed2).transform(input.transform);
    
    let result = match input.op % 3 {
        0 => m1 + m2,
        1 => m1 - m2,
        2 => m1 ^ m2,
        _ => unreachable!(),
    };
    
    // Must not panic, must be valid
    assert!(result.status() == Error::NoError || result.is_empty());
});
```

### 5.2 Polygon Fuzzing

```rust
// fuzz/fuzz_targets/polygon_fuzz.rs
#![no_main]
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary)]
struct FuzzPolygon {
    points: Vec<(f64, f64)>,
}

fuzz_target!(|input: FuzzPolygon| {
    if input.points.len() < 3 {
        return;
    }
    
    let polygon: SimplePolygon = input.points.iter()
        .map(|(x, y)| DVec2::new(*x, *y))
        .collect();
    
    // Must not panic
    let _ = triangulate(&[polygon], -1.0);
});
```

---

## 6. Continuous Integration Test Plan

### 6.1 PR Tests (Fast, <5 min)

- All unit tests
- Basic golden tests (100 samples)
- Basic property tests (1000 iterations)

### 6.2 Nightly Tests (Comprehensive, ~30 min)

- All unit tests
- All golden tests
- Full property tests (100,000 iterations)
- Stress tests

### 6.3 Weekly Tests (Exhaustive, ~4 hours)

- All unit tests
- All golden tests
- Full property tests (1,000,000 iterations)
- All stress tests
- Fuzz testing (1 hour)

### 6.4 Release Tests

- All of the above
- Performance benchmarks vs C++
- Memory safety validation (Miri)
- Documentation tests

---

## Appendix: Test Data Files

### Required Polygon Corpus Files

Copy from `submodules/manifold/test/polygons/`:
- `polygon_corpus.txt`
- `sponge.txt`
- `zebra.txt`
- `zebra3.txt`

### Required Model Files

Copy from `submodules/manifold/test/models/`:
- Any `.stl` or `.glb` files used in tests

---

## 7. Additional Comprehensive Tests (Simplest to Complex)

This section defines additional tests organized from simplest to most complex scenarios to ensure complete correctness and coverage.

### 7.1 Level 0: Primitives and Construction

```rust
mod level0_primitives {
    //=== Empty/Trivial Cases ===
    #[test] fn empty_manifold_default() { }
    #[test] fn empty_manifold_operations() { }
    
    //=== Tetrahedron (Simplest Manifold) ===
    #[test] fn tetrahedron_vertices_4() { }
    #[test] fn tetrahedron_faces_4() { }
    #[test] fn tetrahedron_edges_6() { }
    #[test] fn tetrahedron_volume_known() { }
    #[test] fn tetrahedron_surface_area() { }
    #[test] fn tetrahedron_is_manifold() { }
    #[test] fn tetrahedron_genus_0() { }
    #[test] fn tetrahedron_bounding_box() { }
    
    //=== Cube ===
    #[test] fn cube_unit_vertices_8() { }
    #[test] fn cube_unit_faces_12() { }  // 6 faces × 2 triangles
    #[test] fn cube_unit_volume_1() { }
    #[test] fn cube_unit_surface_area_6() { }
    #[test] fn cube_scaled_volume() { }
    #[test] fn cube_centered_origin() { }
    #[test] fn cube_non_centered() { }
    
    //=== Sphere ===
    #[test] fn sphere_low_res_n4() { }
    #[test] fn sphere_mid_res_n16() { }
    #[test] fn sphere_high_res_n64() { }
    #[test] fn sphere_volume_approaches_4pi_r3() { }
    #[test] fn sphere_surface_area_approaches_4pi_r2() { }
    
    //=== Cylinder ===
    #[test] fn cylinder_n8() { }
    #[test] fn cylinder_n32() { }
    #[test] fn cylinder_volume() { }
    #[test] fn cylinder_with_caps() { }
    #[test] fn cylinder_without_caps() { }
    
    //=== Extrude ===
    #[test] fn extrude_square_to_cube() { }
    #[test] fn extrude_circle_to_cylinder() { }
    #[test] fn extrude_triangle() { }
    #[test] fn extrude_complex_polygon() { }
    #[test] fn extrude_with_hole() { }
    
    //=== Revolve ===
    #[test] fn revolve_rect_to_cylinder() { }
    #[test] fn revolve_triangle_to_cone() { }
    #[test] fn revolve_partial_90() { }
    #[test] fn revolve_partial_180() { }
    #[test] fn revolve_partial_270() { }
}
```

### 7.2 Level 1: Simple Two-Operand Booleans

```rust
mod level1_simple_booleans {
    //=== Identical Operands ===
    #[test] fn union_cube_with_itself() { }
    #[test] fn difference_cube_from_itself() { }
    #[test] fn intersection_cube_with_itself() { }
    
    //=== Disjoint (No Overlap) ===
    #[test] fn union_disjoint_cubes() { }
    #[test] fn difference_disjoint_cubes() { }
    #[test] fn intersection_disjoint_empty() { }
    
    //=== Complete Containment ===
    #[test] fn union_small_inside_large() { }
    #[test] fn difference_small_inside_large_cavity() { }
    #[test] fn intersection_small_inside_large_small() { }
    
    //=== Partial Overlap ===
    #[test] fn union_cubes_overlap_quarter() { }
    #[test] fn difference_cubes_overlap_quarter() { }
    #[test] fn intersection_cubes_overlap_quarter() { }
    #[test] fn union_cubes_overlap_half() { }
    #[test] fn difference_cubes_overlap_half() { }
    #[test] fn intersection_cubes_overlap_half() { }
    
    //=== Different Shapes ===
    #[test] fn union_cube_sphere() { }
    #[test] fn difference_cube_minus_sphere() { }
    #[test] fn difference_sphere_minus_cube() { }
    #[test] fn intersection_cube_sphere() { }
    #[test] fn union_cylinder_sphere() { }
    #[test] fn difference_cylinder_sphere() { }
    #[test] fn intersection_cylinder_sphere() { }
}
```

### 7.3 Level 2: Boundary Conditions (Critical)

```rust
mod level2_boundary_conditions {
    //=== Touching at Single Point (Corner) ===
    #[test] fn union_cubes_touch_corner() { }
    #[test] fn difference_cubes_touch_corner() { }
    #[test] fn intersection_cubes_touch_corner() { }
    
    //=== Touching at Edge ===
    #[test] fn union_cubes_touch_edge() { }
    #[test] fn difference_cubes_touch_edge() { }
    #[test] fn intersection_cubes_touch_edge() { }
    
    //=== Touching at Face (Coplanar) ===
    #[test] fn union_cubes_touch_face() { }
    #[test] fn difference_cubes_touch_face() { }
    #[test] fn intersection_cubes_touch_face() { }
    
    //=== Shared Face (Exactly Coincident) ===
    #[test] fn union_cubes_shared_face() { }
    #[test] fn difference_cubes_shared_face() { }
    #[test] fn intersection_cubes_shared_face() { }
    
    //=== Nearly Touching (Epsilon Gap) ===
    #[test] fn union_cubes_epsilon_gap_1e10() { }
    #[test] fn union_cubes_epsilon_gap_1e12() { }
    #[test] fn union_cubes_epsilon_gap_1e14() { }
    #[test] fn difference_cubes_epsilon_gap() { }
    
    //=== Nearly Coplanar ===
    #[test] fn union_nearly_coplanar_1e10() { }
    #[test] fn union_nearly_coplanar_1e12() { }
    #[test] fn intersection_nearly_coplanar() { }
    
    //=== Edge Intersection ===
    #[test] fn edge_crosses_face_center() { }
    #[test] fn edge_crosses_face_edge() { }
    #[test] fn edge_crosses_face_corner() { }
    
    //=== Vertex on Face ===
    #[test] fn vertex_exactly_on_face() { }
    #[test] fn vertex_epsilon_above_face() { }
    #[test] fn vertex_epsilon_below_face() { }
}
```

### 7.4 Level 3: Chained Operations

```rust
mod level3_chained_operations {
    //=== Simple Chains ===
    #[test] fn chain_2_unions() { }
    #[test] fn chain_3_unions() { }
    #[test] fn chain_5_unions() { }
    #[test] fn chain_10_unions() { }
    
    #[test] fn chain_2_differences() { }
    #[test] fn chain_3_differences() { }
    #[test] fn chain_multiple_holes() { }
    
    #[test] fn chain_2_intersections() { }
    #[test] fn chain_3_intersections() { }
    
    //=== Mixed Operations ===
    #[test] fn chain_union_then_difference() { }
    #[test] fn chain_difference_then_union() { }
    #[test] fn chain_union_then_intersection() { }
    #[test] fn chain_intersection_then_union() { }
    #[test] fn chain_all_three_ops() { }
    
    //=== With Transforms ===
    #[test] fn chain_with_translation() { }
    #[test] fn chain_with_rotation() { }
    #[test] fn chain_with_scale() { }
    #[test] fn chain_with_all_transforms() { }
}
```

### 7.5 Level 4: Nested CSG Trees

```rust
mod level4_nested_csg {
    //=== Shallow Nesting (Depth 2-3) ===
    #[test] fn nested_depth_2_union() { }
    #[test] fn nested_depth_2_difference() { }
    #[test] fn nested_depth_2_mixed() { }
    #[test] fn nested_depth_3_union() { }
    #[test] fn nested_depth_3_mixed() { }
    
    //=== Medium Nesting (Depth 4-6) ===
    #[test] fn nested_depth_4_union() { }
    #[test] fn nested_depth_5_mixed() { }
    #[test] fn nested_depth_6_complex() { }
    
    //=== Deep Nesting (Depth 7-10) ===
    #[test] fn nested_depth_7() { }
    #[test] fn nested_depth_8() { }
    #[test] fn nested_depth_9() { }
    #[test] fn nested_depth_10() { }
    
    //=== Very Deep Nesting (Depth 15-20) - Stress ===
    #[test] fn nested_depth_15() { }
    #[test] fn nested_depth_20() { }
    
    //=== Wide Trees (Many Siblings) ===
    #[test] fn wide_tree_10_children() { }
    #[test] fn wide_tree_20_children() { }
    #[test] fn wide_tree_50_children() { }
    
    //=== Complex Tree Shapes ===
    #[test] fn tree_balanced_binary() { }
    #[test] fn tree_left_leaning() { }
    #[test] fn tree_right_leaning() { }
    #[test] fn tree_random_structure() { }
}
```

### 7.6 Level 5: Symbolic Perturbation Scenarios

```rust
mod level5_perturbation {
    //=== Exact Coincidence ===
    #[test] fn perturb_vertices_exact_same() { }
    #[test] fn perturb_edges_exact_colinear() { }
    #[test] fn perturb_faces_exact_coplanar() { }
    
    //=== Multiple Coincidences ===
    #[test] fn perturb_two_vertices_coincident() { }
    #[test] fn perturb_three_vertices_coincident() { }
    #[test] fn perturb_edge_on_edge() { }
    #[test] fn perturb_edge_through_vertex() { }
    
    //=== Grids (Many Coincidences) ===
    #[test] fn perturb_grid_2x2() { }
    #[test] fn perturb_grid_3x3() { }
    #[test] fn perturb_grid_4x4() { }
    #[test] fn perturb_grid_5x5() { }
    
    //=== Known Difficult Cases from C++ ===
    #[test] fn perturb_cpp_perturb() { }   // From boolean_test.cpp
    #[test] fn perturb_cpp_perturb1() { }  // From boolean_test.cpp
    #[test] fn perturb_cpp_perturb2() { }  // From boolean_test.cpp
    
    //=== Stress Perturbation ===
    #[test] fn perturb_100_coplanar_faces() { }
    #[test] fn perturb_shifted_grid_10x10() { }
}
```

### 7.7 Level 6: Complex Geometry

```rust
mod level6_complex_geometry {
    //=== Complex Shapes ===
    #[test] fn gyroid_basic() { }
    #[test] fn gyroid_intersection() { }
    #[test] fn menger_sponge_depth_1() { }
    #[test] fn menger_sponge_depth_2() { }
    #[test] fn menger_sponge_depth_3() { }
    
    //=== Thin Walls ===
    #[test] fn thin_wall_1e2() { }
    #[test] fn thin_wall_1e4() { }
    #[test] fn thin_wall_1e6() { }
    
    //=== Small Gaps ===
    #[test] fn small_gap_1e2() { }
    #[test] fn small_gap_1e4() { }
    #[test] fn small_gap_1e6() { }
    
    //=== Complex Holes ===
    #[test] fn single_hole_through() { }
    #[test] fn multiple_holes_parallel() { }
    #[test] fn multiple_holes_perpendicular() { }
    #[test] fn holes_intersecting() { }
    
    //=== Concave Shapes ===
    #[test] fn concave_l_shape() { }
    #[test] fn concave_u_shape() { }
    #[test] fn concave_complex() { }
    
    //=== High Genus ===
    #[test] fn torus_genus_1() { }
    #[test] fn double_torus_genus_2() { }
    #[test] fn triple_linked_tori() { }
}
```

### 7.8 Level 7: Stress and Performance

```rust
mod level7_stress {
    //=== Large Triangle Counts ===
    #[test] fn stress_10k_triangles() { }
    #[test] fn stress_100k_triangles() { }
    #[test] fn stress_1m_triangles() { }
    
    //=== Many Operations ===
    #[test] fn stress_100_sequential_ops() { }
    #[test] fn stress_500_sequential_ops() { }
    #[test] fn stress_1000_sequential_ops() { }
    
    //=== Many Small Objects ===
    #[test] fn stress_100_small_spheres() { }
    #[test] fn stress_500_small_spheres() { }
    #[test] fn stress_1000_small_spheres() { }
    
    //=== Memory Stress ===
    #[test] fn stress_repeated_create_destroy() { }
    #[test] fn stress_deep_clone_chain() { }
    
    //=== Parallel Stress ===
    #[test] fn stress_parallel_threshold_boundary() { }
    #[test] fn stress_many_parallel_ops() { }
}
```

### 7.9 Level 8: Determinism Verification

```rust
mod level8_determinism {
    //=== Repeated Operations ===
    #[test] fn determinism_same_result_10_runs() { }
    #[test] fn determinism_same_result_100_runs() { }
    
    //=== Order Independence ===
    #[test] fn determinism_operand_order() { }
    #[test] fn determinism_tree_structure() { }
    
    //=== Thread Count Independence ===
    #[test] fn determinism_single_thread() { }
    #[test] fn determinism_multi_thread() { }
    #[test] fn determinism_max_threads() { }
    
    //=== Platform Independence ===
    #[test] fn determinism_debug_vs_release() { }
    
    //=== Floating Point ===
    #[test] fn determinism_interpolate_identical() { }
    #[test] fn determinism_intersect_identical() { }
    #[test] fn determinism_shadows_identical() { }
}
```

### 7.10 Level 9: Error Handling and Recovery

```rust
mod level9_error_handling {
    //=== Invalid Input ===
    #[test] fn error_nan_vertex() { }
    #[test] fn error_inf_vertex() { }
    #[test] fn error_degenerate_triangle() { }
    #[test] fn error_non_manifold_input() { }
    #[test] fn error_self_intersecting_input() { }
    
    //=== Boundary Errors ===
    #[test] fn error_empty_mesh() { }
    #[test] fn error_single_triangle() { }
    #[test] fn error_open_mesh() { }
    
    //=== Recovery ===
    #[test] fn recovery_from_error() { }
    #[test] fn error_does_not_corrupt_state() { }
}
```

---

## 8. Test Execution Order

Tests should be executed in this order during development:

1. **Level 0**: All primitive tests must pass first
2. **Level 1**: Simple booleans must pass
3. **Level 2**: Boundary conditions must pass
4. **Level 3**: Chained operations must pass
5. **Level 4**: Nested CSG must pass
6. **Level 5**: Perturbation tests must pass
7. **Level 6**: Complex geometry tests
8. **Level 7**: Stress tests
9. **Level 8**: Determinism verification
10. **Level 9**: Error handling

**DO NOT proceed to next level until all tests in current level pass.**

---

*Document Version: 2.0*
*Last Updated: 2026-01-04*
*Total Estimated Tests: 700+ unit tests, 850+ golden tests, 50+ property tests, 10+ fuzz targets = 1,600+ total tests*
