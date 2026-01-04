# Implementation Checklist

This document provides a detailed checklist for tracking implementation progress during the C++ to Rust migration.

> **CRITICAL**: Every item must be implemented from the C++ reference in `submodules/manifold/`. No external algorithms permitted.

---

## Phase 1: Foundation Types (Weeks 1-3)

### 1.1 manifold-types Crate

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `Halfedge` struct | `shared.h:28-35` | `types/halfedge.rs` | ☐ |
| `Halfedge::IsForward()` | `shared.h:37-38` | `types/halfedge.rs` | ☐ |
| `Barycentric` struct | `shared.h:45-48` | `types/barycentric.rs` | ☐ |
| `TriRef` struct | `shared.h:50-55` | `types/tri_ref.rs` | ☐ |
| `TmpEdge` struct | `shared.h:57-62` | `types/tmp_edge.rs` | ☐ |
| `Smoothness` struct | `common.h` | `types/smoothness.rs` | ☐ |
| `MeshGL` struct | `common.h` | `types/mesh_gl.rs` | ☐ |
| `MeshGL64` struct | `common.h` | `types/mesh_gl.rs` | ☐ |
| `Error` enum | `common.h` | `types/error.rs` | ☐ |
| `OpType` enum | `common.h` | `types/op_type.rs` | ☐ |
| `NextHalfedge()` | `shared.h:65-68` | `types/halfedge.rs` | ☐ |
| `GetAxisAlignedProjection()` | `shared.h:103-120` | `types/projection.rs` | ☐ |

### 1.2 manifold-math Crate

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `K_PI`, `K_TWO_PI`, etc. | `common.h:30-35` | `math/constants.rs` | ☐ |
| `K_PRECISION` constant | `common.h:37` | `math/constants.rs` | ☐ |
| `radians()` | `common.h:40` | `math/trig.rs` | ☐ |
| `degrees()` | `common.h:41` | `math/trig.rs` | ☐ |
| `sind()` | `common.h:43-55` | `math/trig.rs` | ☐ |
| `cosd()` | `common.h:57` | `math/trig.rs` | ☐ |
| `smoothstep()` | `common.h:60-63` | `math/interp.rs` | ☐ |
| `Box` struct | `common.h:70-110` | `math/bbox.rs` | ☐ |
| `Box::new()` | `common.h:72-75` | `math/bbox.rs` | ☐ |
| `Box::size()` | `common.h:77` | `math/bbox.rs` | ☐ |
| `Box::center()` | `common.h:78` | `math/bbox.rs` | ☐ |
| `Box::scale()` | `common.h:79-82` | `math/bbox.rs` | ☐ |
| `Box::contains()` | `common.h:84-86` | `math/bbox.rs` | ☐ |
| `Box::does_overlap()` | `common.h:88-92` | `math/bbox.rs` | ☐ |
| `Box::union_()` | `common.h:94-97` | `math/bbox.rs` | ☐ |
| `Box::transform()` | `common.h:99-107` | `math/bbox.rs` | ☐ |
| SVD decomposition | `svd.h` | `math/svd.rs` | ☐ |
| Triangle distance | `tri_dist.h` | `math/tri_dist.rs` | ☐ |

### 1.3 manifold-parallel Crate

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `K_SEQ_THRESHOLD` | `parallel.h:27` | `parallel/policy.rs` | ☐ |
| `ExecutionPolicy` enum | `parallel.h:29-32` | `parallel/policy.rs` | ☐ |
| `autoPolicy()` | `parallel.h:34-40` | `parallel/policy.rs` | ☐ |
| `for_each()` | `parallel.h:42-52` | `parallel/ops.rs` | ☐ |
| `for_each_n()` | `parallel.h:54-64` | `parallel/ops.rs` | ☐ |
| `transform()` | `parallel.h:66-76` | `parallel/ops.rs` | ☐ |
| `reduce()` | `parallel.h:78-95` | `parallel/ops.rs` | ☐ |
| `exclusive_scan()` | `parallel.h:97-115` | `parallel/scan.rs` | ☐ |
| `inclusive_scan()` | `parallel.h:117-130` | `parallel/scan.rs` | ☐ |
| `copy_if()` | `parallel.h:132-155` | `parallel/filter.rs` | ☐ |
| `stable_sort()` | `parallel.h:157-175` | `parallel/sort.rs` | ☐ |
| `mergeRec()` | `parallel.h:74-96` | `parallel/sort.rs` | ☐ |
| `mergeSortRec()` | `parallel.h:98-110` | `parallel/sort.rs` | ☐ |
| `countAt()` | `parallel.h:177-182` | `parallel/iters.rs` | ☐ |
| `sequence()` | `parallel.h:184-190` | `parallel/iters.rs` | ☐ |
| `Permute()` | `parallel.h:192-200` | `parallel/ops.rs` | ☐ |

### 1.4 DisjointSets

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `DisjointSets` struct | `disjoint_sets.h:20-30` | `types/disjoint_sets.rs` | ☐ |
| `parent()` | `disjoint_sets.h:32-35` | `types/disjoint_sets.rs` | ☐ |
| `rank()` | `disjoint_sets.h:37-40` | `types/disjoint_sets.rs` | ☐ |
| `find()` | `disjoint_sets.h:42-55` | `types/disjoint_sets.rs` | ☐ |
| `unite()` | `disjoint_sets.h:57-85` | `types/disjoint_sets.rs` | ☐ |
| `connectedComponents()` | `disjoint_sets.h:87-110` | `types/disjoint_sets.rs` | ☐ |

---

## Phase 2: Spatial Acceleration (Weeks 4-5)

### 2.1 manifold-collider Crate

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `K_INITIAL_LENGTH` | `collider.h:30` | `collider/constants.rs` | ☐ |
| `K_LENGTH_MULTIPLE` | `collider.h:31` | `collider/constants.rs` | ☐ |
| `K_SEQUENTIAL_THRESHOLD` | `collider.h:32` | `collider/constants.rs` | ☐ |
| `K_ROOT` | `collider.h:33` | `collider/constants.rs` | ☐ |
| `IsLeaf()` | `collider.h:35` | `collider/node.rs` | ☐ |
| `IsInternal()` | `collider.h:36` | `collider/node.rs` | ☐ |
| `Node2Internal()` | `collider.h:37` | `collider/node.rs` | ☐ |
| `Internal2Node()` | `collider.h:38` | `collider/node.rs` | ☐ |
| `Node2Leaf()` | `collider.h:39` | `collider/node.rs` | ☐ |
| `Leaf2Node()` | `collider.h:40` | `collider/node.rs` | ☐ |
| `CreateRadixTree` struct | `collider.h:70-130` | `collider/radix_tree.rs` | ☐ |
| `PrefixLength()` | `collider.h:80-87` | `collider/radix_tree.rs` | ☐ |
| `FindSplit()` | `collider.h:89-110` | `collider/radix_tree.rs` | ☐ |
| `Collider` struct | `collider.h:140-200` | `collider/mod.rs` | ☐ |
| `Collider::new()` | `collider.h:200-280` | `collider/mod.rs` | ☐ |
| `Collider::Collisions()` | `collider.h:282-371` | `collider/mod.rs` | ☐ |

### 2.2 Morton Sorting

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| Morton code calculation | `sort.cpp:30-60` | `boolean/sort.rs` | ☐ |
| Morton sorting | `sort.cpp:62-120` | `boolean/sort.rs` | ☐ |

---

## Phase 3: Mesh Implementation (Weeks 6-9)

### 3.1 Impl Struct

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `Impl` struct definition | `impl.h:40-100` | `boolean/impl/mod.rs` | ☐ |
| `MeshRelationD` struct | `impl.h:102-130` | `boolean/impl/relation.rs` | ☐ |
| `Impl::Transform()` | `impl.cpp:50-80` | `boolean/impl/transform.rs` | ☐ |
| `Impl::ApplyTransform()` | `impl.cpp:82-120` | `boolean/impl/transform.rs` | ☐ |
| `Impl::CalculateNormals()` | `impl.cpp:122-180` | `boolean/impl/normals.rs` | ☐ |
| `Impl::CreateFaces()` | `impl.cpp:182-250` | `boolean/impl/faces.rs` | ☐ |
| `Impl::Finish()` | `impl.cpp:252-300` | `boolean/impl/finish.rs` | ☐ |
| `Impl::IsEmpty()` | `impl.cpp:302-310` | `boolean/impl/mod.rs` | ☐ |
| `Impl::NumVert()` | `impl.cpp:312-315` | `boolean/impl/mod.rs` | ☐ |
| `Impl::NumEdge()` | `impl.cpp:317-320` | `boolean/impl/mod.rs` | ☐ |
| `Impl::NumTri()` | `impl.cpp:322-325` | `boolean/impl/mod.rs` | ☐ |

### 3.2 Edge Operations

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `CreateHalfedges()` | `edge_op.cpp:30-100` | `boolean/impl/edge_op.rs` | ☐ |
| `RecordRun()` | `edge_op.cpp:102-150` | `boolean/impl/edge_op.rs` | ☐ |
| Edge collapse | `edge_op.cpp:152-220` | `boolean/impl/edge_op.rs` | ☐ |

### 3.3 Face Operations

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| Face splitting | `face_op.cpp:30-100` | `boolean/impl/face_op.rs` | ☐ |
| Face merging | `face_op.cpp:102-180` | `boolean/impl/face_op.rs` | ☐ |

---

## Phase 4: Boolean Kernel (Weeks 10-16) - CRITICAL

### 4.1 Core Boolean Functions (boolean3.cpp)

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `withSign()` | `boolean3.cpp:38` | `boolean/kernel/helpers.rs` | ☐ |
| `Interpolate()` | `boolean3.cpp:40-53` | `boolean/kernel/interpolate.rs` | ☐ |
| `Intersect()` | `boolean3.cpp:55-72` | `boolean/kernel/intersect.rs` | ☐ |
| `Shadows()` | `boolean3.cpp:74-76` | `boolean/kernel/shadows.rs` | ☐ |
| `Shadow01<expandP, true>` | `boolean3.cpp:78-122` | `boolean/kernel/shadow01.rs` | ☐ |
| `Shadow01<expandP, false>` | `boolean3.cpp:78-122` | `boolean/kernel/shadow01.rs` | ☐ |
| `Kernel11<true>` | `boolean3.cpp:124-180` | `boolean/kernel/kernel11.rs` | ☐ |
| `Kernel11<false>` | `boolean3.cpp:124-180` | `boolean/kernel/kernel11.rs` | ☐ |
| `Kernel02<true, true>` | `boolean3.cpp:182-248` | `boolean/kernel/kernel02.rs` | ☐ |
| `Kernel02<true, false>` | `boolean3.cpp:182-248` | `boolean/kernel/kernel02.rs` | ☐ |
| `Kernel02<false, true>` | `boolean3.cpp:182-248` | `boolean/kernel/kernel02.rs` | ☐ |
| `Kernel02<false, false>` | `boolean3.cpp:182-248` | `boolean/kernel/kernel02.rs` | ☐ |
| `Kernel12<true, true>` | `boolean3.cpp:250-320` | `boolean/kernel/kernel12.rs` | ☐ |
| `Kernel12<true, false>` | `boolean3.cpp:250-320` | `boolean/kernel/kernel12.rs` | ☐ |
| `Kernel12<false, true>` | `boolean3.cpp:250-320` | `boolean/kernel/kernel12.rs` | ☐ |
| `Kernel12<false, false>` | `boolean3.cpp:250-320` | `boolean/kernel/kernel12.rs` | ☐ |
| `Kernel12Recorder` | `boolean3.cpp:322-380` | `boolean/kernel/recorder.rs` | ☐ |
| `Intersect12_<true, true>` | `boolean3.cpp:382-430` | `boolean/kernel/intersect12.rs` | ☐ |
| `Intersect12_<true, false>` | `boolean3.cpp:382-430` | `boolean/kernel/intersect12.rs` | ☐ |
| `Intersect12_<false, true>` | `boolean3.cpp:382-430` | `boolean/kernel/intersect12.rs` | ☐ |
| `Intersect12_<false, false>` | `boolean3.cpp:382-430` | `boolean/kernel/intersect12.rs` | ☐ |
| `Intersect12<true>` | `boolean3.cpp:432-438` | `boolean/kernel/intersect12.rs` | ☐ |
| `Intersect12<false>` | `boolean3.cpp:432-438` | `boolean/kernel/intersect12.rs` | ☐ |
| `Winding03_<true, true>` | `boolean3.cpp:440-520` | `boolean/kernel/winding.rs` | ☐ |
| `Winding03_<true, false>` | `boolean3.cpp:440-520` | `boolean/kernel/winding.rs` | ☐ |
| `Winding03_<false, true>` | `boolean3.cpp:440-520` | `boolean/kernel/winding.rs` | ☐ |
| `Winding03_<false, false>` | `boolean3.cpp:440-520` | `boolean/kernel/winding.rs` | ☐ |
| `Winding03<true>` | `boolean3.cpp:522-528` | `boolean/kernel/winding.rs` | ☐ |
| `Winding03<false>` | `boolean3.cpp:522-528` | `boolean/kernel/winding.rs` | ☐ |
| `Boolean3::Boolean3` | `boolean3.cpp:530-553` | `boolean/kernel/mod.rs` | ☐ |

### 4.2 Boolean Result Assembly (boolean_result.cpp)

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `AddNewEdgeVerts()` | `boolean_result.cpp:40-100` | `boolean/assembly.rs` | ☐ |
| `DuplicateVerts()` | `boolean_result.cpp:102-160` | `boolean/assembly.rs` | ☐ |
| `CreateEdges()` | `boolean_result.cpp:162-220` | `boolean/assembly.rs` | ☐ |
| `CreateFaces()` | `boolean_result.cpp:222-300` | `boolean/assembly.rs` | ☐ |
| `FilterEdges()` | `boolean_result.cpp:302-380` | `boolean/assembly.rs` | ☐ |
| `AssembleResult()` | `boolean_result.cpp:382-450` | `boolean/assembly.rs` | ☐ |

### 4.3 CSG Tree

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `CsgNode` struct | `csg_tree.h:40-80` | `boolean/csg.rs` | ☐ |
| `CsgLeafNode` struct | `csg_tree.h:82-120` | `boolean/csg.rs` | ☐ |
| `CsgOpNode` struct | `csg_tree.h:122-150` | `boolean/csg.rs` | ☐ |
| `CsgNode::Transform()` | `csg_tree.cpp:40-80` | `boolean/csg.rs` | ☐ |
| `CsgNode::GetLeafImpl()` | `csg_tree.cpp:82-130` | `boolean/csg.rs` | ☐ |
| Lazy evaluation | `csg_tree.cpp:132-250` | `boolean/csg.rs` | ☐ |

---

## Phase 5: Constructors & Hull (Weeks 17-19)

### 5.1 Primitive Constructors

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `Tetrahedron()` | `constructors.cpp:30-70` | `manifold/constructors.rs` | ☐ |
| `Cube()` | `constructors.cpp:72-130` | `manifold/constructors.rs` | ☐ |
| `Sphere()` | `constructors.cpp:132-200` | `manifold/constructors.rs` | ☐ |
| `Cylinder()` | `constructors.cpp:202-280` | `manifold/constructors.rs` | ☐ |
| `Extrude()` | `constructors.cpp:282-350` | `manifold/constructors.rs` | ☐ |
| `Revolve()` | `constructors.cpp:352-420` | `manifold/constructors.rs` | ☐ |
| `LevelSet()` | `constructors.cpp:422-500` | `manifold/constructors.rs` | ☐ |

### 5.2 Convex Hull (QuickHull)

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `QuickHull` struct | `quickhull.h:40-100` | `boolean/quickhull.rs` | ☐ |
| `Face` struct | `quickhull.h:102-140` | `boolean/quickhull.rs` | ☐ |
| `HalfEdge` struct | `quickhull.h:142-170` | `boolean/quickhull.rs` | ☐ |
| `buildConvexHull()` | `quickhull.cpp:40-200` | `boolean/quickhull.rs` | ☐ |
| `createInitialSimplex()` | `quickhull.cpp:202-280` | `boolean/quickhull.rs` | ☐ |
| `addPointToFace()` | `quickhull.cpp:282-350` | `boolean/quickhull.rs` | ☐ |
| `findHorizon()` | `quickhull.cpp:352-420` | `boolean/quickhull.rs` | ☐ |

---

## Phase 6: Polygon Triangulation (Weeks 20-22)

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `Tree2D` struct | `tree2d.h:30-80` | `polygon/tree2d.rs` | ☐ |
| `Tree2D::Insert()` | `tree2d.cpp:40-100` | `polygon/tree2d.rs` | ☐ |
| `Tree2D::Query()` | `tree2d.cpp:102-160` | `polygon/tree2d.rs` | ☐ |
| `Triangulate()` | `polygon.cpp:40-200` | `polygon/triangulate.rs` | ☐ |
| `EarClip()` | `polygon.cpp:202-300` | `polygon/triangulate.rs` | ☐ |
| `MonotoneTriangulate()` | `polygon.cpp:302-400` | `polygon/triangulate.rs` | ☐ |

---

## Phase 7: Advanced Features (Weeks 23-26)

### 7.1 Smoothing

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `CalculateTangents()` | `smoothing.cpp:40-120` | `boolean/smoothing.rs` | ☐ |
| `SmoothByNormals()` | `smoothing.cpp:122-200` | `boolean/smoothing.rs` | ☐ |
| `Smooth()` | `smoothing.cpp:202-300` | `boolean/smoothing.rs` | ☐ |
| Bezier subdivision | `smoothing.cpp:302-400` | `boolean/smoothing.rs` | ☐ |

### 7.2 Subdivision

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `Refine()` | `subdivision.cpp:40-150` | `boolean/subdivision.rs` | ☐ |
| `RefineToLength()` | `subdivision.cpp:152-250` | `boolean/subdivision.rs` | ☐ |
| `RefineToTolerance()` | `subdivision.cpp:252-350` | `boolean/subdivision.rs` | ☐ |

### 7.3 Properties

| Item | C++ Source | Rust File | Status |
|------|------------|-----------|--------|
| `InterpolateProp()` | `properties.cpp:40-120` | `boolean/properties.rs` | ☐ |
| `CalculateCurvature()` | `properties.cpp:122-200` | `boolean/properties.rs` | ☐ |
| `MinGap()` | `properties.cpp:202-300` | `boolean/properties.rs` | ☐ |

---

## Test Coverage Tracking

### C++ Test Porting Status

| Test File | Tests | Ported | Passing |
|-----------|-------|--------|---------|
| `boolean_test.cpp` | 39 | ☐ | ☐ |
| `boolean_complex_test.cpp` | 19 | ☐ | ☐ |
| `manifold_test.cpp` | 49 | ☐ | ☐ |
| `hull_test.cpp` | 9 | ☐ | ☐ |
| `polygon_test.cpp` | ~30 (file-based) | ☐ | ☐ |
| `smooth_test.cpp` | 13 | ☐ | ☐ |
| `sdf_test.cpp` | 9 | ☐ | ☐ |
| `cross_section_test.cpp` | 14 | ☐ | ☐ |
| `properties_test.cpp` | 22 | ☐ | ☐ |
| `samples_test.cpp` | 12 | ☐ | ☐ |
| `manifoldc_test.cpp` | 9 | ☐ | ☐ |
| `stl_intersection_test.cpp` | 2 | ☐ | ☐ |
| **Total** | **197+** | **0** | **0** |

### Golden Tests

| Category | Count | Generated | Passing |
|----------|-------|-----------|---------|
| Primitive shapes | 50 | ☐ | ☐ |
| Basic booleans | 100 | ☐ | ☐ |
| Complex booleans | 200 | ☐ | ☐ |
| Edge cases | 150 | ☐ | ☐ |
| Transformations | 100 | ☐ | ☐ |
| Mesh operations | 100 | ☐ | ☐ |
| Hull operations | 50 | ☐ | ☐ |
| Smoothing | 50 | ☐ | ☐ |
| Properties | 50 | ☐ | ☐ |
| **Total** | **850** | **0** | **0** |

---

## Verification Commands

```bash
# Run all tests
cargo test --workspace

# Run specific crate tests
cargo test -p manifold-types
cargo test -p manifold-math
cargo test -p manifold-parallel
cargo test -p manifold-collider
cargo test -p manifold-polygon
cargo test -p manifold-boolean
cargo test -p manifold

# Run golden tests
cargo test --test golden

# Run property tests (slow)
cargo test --test proptest -- --test-threads=1

# Run fuzz tests (requires nightly)
cargo +nightly fuzz run manifold_fuzz -- -max_total_time=3600
cargo +nightly fuzz run polygon_fuzz -- -max_total_time=3600

# Check for memory safety (requires Miri)
cargo +nightly miri test -p manifold-types
```

---

*Document Version: 1.0*
*Last Updated: 2026-01-04*
