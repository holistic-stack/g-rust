# Manifold C++ to Rust Migration Plan

> **Note**: This document provides an overview. For detailed specifications, see:
> - [Comprehensive Migration Plan](./comprehensive_migration_plan.md) - Full technical details
> - [Test Specification](./test_specification.md) - Complete test inventory
> - [File Mapping](./file_mapping.md) - Source-to-destination mapping

## 1. Introduction & Goals
The goal of this project is to perform a **strict 1:1 migration** of the Manifold geometry library from C++ to Rust. 

### 1.1 Principles
*   **Functional Parity**: The Rust implementation must produce bit-for-bit (or epsilon-equivalent) results to the C++ implementation.
*   **Strict Algorithm Mapping**: No alternative algorithms or "better" approaches are allowed. Every branch, loop, and predicate must match the C++ reference.
*   **Numerical Determinism**: Special mathematical functions (e.g., `sun_acos` in `impl.cpp`) must be ported exactly to maintain normal parity. Standard library versions are only allowed if they are functionally identical to the C++ reference's implementation.
*   **Parallelism Parity**: Parallel execution must respect the same thresholds and policies (e.g., `kSeqThreshold = 1e4`) defined in `parallel.h`.
*   **Safety**: Eliminate memory safety issues inherent in C++ while maintaining performance and algorithmic identity.
*   **CRITICAL: C++ Reference Only**: The ONLY authoritative source is `submodules/manifold/`. No external algorithms, libraries, or "improvements" are permitted.

## 2. Architecture Mapping

### 2.1 Data Structures
| C++ Structure | Rust Equivalent | Notes |
| :--- | :--- | :--- |
| `Manifold` | `manifold::Manifold` | Public handle, wraps `Arc<CsgNode>`. |
| `Manifold::Impl` | `manifold_boolean::Impl` | Internal B-Rep representation. |
| `CsgNode` | `manifold_boolean::CsgNode` | Enum or Trait for lazy evaluation tree. |
| `Halfedge` | `manifold_types::Halfedge` | Core topology unit. |
| `vec3` | `glam::DVec3` | Double precision 3D vector. |
| `mat3x4` | `glam::DMat4` (or custom 3x4) | Transformation matrix. |
| `Box` | `manifold_math::Box` | Axis-aligned bounding box. |
| `Vec<T>` | `Vec<T>` | Standard vector. Note: uninitialized resize in C++ must be handled with care in Rust (e.g., `Vec::with_capacity` + `unsafe` or `MaybeUninit` if critical for perf, otherwise zero-init is acceptable if bit-parity is preserved). |
| `VecView<T>` | `&[T]` or `&mut [T]` | Rust slices are native views. |
| `Collider` | `manifold_boolean::Collider` | Spatial acceleration (BVH). |
| `MeshGLP<Precision, I>` | `manifold_types::MeshGeneric<P, I>` | Generic mesh for IO. Match C++ aliases (MeshGL, MeshGL64). |

### 2.2 Numerical Determinism & Deterministic Predicates
*   **`Shadows` logic**: The predicate in `boolean3.cpp` must be ported without simplification. 
*   **`Interpolate` and `Intersect`**: These are the *only* places for floating point operations in the Boolean kernel. No other math libraries (like `robust`) should be introduced here.
*   **`sun_acos`**: The specific `acos` implementation in `impl.cpp` (derived from FreeBSD msun) must be ported to Rust to ensure vertex normals match exactly.

## 3. Crate Structure
To ensure modularity and clear dependency management:

1.  **`manifold-types`**: Core primitive types (`Halfedge`, `TriRef`, `MeshRelation`, etc.).
2.  **`manifold-math`**: Geometric primitives (`Box`, `Ray`, `Plane`), constants, and math utilities.
3.  **`manifold-boolean`**: The core boolean kernel (Winding numbers, Symbolic perturbation, Intersections, BVH).
4.  **`manifold-polygon`**: 2D polygon triangulation and utilities (port of `polygon.cpp` and `tree2d.h`).
5.  **`manifold-meshio`**: Importers/Exporters (port of `meshIO.cpp`, using Rust equivalents for Assimp like `gltf`).
6.  **`manifold`**: The main public API (port of `manifold.cpp`).

## 4. Core Algorithm Migration

### 4.1 Boolean Pipeline (`boolean3.cpp`)
The core of the library is the `Boolean` operation. It must be ported in these stages:
1.  **Intersections**: Port `Intersect12` (Edge-Face intersections).
2.  **Winding Numbers**: Port `Winding03` (Vertex winding numbers).
3.  **Symbolic Perturbation**: Port the logic that handles degenerate cases by infinitesimally expanding/contracting meshes.
4.  **Assembly**: Port the logic that combines the results into a new manifold mesh.

### 4.2 Symbolic Perturbation logic
The `Shadows` function is the critical tie-breaker:
```rust
fn shadows(p: f64, q: f64, dir: f64) -> i32 {
    if p == q {
        if dir < 0.0 { 1 } else { 0 }
    } else if p < q {
        1
    } else {
        0
    }
}
```
This logic must be preserved exactly.

### 4.3 Parallelism
Replace Intel TBB with **Rayon**, but enforce C++ execution policies.
*   **Thresholding**: Implement a wrapper (e.g., `par_if`) that respects `kSeqThreshold = 1e4`.
*   **Operations**:
    *   `tbb::parallel_for` $\rightarrow$ `rayon::iter::IntoParallelIterator::for_each` (if length > threshold).
    *   `tbb::parallel_reduce` $\rightarrow$ `rayon::iter::ParallelIterator::reduce`.
    *   `tbb::parallel_scan` $\rightarrow$ `rayon::iter::ParallelIterator::scan` (careful mapping required).
*   **Ordering**: C++'s `std::map` is ordered. If the C++ code relies on this ordering for determinism (common in geometry), use `BTreeMap` in Rust, not `HashMap`.

### 4.4 Data Structures & Atomics
*   **Disjoint Sets**: Match the thread-safe implementation in `disjoint_sets.h`.
*   **Atomic Counters**: Use `std::sync::atomic` for `meshIDCounter_` and other global counters to match C++ thread-safety exactly.
*   **Sorting**: Use the same Morton code logic in `sort.cpp` for spatial partitioning.

## 5. Testing Strategy (Comprehensive)

The testing strategy follows a multi-level approach to ensure correctness from basic primitives to complex topological edge cases.

> **See [Test Specification](./test_specification.md) for the complete test inventory with 500+ tests.**

### 5.1 Hierarchical Test Levels
1.  **Level 0: Atomic Units**: Test individual math functions (`Interpolate`, `Intersect`, `Shadows`) with edge cases (NaN, Infinity, identical points).
2.  **Level 1: Geometric Primitives**: Validate that `Cube`, `Sphere`, `Cylinder`, and `Tetrahedron` constructors produce valid, closed, and manifold meshes.
3.  **Level 2: Basic Booleans**: Test Union, Difference, and Intersection of two simple primitives in various relative positions:
    *   Disjoint (no contact).
    *   Overlapping.
    *   Fully contained.
    *   Touching at a vertex, edge, or face (coplanar).
4.  **Level 3: Complex CSG**:
    *   Chained operations: `(A + B) - C ^ D`.
    *   Deeply nested trees (depth 10+).
    *   "Stress" scenarios: Union of 100+ overlapping spheres to test robustness of winding number accumulation.
5.  **Level 4: Topological Stress**:
    *   Self-intersecting input handling (generalized winding numbers).
    *   Meshes with extremely small or large faces (epsilon robustness).
    *   "Sliver" triangles and degenerate geometry.

### 5.2 C++ Test Files to Port (ALL REQUIRED)

| C++ Test File | Tests | Priority |
|---------------|-------|----------|
| `boolean_test.cpp` | 35+ | P0 |
| `boolean_complex_test.cpp` | 40+ | P0 |
| `manifold_test.cpp` | 50+ | P0 |
| `hull_test.cpp` | 15+ | P1 |
| `polygon_test.cpp` | 30+ (from corpus) | P1 |
| `smooth_test.cpp` | 20+ | P2 |
| `sdf_test.cpp` | 10+ | P2 |
| `cross_section_test.cpp` | 25+ | P2 |
| `properties_test.cpp` | 15+ | P2 |
| `samples_test.cpp` | 10+ | P3 |
| `manifold_fuzz.cpp` | Fuzz target | P1 |
| `polygon_fuzz.cpp` | Fuzz target | P1 |

### 5.2 Parity Verification (Golden Tests)
*   **Bit-for-Bit Comparison**: For a set of 1000 standard operations, the Rust output must match the C++ output topology exactly.
*   **State Inspection**: Implement a debug feature in both C++ and Rust to export intermediate pipeline states (e.g., intersection lists, winding number volumes) to JSON for automated comparison.

### 5.3 Fuzzing and Property-Based Testing
*   **`proptest`**: Use property-based testing in Rust to verify that `Manifold::is_manifold()` always returns true for any combination of boolean operations.
*   **`cargo-fuzz`**: Continuous fuzzing of the boolean kernel using random geometric transformations and operations.

## 6. Detailed File Mapping

| C++ File | Rust Crate/Module | Status |
| :--- | :--- | :--- |
| `common.h` | `manifold-math::common` | |
| `vec.h`, `vec_view.h` | Native Rust `Vec` and slices | |
| `linalg.h` | `glam` dependency | |
| `impl.h/cpp` | `manifold-boolean::impl` | |
| `boolean3.h/cpp` | `manifold-boolean::kernel` | |
| `polygon.h/cpp` | `manifold-polygon` | |
| `csg_tree.h/cpp` | `manifold-boolean::csg_tree` | |
| `quickhull.h/cpp` | `manifold-boolean::quickhull` | |
| `collider.h` | `manifold-boolean::collider` | |
| `manifold.h/cpp` | `manifold` (top-level) | |

## 7. Migration Roadmap

### Phase 1: Foundation (Weeks 1-2)
*   Setup workspace and crates.
*   Port `manifold-math` and `manifold-types`.
*   Port `collider.h` (BVH).

### Phase 2: Geometry Core (Weeks 3-6)
*   Port `quickhull`.
*   Port mesh implementation (`impl.cpp`).
*   Port primitive constructors (Cube, Sphere, etc.).

### Phase 3: Boolean Kernel (Weeks 7-12)
*   Port intersection logic.
*   Port winding number logic.
*   Port symbolic perturbation.
*   Port assembly logic.

### Phase 4: Refinement & Advanced (Weeks 13-16)
*   Port triangulation (`polygon.cpp`).
*   Port smoothing and subdivision.
*   Port property interpolation.

### Phase 5: Verification (Weeks 17-20)
*   Enable Rayon parallelism.
*   Run comprehensive golden tests.
*   Address performance bottlenecks.

## 8. Risk Mitigation

*   **Floating Point Determinism**: Ensure `f64` operations match C++ `double` exactly. Be wary of compiler-specific optimizations (e.g., FMA) that might differ between GCC/Clang and Rustc.
*   **Safety vs. Parity**: While Rust's safety is a goal, **algorithmic parity is absolute**. If C++ uses `unsafe`-equivalent patterns (like manual memory management for performance), those must be ported to equivalent (safe or unsafe) Rust that produces the *exact same results*.
*   **External Dependencies**: The 2D module depends on **Clipper2**. To maintain parity, we should use a Rust port of Clipper2 that is verified against the version used in C++.
*   **Triangulation Robustness**: Triangulation is notoriously difficult. Use the exact same ear-clipping or sweep-line algorithm as in `polygon.cpp`.
