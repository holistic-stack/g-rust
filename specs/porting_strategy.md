# Manifold Rust Port: Strategy & Plan

## Core Philosophy
1.  **Parity First**: The primary goal is 1:1 algorithmic parity with the C++ reference.
2.  **Rust Idioms**: Use safe Rust where possible, but prioritize parity over "clean" Rust if it compromises the algorithm.
3.  **No New Features**: Do not add features or optimizations not present in the C++ version.

## Architecture
-   **Modules**: `manifold-boolean` is the core.
-   **Parallelism**: Use `rayon` to match TBB execution policies.
-   **Math**: Use `glam` for SIMD-accelerated linear algebra.

## Immediate Priorities
1.  **Fix Compilation**: Resolve name conflicts and type errors in `manifold-boolean`.
2.  **Stub Implementation**: Fill in critical stubs (`edge_ops.rs`, `quickhull.rs`) with functional logic.
3.  **Test Coverage**: Port C++ unit tests to verify parity.

## Dev Workflow
-   Create a plan for each module.
-   Implement logic incrementally.
-   Verify against C++ reference constantly.
