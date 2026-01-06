# PROJECT KNOWLEDGE BASE

**Generated:** 2026-01-06
**Commit:** 0348b58
**Branch:** main

## OVERVIEW
Rust workspace port of the Manifold geometry library. Goal: strict 1:1 algorithmic parity with upstream C++ (submodules/manifold). Heavy focus on deterministic boolean kernels and testing.

## STRUCTURE
```
./
├── crates/             # Workspace crates (boolean kernel, math, types, API, etc.)
├── specs/              # Parity rules, file mapping, migration plans
├── submodules/manifold # C++ reference implementation (read-only)
├── tests/              # Golden JSON manifests + integration tests
└── target/             # Build artifacts (ignore)
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Public API & constructors | crates/manifold | High-level `Manifold` wrapper and primitives |
| Boolean kernel | crates/manifold-boolean/src/kernel | Core CSG + smoothing; parity-critical |
| Shared math/types | crates/manifold-math, crates/manifold-types | Deterministic math, Halfedge mesh basics |
| Parallel primitives | crates/manifold-parallel | Rayon wrappers matching TBB policies |
| Collider/BVH | crates/manifold-collider | AABB tree, single `unsafe` block for slice aliasing |
| Polygon/triangulation | crates/manifold-polygon | 2D cross-section helpers |
| Porting specs | specs/*.md | Golden rules, file mapping, roadmap |
| Golden data | tests/golden | C++ GTest JSON references |

## CONVENTIONS (PROJECT-SPECIFIC)
- Parity first: line-for-line C++ ports; no “improvements”.
- Determinism: use `BTreeMap` where C++ uses `std::map`; stable sorts after parallel work.
- Numeric: `f64` only where C++ uses double; no epsilon unless present upstream.
- Parallel thresholds: respect C++ constants (e.g., `kSeqThreshold=10000`).
- Dependencies: only glam, rayon, serde, thiserror, proptest, rstest, criterion; geometry crates forbidden.

## ANTI-PATTERNS (DO NOT)
- Do not swap algorithms or reorder floating-point ops.
- Do not introduce new deps (geo/nalgebra/parry/etc.).
- Do not change execution policy or thresholds without C++ proof.
- Avoid refactors while fixing bugs; keep minimal diffs.

## UNIQUE STYLES
- Dimensional suffixing in kernel names (`01`, `12`, `03`) mirrors C++ template arity.
- Manual monomorphization via `_true_false` style flags for hot paths.
- Determinism helpers: canonical signature hashing in `crates/manifold-boolean/src/tests/util.rs`.

## COMMANDS
```bash
cargo fmt
cargo clippy --workspace
cargo test --workspace
# Boolean-only: cargo test -p manifold-boolean
```

## NOTES
- C++ reference lives in submodules/manifold; consult when porting.
- Specs folder is authoritative for file mapping and critical algorithms.
- Only known `unsafe` is in `crates/manifold-collider/src/lib.rs` (slice view in parallel loop).
