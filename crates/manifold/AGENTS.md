# AGENTS — manifold

## OVERVIEW
Top-level public API crate. Wraps ManifoldImpl from boolean crate, exposes primitives (cube/sphere/etc.), transforms, and boolean operations.

## STRUCTURE
```
crates/manifold/
├── src/lib.rs        # Public API, re-exports, constructors
└── tests/            # API-level integration (boolean_tests, primitive_tests, math_tests)
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Public API surface | src/lib.rs | Wraps ManifoldImpl; constructors, boolean wrappers |
| Integration tests | tests/*.rs | Ported/ported-in-progress; currently failing, see diagnostics |

## CONVENTIONS
- Delegates heavy lifting to manifold-boolean; no business logic here.
- Constructors mirror C++ `constructors.cpp`; keep parameter parity.
- Error handling via thiserror; avoid panics unless C++ asserts.

## ANTI-PATTERNS
- Do not add new geometry features beyond C++ reference.
- Avoid duplicating math/helpers already in lower crates.

## NOTES
- Tests currently contain syntax/method-resolution errors; see diagnostics before enabling CI.
- Many methods expected by tests (get_mesh_gl, original_id, calculate_normals) are not implemented yet; check specs before stubbing.
