# AGENTS — manifold-types

## OVERVIEW
Shared geometry types crate. Holds Halfedge mesh primitives, TriRef, OpType, error types. Used by all other crates.

## STRUCTURE
```
crates/manifold-types/
├── src/lib.rs     # Re-exports; type definitions
└── src/*          # Modules: halfedge, tri_ref, op_type, bbox, errors
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Halfedge mesh types | src/halfedge.rs | Core mesh connectivity
| Triangle references | src/tri_ref.rs | Tri indexing helpers
| Operation metadata | src/op_type.rs  | OpType enum (Add/Sub/Intersect)
| Bounding boxes | src/bbox.rs       | AABB helpers
| Errors | src/errors.rs            | thiserror-based

## CONVENTIONS
- Types mirror C++ `shared.h`; keep field order and signedness.
- Prefer `i32` for indices when matching C++.
- Deterministic maps: use BTreeMap when porting std::map.

## ANTI-PATTERNS
- Do not add new geometry features or alter field ordering.
- Avoid HashMap where order matters.

## NOTES
- Crate is dependency root for most workspace crates; changing types has wide blast radius.
