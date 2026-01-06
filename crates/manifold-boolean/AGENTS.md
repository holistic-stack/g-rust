# AGENTS — manifold-boolean

## OVERVIEW
Core Boolean/kernel crate: CSG, smoothing, topology repair. Parity-critical with C++ `boolean3.cpp`/`impl.cpp`.

## STRUCTURE
```
crates/manifold-boolean/
├── src/
│   ├── kernel/         # Intersection kernels (01/11/12/03), smoothing, boolean_result
│   ├── impl/           # ManifoldImpl plumbing, edge ops stubs
│   ├── quickhull.rs    # Convex hull port
│   ├── helpers.rs      # Atomic math utils (intersect/interpolate)
│   └── tests/          # Tranche-based C++ test ports + determinism helpers
└── Cargo.toml
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Kernel intersections | src/kernel | Dimensional suffixes (01,11,12,03); manual monomorphization `_true_false` |
| Smoothing | src/kernel/smoothing.rs | Largest file (~960 lines); Bezier interp/tangent distribution |
| Boolean assembly | src/kernel/boolean_result.rs | Output sizing, edge stitching determinism |
| ManifoldImpl core | src/impl | Halfedge mesh ops; edge_ops stubs noted in specs |
| Helpers | src/helpers.rs | Shared atomic math; slated to move to math crate per coupling notes |
| Tests | src/tests | Tranche1/1b, util.rs (canonical signature, determinism), triangulation/determinism |

## CONVENTIONS
- Dimensional suffixes: 0=vertex,1=edge,2=face,3=volume.
- Manual monomorphization: `_true_false` variants for template bools; no runtime flags.
- Determinism: stable sort after parallel ops; prefer BTreeMap for ordered maps.
- Types mirror C++ (indices often i32, not usize) to ease parity.
- Coordinate swapping to reuse planar logic across axes.

## ANTI-PATTERNS
- No algorithmic “improvements” or reordered floating-point ops.
- Do not change parallel thresholds/policies (`kSeqThreshold=10000`).
- Avoid HashMap in parity-critical paths; keep order deterministic.
- No external geometry crates.

## NOTES
- Known stubs: impl/edge_ops.rs; quickhull completeness tracked in specs.
- Only root unsafe is in collider crate; this crate is safe—keep it so unless C++ forces unsafe.
- Reference comments map directly to `submodules/manifold` lines; consult before edits.
