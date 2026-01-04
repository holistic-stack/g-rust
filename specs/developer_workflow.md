# Developer Workflow for Manifold C++ to Rust Migration

## Overview

This document defines the **mandatory workflow** that every developer must follow when porting C++ code from `submodules/manifold/` to Rust. Following this workflow ensures algorithm parity, correctness, and prevents deviation from the reference implementation.

> **CRITICAL**: No code may be merged without completing ALL workflow steps. No external algorithms or "improvements" are permitted.

---

## 1. Pre-Implementation Requirements

### 1.1 Before Starting ANY Implementation

Every developer **MUST** complete these steps before writing any Rust code:

| Step | Action | Verification |
|------|--------|--------------|
| 1 | Read the C++ source file completely | Sign-off in PR description |
| 2 | Document all functions to be ported | List in PR description |
| 3 | Identify all dependencies | Dependency checklist complete |
| 4 | Read related C++ tests | Test list in PR description |
| 5 | Check if dependent Rust code exists | Dependency check passed |

### 1.2 Mandatory Reading List

Before porting ANY file, developers must read and understand:

1. **This workflow document** (specs/developer_workflow.md)
2. **Critical algorithms** (specs/critical_algorithms.md)
3. **The specific C++ source file** in full
4. **The C++ header file** for any structs/classes
5. **Related C++ test file(s)**

---

## 2. Function-by-Function Porting Workflow

### 2.1 The Golden Rule

> **PORT CHARACTER-BY-CHARACTER. DO NOT "IMPROVE" OR "OPTIMIZE".**

### 2.2 Step-by-Step Process for Each Function

```
┌─────────────────────────────────────────────────────────────────┐
│  STEP 1: DOCUMENT                                               │
│  • Copy C++ function signature to comment                       │
│  • Note exact line numbers in source file                       │
│  • Document all parameters and return types                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 2: TRANSLATE SIGNATURE                                    │
│  • Convert C++ types to Rust equivalents                        │
│  • Preserve parameter order exactly                             │
│  • Use same naming (snake_case conversion only)                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 3: TRANSLATE BODY LINE-BY-LINE                            │
│  • One C++ line → One Rust statement (where possible)           │
│  • Same variable names (converted to Rust style)                │
│  • Same order of operations                                     │
│  • Same branching structure                                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 4: ADD REFERENCE COMMENT                                  │
│  • Document C++ source file and line numbers                    │
│  • Add original C++ code in doc comment if complex              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 5: WRITE UNIT TEST                                        │
│  • Port corresponding C++ test first                            │
│  • Test must pass before continuing                             │
│  • Add edge case tests                                          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 6: VERIFY                                                 │
│  • Compare output with C++ for same input                       │
│  • Check floating-point results bit-for-bit if critical         │
│  • Run regression tests                                         │
└─────────────────────────────────────────────────────────────────┘
```

### 2.3 Function Documentation Template

Every ported function MUST have this documentation format:

```rust
/// Port of C++ function `FunctionName` from `source_file.cpp` lines X-Y.
///
/// # C++ Reference
/// ```cpp
/// // Original C++ signature (copy exactly)
/// ReturnType FunctionName(Param1 p1, Param2 p2);
/// ```
///
/// # Algorithm
/// Brief description of what the function does (from C++ comments or understanding)
///
/// # Parameters
/// * `param1` - Description (match C++ parameter name)
/// * `param2` - Description
///
/// # Returns
/// Description of return value
///
/// # Panics
/// Conditions under which this function panics (matching C++ assertions)
///
/// # Safety
/// If unsafe, explain why and what invariants must be maintained
pub fn function_name(param1: Type1, param2: Type2) -> ReturnType {
    // Implementation
}
```

---

## 3. Type Translation Rules

### 3.1 Mandatory Type Mappings

| C++ Type | Rust Type | Notes |
|----------|-----------|-------|
| `double` | `f64` | NEVER use f32 |
| `float` | `f32` | Only if C++ uses float |
| `int` | `i32` | Match signedness |
| `int32_t` | `i32` | Exact match |
| `int64_t` | `i64` | Exact match |
| `uint32_t` | `u32` | Exact match |
| `uint64_t` | `u64` | Exact match |
| `size_t` | `usize` | Platform-dependent OK |
| `bool` | `bool` | Exact match |
| `vec2` | `glam::DVec2` | Double precision |
| `vec3` | `glam::DVec3` | Double precision |
| `vec4` | `glam::DVec4` | Double precision |
| `mat3` | `glam::DMat3` | Double precision |
| `mat4` | `glam::DMat4` | Double precision |
| `mat3x4` | Custom struct | Port from linalg.h |
| `std::vector<T>` | `Vec<T>` | Capacity behavior may differ |
| `std::array<T,N>` | `[T; N]` | Exact match |
| `std::pair<A,B>` | `(A, B)` | Tuple |
| `std::optional<T>` | `Option<T>` | Exact match |
| `std::map<K,V>` | `BTreeMap<K,V>` | MUST be ordered |
| `std::unordered_map<K,V>` | `HashMap<K,V>` | Only if order doesn't matter |
| `std::shared_ptr<T>` | `Arc<T>` | Thread-safe |
| `std::unique_ptr<T>` | `Box<T>` | Owned |
| `T*` (non-owning) | `&T` or `&mut T` | Borrow |
| `T&` | `&T` | Immutable reference |
| `const T&` | `&T` | Immutable reference |

### 3.2 Template Translation Rules

| C++ Pattern | Rust Pattern |
|-------------|--------------|
| `template<typename T>` | Generic `<T>` |
| `template<bool B>` | Const generic `<const B: bool>` |
| `template<int N>` | Const generic `<const N: i32>` |
| SFINAE | Trait bounds |
| `std::enable_if` | `where` clauses |

### 3.3 Parallel Primitive Mappings

| C++ (TBB) | Rust (Rayon) | Notes |
|-----------|--------------|-------|
| `tbb::parallel_for` | `par_iter().for_each()` | With policy check |
| `tbb::parallel_reduce` | `par_iter().reduce()` | With policy check |
| `tbb::parallel_sort` | `par_sort()` | Stable sort required |
| `tbb::parallel_invoke` | `rayon::join()` | Fork-join |
| Sequential fallback | Use policy enum | Match threshold |

---

## 4. Acceptance Criteria Per Phase

### 4.1 Phase 1: Foundation Types (manifold-types, manifold-math)

| Criterion | Verification Method | Pass Condition |
|-----------|--------------------|--------------------|
| All structs match C++ | Side-by-side comparison | Field names, types, order match |
| All methods implemented | Checklist review | 100% coverage |
| Unit tests pass | `cargo test` | All green |
| No external dependencies | `cargo tree` inspection | Only glam, no algorithm crates |

**Sign-off required from:** Tech Lead

### 4.2 Phase 2: Parallel Primitives (manifold-parallel)

| Criterion | Verification Method | Pass Condition |
|-----------|--------------------|--------------------|
| All functions from parallel.h | Line-by-line comparison | 100% coverage |
| Threshold matches | Constant comparison | K_SEQ_THRESHOLD = 10000 |
| Policy logic matches | Code review | Exact branch structure |
| Deterministic output | Golden tests | Bit-identical results |

**Sign-off required from:** Tech Lead + Reviewer

### 4.3 Phase 3: Spatial Acceleration (manifold-collider)

| Criterion | Verification Method | Pass Condition |
|-----------|--------------------|--------------------|
| BVH construction matches | Golden tests | Same tree structure |
| Collision detection matches | Golden tests | Same collision pairs |
| Morton codes match | Unit tests | Bit-identical |

**Sign-off required from:** Tech Lead + Reviewer

### 4.4 Phase 4: Boolean Kernel (manifold-boolean) - CRITICAL

| Criterion | Verification Method | Pass Condition |
|-----------|--------------------|--------------------|
| `Interpolate` exact | Golden test + bit comparison | Bit-identical for 1000+ cases |
| `Intersect` exact | Golden test + bit comparison | Bit-identical for 1000+ cases |
| `Shadows` exact | Unit tests | Exact equality semantics |
| All Kernel functions | Code review line-by-line | Character-by-character match |
| All 31 boolean_test.cpp | `cargo test` | All pass |
| All 22 boolean_complex_test.cpp | `cargo test` | All pass |
| Symbolic perturbation | Stress tests | No failures |

**Sign-off required from:** Tech Lead + 2 Reviewers + Golden test suite

### 4.5 Phase 5-7: Constructors, Hull, Polygon, Advanced

| Criterion | Verification Method | Pass Condition |
|-----------|--------------------|--------------------|
| All C++ tests ported | Test checklist | 100% coverage |
| All tests pass | `cargo test` | All green |
| Golden tests match | Golden test suite | 100% match |

**Sign-off required from:** Tech Lead + Reviewer

### 4.6 Phase 8-9: API, Verification, Release

| Criterion | Verification Method | Pass Condition |
|-----------|--------------------|--------------------|
| All 197 C++ tests pass | `cargo test` | 100% pass |
| All 850+ golden tests pass | Golden test suite | 100% match |
| Fuzz testing | 1 week continuous | No crashes |
| Property tests | 1M iterations | All pass |
| Performance | Benchmarks | Within 2x of C++ |
| Documentation | Doc coverage | 100% public items |

**Final sign-off required from:** Project Lead + All Tech Leads

---

## 5. Checks and Balances

### 5.1 Automated Checks (CI/CD)

Every PR must pass these automated checks:

```yaml
# .github/workflows/ci.yml (conceptual)
checks:
  - name: "Build"
    command: "cargo build --workspace"
    
  - name: "Unit Tests"
    command: "cargo test --workspace"
    
  - name: "Clippy Lints"
    command: "cargo clippy --workspace -- -D warnings"
    
  - name: "Format Check"
    command: "cargo fmt --check"
    
  - name: "No External Algorithms"
    command: "./scripts/check_no_external_algorithms.sh"
    
  - name: "Golden Tests"
    command: "cargo test --test golden"
    
  - name: "Documentation Coverage"
    command: "cargo doc --no-deps"
```

### 5.2 Manual Review Checklist

Every PR review MUST verify:

- [ ] C++ source file and line numbers documented
- [ ] Function signature matches C++ exactly
- [ ] Variable names match (with Rust style conversion)
- [ ] Order of operations preserved
- [ ] Branching structure preserved
- [ ] No "optimizations" or "improvements" added
- [ ] No external algorithm crates used
- [ ] Related C++ test ported
- [ ] Test passes
- [ ] Golden test comparison (if applicable)

### 5.3 Code Review Requirements

| Change Type | Minimum Reviewers | Special Requirements |
|-------------|-------------------|----------------------|
| Types/Math | 1 | Type comparison sign-off |
| Parallel primitives | 2 | Determinism verification |
| Boolean kernel | 2 + Golden tests | Line-by-line review |
| Other | 1 | Standard review |

---

## 6. Guardrails

### 6.1 Prohibited Actions

The following are **STRICTLY PROHIBITED**:

1. **NO external geometry algorithms**
   - Do NOT use `geo`, `cgmath`, `nalgebra` for algorithms
   - Do NOT use `parry`, `ncollide` or any collision library
   - Do NOT use any triangulation library except for reference comparison
   
2. **NO "improvements" to C++ algorithms**
   - Do NOT optimize loop structures
   - Do NOT change data structure choices
   - Do NOT add SIMD unless C++ uses it
   
3. **NO floating-point "fixes"**
   - Do NOT add epsilon comparisons where C++ uses exact equality
   - Do NOT change order of operations
   - Do NOT use fast-math optimizations
   
4. **NO changing parallelism behavior**
   - Do NOT change threshold values
   - Do NOT change execution policy logic
   - Do NOT use different parallel constructs

### 6.2 Automated Guardrail Script

```bash
#!/bin/bash
# scripts/check_no_external_algorithms.sh

FORBIDDEN_CRATES=(
    "geo"
    "cgmath"
    "parry"
    "ncollide"
    "delaunator"
    "spade"
    "robust"
    "predicates"
)

for crate in "${FORBIDDEN_CRATES[@]}"; do
    if grep -r "use $crate" crates/; then
        echo "ERROR: Forbidden crate '$crate' found!"
        exit 1
    fi
    if grep "$crate" Cargo.toml Cargo.lock 2>/dev/null; then
        echo "ERROR: Forbidden crate '$crate' in dependencies!"
        exit 1
    fi
done

echo "OK: No forbidden crates found"
```

### 6.3 Dependency Whitelist

Only these external crates are permitted:

| Crate | Purpose | Version |
|-------|---------|---------|
| `glam` | Linear algebra types | Latest |
| `rayon` | Parallelism | Latest |
| `thiserror` | Error types | Latest |
| `serde` | Serialization (optional) | Latest |
| `proptest` | Property testing (dev) | Latest |
| `criterion` | Benchmarking (dev) | Latest |
| `rstest` | Test fixtures (dev) | Latest |

Any additional dependencies require **Project Lead approval**.

---

## 7. Testing Requirements

### 7.1 Test Pyramid

```
                    ┌───────────────────┐
                    │   Golden Tests    │  ← 850+ comparisons with C++
                    │   (Parity)        │
                    └─────────┬─────────┘
                              │
                    ┌─────────▼─────────┐
                    │  Integration      │  ← 197 C++ tests ported
                    │  Tests            │
                    └─────────┬─────────┘
                              │
            ┌─────────────────▼─────────────────┐
            │         Unit Tests                │  ← Every function tested
            │         (Function level)          │
            └─────────────────┬─────────────────┘
                              │
    ┌─────────────────────────▼─────────────────────────┐
    │              Property Tests + Fuzz Tests          │  ← 1M+ iterations
    │              (Invariant verification)             │
    └───────────────────────────────────────────────────┘
```

### 7.2 Test Coverage Requirements

| Category | Minimum Coverage |
|----------|-----------------|
| manifold-types | 100% |
| manifold-math | 100% |
| manifold-parallel | 95% |
| manifold-collider | 95% |
| manifold-boolean | 100% |
| manifold-polygon | 95% |
| manifold | 90% |

### 7.3 Golden Test Generation

Golden tests are generated by running C++ with specific inputs and capturing outputs:

```bash
# Generate golden test data from C++ implementation
cd submodules/manifold/build
./test/manifold_test --gtest_output=json:golden_tests.json

# Convert to Rust test format
python scripts/convert_golden_tests.py golden_tests.json > tests/golden/data.rs
```

### 7.4 Mandatory Test Categories

Every module must have:

1. **Basic functionality tests** - Simple cases
2. **Edge case tests** - Boundary conditions, empty inputs
3. **Error case tests** - Invalid inputs, expected failures
4. **Regression tests** - Bugs found during development
5. **Golden tests** - Comparison with C++ output
6. **Stress tests** - Large inputs, many operations

---

## 8. Debugging and Verification Tools

### 8.1 Comparing Rust and C++ Output

```rust
// Debug helper for comparing floating-point arrays
fn compare_arrays(rust: &[f64], cpp: &[f64], name: &str) {
    assert_eq!(rust.len(), cpp.len(), "{} length mismatch", name);
    for (i, (r, c)) in rust.iter().zip(cpp.iter()).enumerate() {
        if (r - c).abs() > 1e-15 {
            eprintln!("MISMATCH {name}[{i}]: rust={r}, cpp={c}, diff={}", r - c);
        }
    }
}
```

### 8.2 Debug Build Features

```toml
# Cargo.toml
[features]
debug-compare = []  # Enable C++ comparison output
trace-ops = []      # Trace all operations
verbose-assert = [] # Extra assertion messages
```

### 8.3 Verification Script

```bash
#!/bin/bash
# scripts/verify_parity.sh

# Run both C++ and Rust with same input
./cpp_test_runner $INPUT > cpp_output.json
./rust_test_runner $INPUT > rust_output.json

# Compare outputs
python scripts/compare_outputs.py cpp_output.json rust_output.json
```

---

## 9. Pull Request Template

```markdown
## C++ Source Reference
- **File**: `src/filename.cpp`
- **Lines**: X-Y
- **Functions ported**: 
  - `FunctionA` (lines X-Y)
  - `FunctionB` (lines X-Y)

## Pre-Implementation Checklist
- [ ] Read entire C++ source file
- [ ] Read related C++ tests
- [ ] Verified dependencies are available
- [ ] Updated implementation_checklist.md

## Implementation Checklist
- [ ] All functions ported character-by-character
- [ ] All doc comments reference C++ source
- [ ] No external algorithms used
- [ ] No "improvements" made

## Testing Checklist
- [ ] Related C++ tests ported
- [ ] All unit tests pass
- [ ] Golden tests pass (if applicable)
- [ ] Edge cases covered

## Review Checklist (for reviewers)
- [ ] Verified C++ line numbers are correct
- [ ] Verified function signatures match
- [ ] Verified order of operations matches
- [ ] Verified no external algorithms
- [ ] Verified tests are comprehensive
```

---

## 10. Escalation Path

### 10.1 When to Escalate

Escalate to Tech Lead when:
- C++ code is ambiguous or undocumented
- C++ test doesn't exist for a function
- Rust equivalent behavior is unclear
- Performance differs significantly from C++

### 10.2 Escalation Process

1. Document the issue in detail
2. Show C++ code and attempted Rust port
3. Explain specific uncertainty
4. Wait for guidance before proceeding

### 10.3 Decision Authority

| Decision Type | Authority |
|---------------|-----------|
| Type mapping questions | Tech Lead |
| Algorithm ambiguity | Tech Lead + C++ author (if available) |
| Dependency additions | Project Lead |
| Architecture changes | Project Lead + All Tech Leads |

---

## 11. Progress Tracking

### 11.1 Daily Standup Items

- Functions ported yesterday
- Functions planned today
- Blockers (C++ ambiguity, test failures, etc.)

### 11.2 Weekly Review

- Phase progress update
- Test coverage metrics
- Golden test pass rate
- Blockers and escalations

### 11.3 Phase Completion

A phase is complete only when:
- [ ] All items in implementation_checklist.md marked done
- [ ] All related tests pass
- [ ] Code review approved
- [ ] Sign-off received from required parties

---

*Document Version: 1.0*
*Last Updated: 2026-01-04*
