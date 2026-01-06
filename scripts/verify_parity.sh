#!/bin/bash
set -e

echo "Running parity checks..."

# 1. Format Check
echo "Checking format..."
cargo fmt --check

# 2. Clippy Check
echo "Running clippy..."
cargo clippy --workspace -- -D warnings

# 3. Compile Check
echo "Checking compilation..."
cargo check --workspace

# 4. Forbidden Patterns
echo "Checking for forbidden patterns..."
if grep -r "HashMap" crates/manifold-boolean/src/kernel/boolean_result.rs; then
    echo "Error: HashMap found in boolean_result.rs. Use BTreeMap for determinism."
    exit 1
fi

if grep -r "HashMap" crates/manifold-boolean/src/kernel/mod.rs; then
    # We allow it if it's explicitly justified, but ideally we want BTreeMap in MeshRelationD
    echo "Warning: HashMap found in kernel/mod.rs. Verify it's not affecting MeshRelationD output order."
fi

# Check for massive doc blocks in quickhull (heuristic)
if [ -f crates/manifold-boolean/src/kernel/quickhull.rs ]; then
    DOC_LINES=$(grep "///" crates/manifold-boolean/src/kernel/quickhull.rs | wc -l)
    if [ "$DOC_LINES" -gt 100 ]; then
        echo "Error: quickhull.rs seems to have excessive documentation lines ($DOC_LINES). Replace with pointers to C++."
        exit 1
    fi
fi

echo "Parity checks passed!"
