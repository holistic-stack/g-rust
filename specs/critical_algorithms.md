# Critical Algorithm Implementations

This document specifies the exact algorithms that **MUST** be ported character-by-character from the C++ reference. No modifications, optimizations, or alternative implementations are permitted.

---

## 1. Boolean Kernel Functions (boolean3.cpp)

### 1.1 Interpolate Function (Lines 44-56)

**Purpose**: Linear interpolation along x-axis, choosing the endpoint closer to the target.

**C++ Reference**:
```cpp
inline vec2 Interpolate(vec3 aL, vec3 aR, double x) {
  const double dxL = x - aL.x;
  const double dxR = x - aR.x;
  DEBUG_ASSERT(dxL * dxR <= 0, logicErr,
               "Boolean manifold error: not in domain");
  const bool useL = fabs(dxL) < fabs(dxR);
  const vec3 dLR = aR - aL;
  const double lambda = (useL ? dxL : dxR) / dLR.x;
  if (!std::isfinite(lambda) || !std::isfinite(dLR.y) || !std::isfinite(dLR.z))
    return vec2(aL.y, aL.z);
  vec2 yz;
  yz[0] = lambda * dLR.y + (useL ? aL.y : aR.y);
  yz[1] = lambda * dLR.z + (useL ? aL.z : aR.z);
  return yz;
}
```

**Rust Port** (EXACT):
```rust
#[inline]
pub fn interpolate(a_l: DVec3, a_r: DVec3, x: f64) -> DVec2 {
    let dx_l = x - a_l.x;
    let dx_r = x - a_r.x;
    debug_assert!(
        dx_l * dx_r <= 0.0,
        "Boolean manifold error: not in domain"
    );
    let use_l = dx_l.abs() < dx_r.abs();
    let d_lr = a_r - a_l;
    let lambda = if use_l { dx_l } else { dx_r } / d_lr.x;
    if !lambda.is_finite() || !d_lr.y.is_finite() || !d_lr.z.is_finite() {
        return DVec2::new(a_l.y, a_l.z);
    }
    DVec2::new(
        lambda * d_lr.y + if use_l { a_l.y } else { a_r.y },
        lambda * d_lr.z + if use_l { a_l.z } else { a_r.z },
    )
}
```

### 1.2 Intersect Function (Lines 58-79)

**Purpose**: Find intersection point of two line segments in 3D projected to y-axis.

**C++ Reference**:
```cpp
vec4 Intersect(const vec3& aL, const vec3& aR, const vec3& bL, const vec3& bR) {
  const double dyL = bL.y - aL.y;
  const double dyR = bR.y - aR.y;
  DEBUG_ASSERT(dyL * dyR <= 0, logicErr,
               "Boolean manifold error: no intersection");
  const bool useL = fabs(dyL) < fabs(dyR);
  const double dx = aR.x - aL.x;
  double lambda = (useL ? dyL : dyR) / (dyL - dyR);
  if (!std::isfinite(lambda)) lambda = 0.0;
  vec4 xyzz;
  xyzz.x = lambda * dx + (useL ? aL.x : aR.x);
  const double aDy = aR.y - aL.y;
  const double bDy = bR.y - bL.y;
  const bool useA = fabs(aDy) < fabs(bDy);
  xyzz.y = lambda * (useA ? aDy : bDy) +
           (useL ? (useA ? aL.y : bL.y) : (useA ? aR.y : bR.y));
  xyzz.z = lambda * (aR.z - aL.z) + (useL ? aL.z : aR.z);
  xyzz.w = lambda * (bR.z - bL.z) + (useL ? bL.z : bR.z);
  return xyzz;
}
```

**Rust Port** (EXACT):
```rust
pub fn intersect(a_l: DVec3, a_r: DVec3, b_l: DVec3, b_r: DVec3) -> DVec4 {
    let dy_l = b_l.y - a_l.y;
    let dy_r = b_r.y - a_r.y;
    debug_assert!(
        dy_l * dy_r <= 0.0,
        "Boolean manifold error: no intersection"
    );
    let use_l = dy_l.abs() < dy_r.abs();
    let dx = a_r.x - a_l.x;
    let mut lambda = if use_l { dy_l } else { dy_r } / (dy_l - dy_r);
    if !lambda.is_finite() {
        lambda = 0.0;
    }
    
    let x = lambda * dx + if use_l { a_l.x } else { a_r.x };
    
    let a_dy = a_r.y - a_l.y;
    let b_dy = b_r.y - b_l.y;
    let use_a = a_dy.abs() < b_dy.abs();
    
    let y = lambda * if use_a { a_dy } else { b_dy } +
            if use_l {
                if use_a { a_l.y } else { b_l.y }
            } else {
                if use_a { a_r.y } else { b_r.y }
            };
    
    let z = lambda * (a_r.z - a_l.z) + if use_l { a_l.z } else { a_r.z };
    let w = lambda * (b_r.z - b_l.z) + if use_l { b_l.z } else { b_r.z };
    
    DVec4::new(x, y, z, w)
}
```

### 1.3 Shadows Predicate (Lines 81-83)

**Purpose**: Symbolic perturbation tie-breaker for degenerate cases.

**C++ Reference**:
```cpp
inline bool Shadows(double p, double q, double dir) {
  return p == q ? dir < 0 : p < q;
}
```

**Rust Port** (EXACT):
```rust
#[inline]
pub fn shadows(p: f64, q: f64, dir: f64) -> bool {
    if p == q {
        dir < 0.0
    } else {
        p < q
    }
}
```

**CRITICAL**: The `p == q` comparison MUST use exact floating-point equality, NOT epsilon comparison.

### 1.4 withSign Helper (Lines 37-42)

**C++ Reference**:
```cpp
inline double withSign(bool pos, double v) { return pos ? v : -v; }
```

**Rust Port** (EXACT):
```rust
#[inline]
pub fn with_sign(pos: bool, v: f64) -> f64 {
    if pos { v } else { -v }
}
```

---

## 2. Parallel Primitives (parallel.h)

### 2.1 Threshold Constant

```cpp
constexpr size_t kSeqThreshold = 1e4;
```

```rust
pub const K_SEQ_THRESHOLD: usize = 10_000;
```

### 2.2 Auto Policy Selection

**C++ Reference**:
```cpp
inline constexpr ExecutionPolicy autoPolicy(size_t size,
                                            size_t threshold = kSeqThreshold) {
  if (size <= threshold) {
    return ExecutionPolicy::Seq;
  }
  return ExecutionPolicy::Par;
}
```

**Rust Port**:
```rust
#[inline]
pub fn auto_policy(size: usize, threshold: usize) -> ExecutionPolicy {
    if size <= threshold {
        ExecutionPolicy::Seq
    } else {
        ExecutionPolicy::Par
    }
}

#[inline]
pub fn auto_policy_default(size: usize) -> ExecutionPolicy {
    auto_policy(size, K_SEQ_THRESHOLD)
}
```

### 2.3 Merge Sort (Lines 74-110)

This is used for stable parallel sorting. Must match exactly.

**C++ Reference**:
```cpp
template <typename SrcIter, typename DestIter, typename Comp>
void mergeRec(SrcIter src, DestIter dest, size_t p1, size_t r1, size_t p2,
              size_t r2, size_t p3, Comp comp) {
  size_t length1 = r1 - p1;
  size_t length2 = r2 - p2;
  if (length1 < length2) {
    std::swap(p1, p2);
    std::swap(r1, r2);
    std::swap(length1, length2);
  }
  if (length1 == 0) return;
  if (length1 + length2 <= kSeqThreshold) {
    std::merge(src + p1, src + r1, src + p2, src + r2, dest + p3, comp);
  } else {
    size_t q1 = p1 + length1 / 2;
    size_t q2 =
        std::distance(src, std::lower_bound(src + p2, src + r2, src[q1], comp));
    size_t q3 = p3 + (q1 - p1) + (q2 - p2);
    dest[q3] = src[q1];
    tbb::parallel_invoke(
        [=] { mergeRec(src, dest, p1, q1, p2, q2, p3, comp); },
        [=] { mergeRec(src, dest, q1 + 1, r1, q2, r2, q3 + 1, comp); });
  }
}
```

---

## 3. Trigonometric Functions (common.h)

### 3.1 sind - Sine in Degrees

**Purpose**: Exact sine for multiples of 90°

**C++ Reference**:
```cpp
inline double sind(double x) {
  if (!la::isfinite(x)) return sin(x);
  if (x < 0.0) return -sind(-x);
  int quo;
  x = remquo(fabs(x), 90.0, &quo);
  switch (quo % 4) {
    case 0:
      return sin(radians(x));
    case 1:
      return cos(radians(x));
    case 2:
      return -sin(radians(x));
    case 3:
      return -cos(radians(x));
  }
  return 0.0;
}
```

**Rust Port** (requires `remquo` implementation):
```rust
pub fn sind(x: f64) -> f64 {
    if !x.is_finite() {
        return x.sin();
    }
    if x < 0.0 {
        return -sind(-x);
    }
    let (rem, quo) = remquo(x.abs(), 90.0);
    match quo % 4 {
        0 => radians(rem).sin(),
        1 => radians(rem).cos(),
        2 => -radians(rem).sin(),
        3 => -radians(rem).cos(),
        _ => 0.0,
    }
}

pub fn cosd(x: f64) -> f64 {
    sind(x + 90.0)
}

// Must implement remquo to match C behavior
fn remquo(x: f64, y: f64) -> (f64, i32) {
    let quo = (x / y).round() as i32;
    let rem = x - (quo as f64) * y;
    (rem, quo)
}
```

---

## 4. DisjointSets (disjoint_sets.h)

### 4.1 Find with Path Compression

**C++ Reference**:
```cpp
uint32_t find(uint32_t id) const {
  while (id != parent(id)) {
    uint64_t value = mData[id];
    uint32_t new_parent = parent((uint32_t)value);
    uint64_t new_value = (value & 0xFFFFFFFF00000000ULL) | new_parent;
    /* Try to update parent (may fail, that's ok) */
    if (value != new_value) mData[id].compare_exchange_weak(value, new_value);
    id = new_parent;
  }
  return id;
}
```

**Rust Port** (atomic operations must match):
```rust
pub fn find(&self, mut id: u32) -> u32 {
    while id != self.parent(id) {
        let value = self.data[id as usize].load(Ordering::Relaxed);
        let new_parent = self.parent(value as u32);
        let new_value = (value & 0xFFFFFFFF00000000) | (new_parent as u64);
        if value != new_value {
            let _ = self.data[id as usize].compare_exchange_weak(
                value,
                new_value,
                Ordering::Relaxed,
                Ordering::Relaxed,
            );
        }
        id = new_parent;
    }
    id
}
```

### 4.2 Unite with Rank

**C++ Reference**:
```cpp
uint32_t unite(uint32_t id1, uint32_t id2) {
  for (;;) {
    id1 = find(id1);
    id2 = find(id2);

    if (id1 == id2) return id1;

    uint32_t r1 = rank(id1), r2 = rank(id2);

    if (r1 > r2 || (r1 == r2 && id1 < id2)) {
      std::swap(r1, r2);
      std::swap(id1, id2);
    }

    uint64_t oldEntry = ((uint64_t)r1 << 32) | id1;
    uint64_t newEntry = ((uint64_t)r1 << 32) | id2;

    if (!mData[id1].compare_exchange_strong(oldEntry, newEntry)) continue;

    if (r1 == r2) {
      oldEntry = ((uint64_t)r2 << 32) | id2;
      newEntry = ((uint64_t)(r2 + 1) << 32) | id2;
      /* Try to update the rank (may fail, retry if rank = 0) */
      if (!mData[id2].compare_exchange_strong(oldEntry, newEntry) && r2 == 0)
        continue;
    }

    break;
  }
  return id2;
}
```

---

## 5. Collider BVH (collider.h)

### 5.1 Node Type Helpers

```cpp
constexpr inline bool IsLeaf(int node) { return node % 2 == 0; }
constexpr inline bool IsInternal(int node) { return node % 2 == 1; }
constexpr inline int Node2Internal(int node) { return (node - 1) / 2; }
constexpr inline int Internal2Node(int internal) { return internal * 2 + 1; }
constexpr inline int Node2Leaf(int node) { return node / 2; }
constexpr inline int Leaf2Node(int leaf) { return leaf * 2; }
```

```rust
#[inline]
pub const fn is_leaf(node: i32) -> bool { node % 2 == 0 }
#[inline]
pub const fn is_internal(node: i32) -> bool { node % 2 == 1 }
#[inline]
pub const fn node_to_internal(node: i32) -> i32 { (node - 1) / 2 }
#[inline]
pub const fn internal_to_node(internal: i32) -> i32 { internal * 2 + 1 }
#[inline]
pub const fn node_to_leaf(node: i32) -> i32 { node / 2 }
#[inline]
pub const fn leaf_to_node(leaf: i32) -> i32 { leaf * 2 }
```

### 5.2 Prefix Length (Morton Code)

```cpp
int PrefixLength(uint32_t a, uint32_t b) const {
#ifdef _MSC_VER
  return clz(a ^ b);
#else
  return __builtin_clz(a ^ b);
#endif
}
```

```rust
#[inline]
fn prefix_length(a: u32, b: u32) -> i32 {
    (a ^ b).leading_zeros() as i32
}
```

---

## 6. Box Operations (common.h)

### 6.1 DoesOverlap

**C++ Reference**:
```cpp
constexpr bool DoesOverlap(const Box& box) const {
  return min.x <= box.max.x && min.y <= box.max.y && min.z <= box.max.z &&
         max.x >= box.min.x && max.y >= box.min.y && max.z >= box.min.z;
}
```

**Rust Port**:
```rust
#[inline]
pub fn does_overlap(&self, other: &Box) -> bool {
    self.min.x <= other.max.x && self.min.y <= other.max.y && self.min.z <= other.max.z &&
    self.max.x >= other.min.x && self.max.y >= other.min.y && self.max.z >= other.min.z
}
```

---

## 7. Halfedge Operations (shared.h)

### 7.1 NextHalfedge

```cpp
inline int NextHalfedge(int current) {
  ++current;
  if (current % 3 == 0) current -= 3;
  return current;
}
```

```rust
#[inline]
pub fn next_halfedge(current: i32) -> i32 {
    let next = current + 1;
    if next % 3 == 0 { next - 3 } else { next }
}
```

### 7.2 GetAxisAlignedProjection

```cpp
inline mat2x3 GetAxisAlignedProjection(vec3 normal) {
  vec3 absNormal = la::abs(normal);
  double xyzMax;
  mat3x2 projection;
  if (absNormal.z > absNormal.x && absNormal.z > absNormal.y) {
    projection = mat3x2({1.0, 0.0, 0.0},  //
                        {0.0, 1.0, 0.0});
    xyzMax = normal.z;
  } else if (absNormal.y > absNormal.x) {
    projection = mat3x2({0.0, 0.0, 1.0},  //
                        {1.0, 0.0, 0.0});
    xyzMax = normal.y;
  } else {
    projection = mat3x2({0.0, 1.0, 0.0},  //
                        {0.0, 0.0, 1.0});
    xyzMax = normal.x;
  }
  if (xyzMax < 0) projection[0] *= -1.0;
  return la::transpose(projection);
}
```

---

## Verification Checklist

For each function ported:

- [ ] Line-by-line comparison with C++ reference
- [ ] Same variable names (translated to Rust style)
- [ ] Same order of operations
- [ ] Same branching logic
- [ ] Same edge case handling
- [ ] Unit tests match C++ test expectations
- [ ] Golden test comparison passes

---

*Document Version: 1.0*
*Reference: submodules/manifold/src/ @ HEAD*
