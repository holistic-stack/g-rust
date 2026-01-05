// Copyright 2024 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Thread-safe Union-Find data structure
// Port of C++ DisjointSets from src/disjoint_sets.h
// Original from: https://github.com/wjakob/dset

use std::sync::atomic::{AtomicU64, Ordering};

/// Thread-safe Union-Find data structure.
///
/// Uses atomic operations to store both parent and rank in a single u64:
/// - Upper 32 bits: rank
/// - Lower 32 bits: parent
///
/// C++ Reference: disjoint_sets.h:29-121
pub struct DisjointSets {
    /// Stores parent in lower 32 bits, rank in upper 32 bits
    data: Vec<AtomicU64>,
}

impl DisjointSets {
    /// Create a new DisjointSets with `size` elements.
    ///
    /// Each element starts as its own set (parent points to itself).
    ///
    /// C++ Reference: disjoint_sets.h:31-33
    pub fn new(size: u32) -> Self {
        let mut data = Vec::with_capacity(size as usize);
        for i in 0..size {
            // Each element starts as its own parent, rank 0
            data.push(AtomicU64::new(i as u64));
        }
        Self { data }
    }

    /// Extract parent from packed data (lower 32 bits).
    ///
    /// C++ Reference: disjoint_sets.h:94
    #[inline]
    fn parent(&self, id: u32) -> u32 {
        (self.data[id as usize].load(Ordering::Relaxed) as u32) & 0xFFFFFFFF
    }

    /// Extract rank from packed data (upper 32 bits).
    ///
    /// C++ Reference: disjoint_sets.h:90-92
    #[inline]
    fn rank(&self, id: u32) -> u32 {
        ((self.data[id as usize].load(Ordering::Relaxed) >> 32) as u32)
    }

    /// Find with path compression (atomic).
    ///
    /// Returns the root of the set containing `id`.
    ///
    /// C++ Reference: disjoint_sets.h:35-45
    pub fn find(&self, id: u32) -> u32 {
        let mut id = id;
        while id != self.parent(id) {
            let value = self.data[id as usize].load(Ordering::Relaxed);
            let new_parent = self.parent(value as u32);
            // Pack new parent with existing rank
            let new_value = (value & 0xFFFFFFFF00000000) | (new_parent as u64);
            // Try to update parent (may fail, that's ok)
            let _ = self.data[id as usize].compare_exchange_weak(
                value,
                new_value,
                Ordering::Relaxed,
                Ordering::Relaxed,
            );
            id = new_parent;
        }
        id
    }

    /// Check if two elements are in the same set.
    ///
    /// C++ Reference: disjoint_sets.h:47-54
    pub fn same(&self, id1: u32, id2: u32) -> bool {
        loop {
            let root1 = self.find(id1);
            let root2 = self.find(id2);
            if root1 == root2 {
                return true;
            }
            if self.parent(root1) == root1 {
                return false;
            }
        }
    }

    /// Unite two sets with rank heuristic (atomic).
    ///
    /// Returns the root of the united set.
    ///
    /// C++ Reference: disjoint_sets.h:56-86
    pub fn unite(&self, id1: u32, id2: u32) -> u32 {
        loop {
            let mut id1 = self.find(id1);
            let mut id2 = self.find(id2);

            if id1 == id2 {
                return id1;
            }

            let mut r1 = self.rank(id1);
            let mut r2 = self.rank(id2);

            // Ensure r1 >= r2 and id1 < id2 for deterministic merge
            if r1 > r2 || (r1 == r2 && id1 < id2) {
                std::mem::swap(&mut r1, &mut r2);
                std::mem::swap(&mut id1, &mut id2);
            }

            // Pack rank and parent
            let old_entry = ((r1 as u64) << 32) | (id1 as u64);
            let new_entry = ((r1 as u64) << 32) | (id2 as u64);

            if self.data[id1 as usize]
                .compare_exchange(old_entry, new_entry, Ordering::SeqCst, Ordering::Relaxed)
                .is_ok()
            {
                if r1 == r2 {
                    let old_entry2 = ((r2 as u64) << 32) | (id2 as u64);
                    let new_entry2 = (((r2 + 1) as u64) << 32) | (id2 as u64);
                    if self.data[id2 as usize]
                        .compare_exchange(
                            old_entry2,
                            new_entry2,
                            Ordering::SeqCst,
                            Ordering::Relaxed,
                        )
                        .is_err()
                        && r2 == 0
                    {
                        continue;
                    }
                }
                break;
            }
        }
        id2
    }

    /// Count connected components.
    ///
    /// Fills `components` with component ID for each element.
    /// Returns the number of connected components.
    ///
    /// C++ Reference: disjoint_sets.h:96-118
    pub fn connected_components(&self, components: &mut [i32]) -> usize {
        let num_nodes = components.len();
        let mut lonely_nodes = 0i32;
        let mut to_label = std::collections::HashMap::new();

        for i in 0..num_nodes {
            // Optimize for connected component of size 1
            // No need to put them into the hashmap
            let i_parent = self.find(i as u32) as usize;
            if self.rank(i_parent as u32) == 0 {
                components[i] = (to_label.len() as i32) + lonely_nodes;
                lonely_nodes += 1;
                continue;
            }

            match to_label.get(&i_parent) {
                Some(&label) => {
                    components[i] = label as i32;
                }
                None => {
                    let label = (to_label.len() as i32) + lonely_nodes;
                    to_label.insert(i_parent, label);
                    components[i] = label;
                }
            }
        }

        to_label.len() + (lonely_nodes as usize)
    }

    /// Get the size of the DisjointSets.
    ///
    /// C++ Reference: disjoint_sets.h:88
    pub fn size(&self) -> u32 {
        self.data.len() as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_disjoint_sets_basic() {
        let ds = DisjointSets::new(5);

        // Initially each element is its own set
        for i in 0..5 {
            assert!(ds.same(i, i));
            for j in (i + 1)..5 {
                assert!(!ds.same(i, j));
            }
        }
    }

    #[test]
    fn test_disjoint_sets_unite() {
        let ds = DisjointSets::new(10);

        // Unite some elements
        ds.unite(0, 1);
        assert!(ds.same(0, 1));
        assert!(!ds.same(0, 2));

        ds.unite(2, 3);
        assert!(ds.same(2, 3));
        assert!(!ds.same(1, 2));

        // Unite two sets
        ds.unite(1, 2);
        assert!(ds.same(0, 1));
        assert!(ds.same(2, 3));
        assert!(ds.same(0, 3));

        // Element 4 is still separate
        assert!(!ds.same(0, 4));
    }

    #[test]
    fn test_disjoint_sets_path_compression() {
        let ds = DisjointSets::new(100);

        // Create a long chain
        for i in 0..99 {
            ds.unite(i, i + 1);
        }

        // All should be connected
        for i in 0..100 {
            for j in 0..100 {
                assert!(ds.same(i, j));
            }
        }

        // Find should compress paths
        let root = ds.find(0);
        assert_eq!(root, ds.find(99));
    }

    #[test]
    fn test_disjoint_sets_connected_components() {
        let ds = DisjointSets::new(7);

        // Create three separate components
        ds.unite(0, 1);
        ds.unite(1, 2);
        // Component 1: {0, 1, 2}

        ds.unite(3, 4);
        // Component 2: {3, 4}

        // Component 3: {5, 6}
        ds.unite(5, 6);

        let mut components = vec![0i32; 7];
        let num_components = ds.connected_components(&mut components);

        // Should have 3 components
        assert_eq!(num_components, 3);

        assert_eq!(components[0], components[1]);
        assert_eq!(components[1], components[2]);
        assert_eq!(components[3], components[4]);
        assert_eq!(components[5], components[6]);

        assert_ne!(components[0], components[3]);
        assert_ne!(components[0], components[5]);
        assert_ne!(components[3], components[5]);
    }

    #[test]
    fn test_disjoint_sets_size() {
        let ds = DisjointSets::new(42);
        assert_eq!(ds.size(), 42);
    }
}
