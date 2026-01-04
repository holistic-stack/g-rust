// Copyright 2021 The Manifold Authors.
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

pub mod build_boxes;
pub mod find_collision;
pub mod node;
pub mod radix_tree;

use build_boxes::BuildInternalBoxes;
use find_collision::{FindCollision, Query, Recorder};
use glam::{DMat4, DVec3};
use manifold_math::Box;
use manifold_parallel::{auto_policy, for_each, for_each_n, ExecutionPolicy};
use node::*;
use radix_tree::CreateRadixTree;
use std::sync::atomic::AtomicI32;

pub struct Collider {
    node_bbox: Vec<Box>,
    node_parent: Vec<i32>,
    internal_children: Vec<(i32, i32)>,
}

impl Collider {
    pub fn new(leaf_bb: &[Box], leaf_morton: &[u32]) -> Self {
        assert_eq!(leaf_bb.len(), leaf_morton.len());
        let num_leaves = leaf_bb.len() as i32;
        if num_leaves == 0 {
            return Self {
                node_bbox: Vec::new(),
                node_parent: Vec::new(),
                internal_children: Vec::new(),
            };
        }
        let num_nodes = 2 * num_leaves - 1;
        let node_bbox = vec![Box::default(); num_nodes as usize];
        let mut node_parent = vec![-1i32; num_nodes as usize];
        let mut internal_children = vec![(-1i32, -1i32); (num_leaves - 1) as usize];

        let num_internal = num_leaves - 1;
        if num_internal > 0 {
            let mut radix_tree = CreateRadixTree {
                node_parent: &mut node_parent,
                internal_children: &mut internal_children,
                leaf_morton,
            };

            for i in 0..num_internal {
                radix_tree.execute(i);
            }
        }

        let mut collider = Self {
            node_bbox,
            node_parent,
            internal_children,
        };
        collider.update_boxes(leaf_bb);
        collider
    }

    pub fn transform(&mut self, transform: DMat4) -> bool {
        let mut axis_aligned = true;
        let cols = transform.to_cols_array_2d();
        for i in 0..3 {
            let mut count = 0;
            for j in 0..3 {
                if cols[j][i] == 0.0 {
                    count += 1;
                }
            }
            if count != 2 {
                axis_aligned = false;
            }
        }

        if axis_aligned {
            let policy = auto_policy(self.node_bbox.len(), 100_000);
            for_each(policy, &mut self.node_bbox, |box_| {
                *box_ = box_.transform(transform);
            });
        }
        axis_aligned
    }

    pub fn update_boxes(&mut self, leaf_bb: &[Box]) {
        assert_eq!(leaf_bb.len() as i32, self.num_leaves());
        // Copy in leaf node boxes
        for (i, bb) in leaf_bb.iter().enumerate() {
            self.node_bbox[leaf_to_node(i as i32) as usize] = *bb;
        }

        let num_internal = self.num_internal();
        if num_internal == 0 {
            return;
        }

        let counter: Vec<AtomicI32> = (0..num_internal).map(|_| AtomicI32::new(0)).collect();
        let num_leaves = leaf_bb.len() as i32;

        // Build internal boxes
        let node_bbox_ptr = self.node_bbox.as_mut_ptr() as usize;
        let node_bbox_len = self.node_bbox.len();

        for_each_n(
            auto_policy(num_leaves as usize, 1000),
            0,
            num_leaves as usize,
            |leaf| {
                let node_bbox = unsafe {
                    std::slice::from_raw_parts_mut(node_bbox_ptr as *mut Box, node_bbox_len)
                };
                let mut build_boxes = BuildInternalBoxes {
                    node_bbox,
                    counter: &counter,
                    node_parent: &self.node_parent,
                    internal_children: &self.internal_children,
                };
                build_boxes.execute(leaf as i32);
            },
        );
    }

    pub fn collisions<const SELF_COLLISION: bool, Q, R>(
        &self,
        query: &Q,
        n: i32,
        recorder: &R,
        parallel: bool,
    ) where
        Q: Query,
        R: Recorder,
    {
        if self.internal_children.is_empty() {
            return;
        }
        let policy = if parallel {
            auto_policy(n as usize, K_SEQUENTIAL_THRESHOLD as usize)
        } else {
            ExecutionPolicy::Seq
        };

        for_each_n(policy, 0, n as usize, |i| {
            let find_collision = FindCollision::<Q, R, SELF_COLLISION> {
                query,
                node_bbox: &self.node_bbox,
                internal_children: &self.internal_children,
                recorder,
            };
            find_collision.execute(i as i32);
        });
    }

    pub fn num_internal(&self) -> i32 {
        self.internal_children.len() as i32
    }

    pub fn num_leaves(&self) -> i32 {
        if self.internal_children.is_empty() {
            0
        } else {
            self.num_internal() + 1
        }
    }

    pub fn spread_bits_3(mut v: u32) -> u32 {
        v = 0xFF0000FF & (v.wrapping_mul(0x00010001));
        v = 0x0F00F00F & (v.wrapping_mul(0x00000101));
        v = 0xC30C30C3 & (v.wrapping_mul(0x00000011));
        v = 0x49249249 & (v.wrapping_mul(0x00000005));
        v
    }

    pub fn morton_code(position: DVec3, bbox: Box) -> u32 {
        let xyz = (position - bbox.min) / (bbox.max - bbox.min);
        let x = f64::min(1023.0, f64::max(0.0, 1024.0 * xyz.x)) as u32;
        let y = f64::min(1023.0, f64::max(0.0, 1024.0 * xyz.y)) as u32;
        let z = f64::min(1023.0, f64::max(0.0, 1024.0 * xyz.z)) as u32;
        Self::spread_bits_3(x) * 4 + Self::spread_bits_3(y) * 2 + Self::spread_bits_3(z)
    }
}
