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

use crate::node::*;
use manifold_math::Box;
use std::sync::atomic::{AtomicI32, Ordering};

pub struct BuildInternalBoxes<'a> {
    pub node_bbox: &'a mut [Box],
    pub counter: &'a [AtomicI32],
    pub node_parent: &'a [i32],
    pub internal_children: &'a [(i32, i32)],
}

impl<'a> BuildInternalBoxes<'a> {
    pub fn execute(&mut self, leaf: i32) {
        let mut node = leaf_to_node(leaf);
        loop {
            node = self.node_parent[node as usize];
            if node == -1 {
                break;
            }
            let internal = node_to_internal(node);
            if self.counter[internal as usize].fetch_add(1, Ordering::SeqCst) == 0 {
                return;
            }
            let (child1, child2) = self.internal_children[internal as usize];
            self.node_bbox[node as usize] =
                self.node_bbox[child1 as usize].union_box(&self.node_bbox[child2 as usize]);
            if node == K_ROOT {
                break;
            }
        }
    }
}
