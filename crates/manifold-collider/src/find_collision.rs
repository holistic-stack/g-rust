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

pub trait Query: Sync + Send {
    fn get_bbox(&self, i: i32) -> Box;
}

pub trait Recorder: Sync + Send {
    fn record(&self, query_idx: i32, leaf_idx: i32);
}

pub struct FindCollision<'a, Q, R, const SELF_COLLISION: bool>
where
    Q: Query,
    R: Recorder,
{
    pub query: &'a Q,
    pub node_bbox: &'a [Box],
    pub internal_children: &'a [(i32, i32)],
    pub recorder: &'a R,
}

impl<'a, Q, R, const SELF_COLLISION: bool> FindCollision<'a, Q, R, SELF_COLLISION>
where
    Q: Query,
    R: Recorder,
{
    pub fn record_collision(&self, node: i32, query_idx: i32) -> bool {
        let overlaps = self.node_bbox[node as usize].does_overlap(&self.query.get_bbox(query_idx));
        if overlaps && is_leaf(node) {
            let leaf_idx = node_to_leaf(node);
            if !SELF_COLLISION || leaf_idx != query_idx {
                self.recorder.record(query_idx, leaf_idx);
            }
        }
        overlaps && is_internal(node)
    }

    pub fn execute(&self, query_idx: i32) {
        let mut stack = [0i32; 64];
        let mut top = -1;
        let mut node = K_ROOT;

        loop {
            let internal = node_to_internal(node);
            let (child1, child2) = self.internal_children[internal as usize];

            let traverse1 = self.record_collision(child1, query_idx);
            let traverse2 = self.record_collision(child2, query_idx);

            if !traverse1 && !traverse2 {
                if top < 0 {
                    break;
                }
                node = stack[top as usize];
                top -= 1;
            } else {
                node = if traverse1 { child1 } else { child2 };
                if traverse1 && traverse2 {
                    top += 1;
                    stack[top as usize] = child2;
                }
            }
        }
    }
}
