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

pub struct CreateRadixTree<'a> {
    pub node_parent: &'a mut [i32],
    pub internal_children: &'a mut [(i32, i32)],
    pub leaf_morton: &'a [u32],
}

impl<'a> CreateRadixTree<'a> {
    pub fn prefix_length_u32(&self, a: u32, b: u32) -> i32 {
        (a ^ b).leading_zeros() as i32
    }

    pub fn prefix_length(&self, i: i32, j: i32) -> i32 {
        if j < 0 || j >= self.leaf_morton.len() as i32 {
            -1
        } else {
            if self.leaf_morton[i as usize] == self.leaf_morton[j as usize] {
                32 + self.prefix_length_u32(i as u32, j as u32)
            } else {
                self.prefix_length_u32(self.leaf_morton[i as usize], self.leaf_morton[j as usize])
            }
        }
    }

    pub fn range_end(&self, i: i32) -> i32 {
        let mut dir = self.prefix_length(i, i + 1) - self.prefix_length(i, i - 1);
        dir = if dir > 0 {
            1
        } else if dir < 0 {
            -1
        } else {
            0
        };

        let common_prefix = self.prefix_length(i, i - dir);
        let mut max_length = K_INITIAL_LENGTH;
        while self.prefix_length(i, i + dir * max_length) > common_prefix {
            max_length *= K_LENGTH_MULTIPLE;
        }

        let mut length = 0;
        let mut step = max_length / 2;
        while step > 0 {
            if self.prefix_length(i, i + dir * (length + step)) > common_prefix {
                length += step;
            }
            step /= 2;
        }
        i + dir * length
    }

    pub fn find_split(&self, first: i32, last: i32) -> i32 {
        let common_prefix = self.prefix_length(first, last);
        let mut split = first;
        let mut step = last - first;
        loop {
            step = (step + 1) >> 1;
            let new_split = split + step;
            if new_split < last {
                let split_prefix = self.prefix_length(first, new_split);
                if split_prefix > common_prefix {
                    split = new_split;
                }
            }
            if step <= 1 {
                break;
            }
        }
        split
    }

    pub fn execute(&mut self, internal: i32) {
        let first = internal;
        let mut last = self.range_end(first);
        let mut first_actual = first;
        if first_actual > last {
            std::mem::swap(&mut first_actual, &mut last);
        }
        let split = self.find_split(first_actual, last);
        let child1 = if split == first_actual {
            leaf_to_node(split)
        } else {
            internal_to_node(split)
        };
        let split_plus_1 = split + 1;
        let child2 = if split_plus_1 == last {
            leaf_to_node(split_plus_1)
        } else {
            internal_to_node(split_plus_1)
        };

        self.internal_children[internal as usize] = (child1, child2);
        let node = internal_to_node(internal);
        self.node_parent[child1 as usize] = node;
        self.node_parent[child2 as usize] = node;
    }
}
