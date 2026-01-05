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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Halfedge {
    pub start_vert: i32,
    pub end_vert: i32,
    pub paired_halfedge: i32,
    pub prop_vert: i32,
}

impl Halfedge {
    pub fn is_forward(&self) -> bool {
        self.start_vert < self.end_vert
    }
}

pub fn next_halfedge(current: i32) -> i32 {
    let mut next = current + 1;
    if next % 3 == 0 {
        next -= 3;
    }
    next
}
