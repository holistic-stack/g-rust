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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct TmpEdge {
    pub first: i32,
    pub second: i32,
    pub halfedge_idx: i32,
}

impl TmpEdge {
    pub fn new(start: i32, end: i32, idx: i32) -> Self {
        Self {
            first: i32::min(start, end),
            second: i32::max(start, end),
            halfedge_idx: idx,
        }
    }
}
