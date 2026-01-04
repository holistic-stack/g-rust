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

pub const K_INITIAL_LENGTH: i32 = 128;
pub const K_LENGTH_MULTIPLE: i32 = 4;
pub const K_SEQUENTIAL_THRESHOLD: i32 = 512;
pub const K_ROOT: i32 = 1;

#[inline]
pub fn is_leaf(node: i32) -> bool {
    node % 2 == 0
}

#[inline]
pub fn is_internal(node: i32) -> bool {
    node % 2 == 1
}

#[inline]
pub fn node_to_internal(node: i32) -> i32 {
    (node - 1) / 2
}

#[inline]
pub fn internal_to_node(internal: i32) -> i32 {
    internal * 2 + 1
}

#[inline]
pub fn node_to_leaf(node: i32) -> i32 {
    node / 2
}

#[inline]
pub fn leaf_to_node(leaf: i32) -> i32 {
    leaf * 2
}
