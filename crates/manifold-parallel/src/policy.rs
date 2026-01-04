// Copyright 2022 The Manifold Authors.
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

pub const K_SEQ_THRESHOLD: usize = 10_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionPolicy {
    Par,
    Seq,
}

pub fn auto_policy(size: usize, threshold: usize) -> ExecutionPolicy {
    if size <= threshold {
        ExecutionPolicy::Seq
    } else {
        ExecutionPolicy::Par
    }
}
