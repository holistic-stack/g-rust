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

use glam::{DVec4, IVec3, IVec4};

pub struct Partition {
    pub idx: IVec4,
    pub sorted_divisions: IVec4,
    pub vert_bary: Vec<DVec4>,
    pub tri_vert: Vec<IVec3>,
}

impl Partition {
    pub fn new() -> Self {
        Self {
            idx: IVec4::ZERO,
            sorted_divisions: IVec4::ZERO,
            vert_bary: Vec::new(),
            tri_vert: Vec::new(),
        }
    }

    pub fn interior_offset(&self) -> i32 {
        self.sorted_divisions[0]
            + self.sorted_divisions[1]
            + self.sorted_divisions[2]
            + self.sorted_divisions[3]
    }

    pub fn num_interior(&self) -> usize {
        self.vert_bary.len() - self.interior_offset() as usize
    }

    pub fn get_partition(divisions: IVec4) -> Self {
        if divisions[0] == 0 {
            return Self::new();
        }

        let mut sorted_div = divisions;
        let mut tri_idx = IVec4::new(0, 1, 2, 3);
        if divisions[3] == 0 {
            if sorted_div[2] > sorted_div[1] {
                sorted_div = IVec4::new(sorted_div[0], sorted_div[2], sorted_div[1], 0);
                tri_idx = IVec4::new(tri_idx[0], tri_idx[2], tri_idx[1], 3);
            }
            if sorted_div[1] > sorted_div[0] {
                sorted_div = IVec4::new(sorted_div[1], sorted_div[0], sorted_div[2], 0);
                tri_idx = IVec4::new(tri_idx[1], tri_idx[0], tri_idx[2], 3);
                if sorted_div[2] > sorted_div[1] {
                    sorted_div = IVec4::new(sorted_div[0], sorted_div[2], sorted_div[1], 0);
                    tri_idx = IVec4::new(tri_idx[0], tri_idx[2], tri_idx[1], 3);
                }
            }
        } else {
            let mut min_idx = 0;
            let mut min = divisions[min_idx];
            let mut next = divisions[1];
            for i in 1..4 {
                let n = divisions[(i + 1) % 4];
                if divisions[i] < min || (divisions[i] == min && n < next) {
                    min_idx = i;
                    min = divisions[i];
                    next = n;
                }
            }
            let tmp = sorted_div;
            for i in 0..4 {
                tri_idx[i] = ((i + min_idx) % 4) as i32;
                sorted_div[i] = tmp[tri_idx[i] as usize];
            }
        }

        let mut partition = Self::get_cached_partition(sorted_div);
        partition.idx = tri_idx;
        partition
    }

    fn get_cached_partition(n: IVec4) -> Self {
        let mut partition = Self::new();
        partition.sorted_divisions = n;
        if n[3] > 0 {
            partition.vert_bary.push(DVec4::new(1.0, 0.0, 0.0, 0.0));
            partition.vert_bary.push(DVec4::new(0.0, 1.0, 0.0, 0.0));
            partition.vert_bary.push(DVec4::new(0.0, 0.0, 1.0, 0.0));
            partition.vert_bary.push(DVec4::new(0.0, 0.0, 0.0, 1.0));
        } else {
            partition.vert_bary.push(DVec4::new(1.0, 0.0, 0.0, 0.0));
            partition.vert_bary.push(DVec4::new(0.0, 1.0, 0.0, 0.0));
            partition.vert_bary.push(DVec4::new(0.0, 0.0, 1.0, 0.0));
            if n[0] == 1 && n[1] == 1 && n[2] == 1 {
                partition.tri_vert.push(IVec3::new(0, 1, 2));
            }
        }
        partition
    }
}
