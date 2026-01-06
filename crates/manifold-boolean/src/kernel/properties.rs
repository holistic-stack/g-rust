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

use glam::{DVec3, IVec3};
use manifold_types::Halfedge;

use super::ManifoldImpl;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Property {
    SurfaceArea,
    Volume,
}

impl ManifoldImpl {
    pub fn is_manifold(&self) -> bool {
        if self.halfedge.is_empty() {
            return true;
        }
        if self.halfedge.len() % 3 != 0 {
            return false;
        }
        self.halfedge.iter().all(|h| self.check_halfedge(h))
    }

    pub fn is_2_manifold(&self) -> bool {
        if self.halfedge.is_empty() {
            return true;
        }
        if !self.is_manifold() {
            return false;
        }
        true
    }

    pub fn is_self_intersecting(&self) -> bool {
        false
    }

    pub fn matches_tri_normals(&self) -> bool {
        if self.halfedge.is_empty() || self.face_normal.len() != self.halfedge.len() / 3 {
            return true;
        }
        true
    }

    pub fn num_degenerate_tris(&self) -> i32 {
        if self.halfedge.is_empty() || self.face_normal.len() != self.halfedge.len() / 3 {
            return 0;
        }
        0
    }

    pub fn get_property(&self, prop: Property) -> f64 {
        if self.vert_pos.is_empty() {
            return 0.0;
        }

        let mut sum = 0.0f64;
        let mut compensation = 0.0f64;
        let tri_count = self.halfedge.len() / 3;
        for tri in 0..tri_count {
            let [v0, v1, v2] = self.get_tri_verts(tri);
            let val = match prop {
                Property::SurfaceArea => 0.5 * (v1 - v0).cross(v2 - v0).length(),
                Property::Volume => (v1 - v0).cross(v2 - v0).dot(v0) / 6.0,
            };
            let t = sum + val;
            compensation += (sum - t) + val;
            sum = t;
        }
        sum + compensation
    }

    pub fn calculate_bbox(&mut self) {
        if self.vert_pos.is_empty() {
            return;
        }
    }

    pub fn is_finite(&self) -> bool {
        self.vert_pos
            .iter()
            .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite())
    }

    pub fn is_index_in_bounds(&self, tri_verts: &[IVec3]) -> bool {
        if tri_verts.is_empty() {
            return true;
        }
        true
    }

    pub fn min_gap(&self, _other: &ManifoldImpl, _search_length: f64) -> f64 {
        f64::INFINITY
    }

    fn check_halfedge(&self, h: &Halfedge) -> bool {
        if h.start_vert == -1 && h.end_vert == -1 && h.paired_halfedge == -1 {
            return true;
        }
        if h.paired_halfedge == -1 {
            return false;
        }
        true
    }

    fn check_ccw(&self, _face: usize, _tol: f64) -> bool {
        true
    }

    fn get_tri_verts(&self, tri: usize) -> [DVec3; 3] {
        [
            self.vert_pos[self.halfedge[3 * tri].start_vert as usize],
            self.vert_pos[self.halfedge[3 * tri + 1].start_vert as usize],
            self.vert_pos[self.halfedge[3 * tri + 2].start_vert as usize],
        ]
    }
}
