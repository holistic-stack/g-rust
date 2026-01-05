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

use glam::{DVec2, DVec3, DVec4, IVec3};
use manifold_collider::Collider;
use manifold_math::{morton_code, Box as GeoBox, K_NO_CODE};
use manifold_types::{next_halfedge, Halfedge, OpType, TriRef};

pub mod boolean_result;
pub mod kernel02;
pub mod kernel11;
pub mod kernel12;
pub mod shadow01;
pub mod subdivision;
pub mod winding03;

pub use boolean_result::create_boolean_result;
pub use subdivision::Partition;

/// Intersections between manifolds P and Q.
#[derive(Debug, Clone)]
pub struct Intersections {
    pub p1q2: Vec<[i32; 2]>,
    pub x12: Vec<i32>,
    pub v12: Vec<DVec3>,
}

impl Intersections {
    pub fn new() -> Self {
        Self {
            p1q2: Vec::new(),
            x12: Vec::new(),
            v12: Vec::new(),
        }
    }
}

pub struct Boolean3<'a> {
    pub in_p: &'a ManifoldImpl,
    pub in_q: &'a ManifoldImpl,
    pub expand_p: bool,
    pub xv12_: Intersections,
    pub xv21_: Intersections,
    pub w03_: Vec<i32>,
    pub w30_: Vec<i32>,
    pub valid: bool,
}

impl<'a> Boolean3<'a> {
    pub fn new(in_p: &'a ManifoldImpl, in_q: &'a ManifoldImpl, op: OpType) -> Self {
        let expand_p = matches!(op, OpType::Add);
        let xv12 = Intersections::new();
        let xv21 = Intersections::new();
        let w03 = winding03::winding03_true_true(in_p, in_q, &xv12, expand_p);
        let w30 = winding03::winding03_false_true(in_p, in_q, &xv21, expand_p);

        Self {
            in_p,
            in_q,
            expand_p,
            xv12_: xv12,
            xv21_: xv21,
            w03_: w03,
            w30_: w30,
            valid: true,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Relation {
    pub original_id: i32,
    pub transform: glam::DMat4,
    pub back_side: bool,
}

#[derive(Debug, Clone, Default)]
pub struct MeshRelationD {
    pub original_id: i32,
    pub mesh_id_transform: std::collections::HashMap<i32, Relation>,
    pub tri_ref: Vec<TriRef>,
}

#[derive(Debug, Clone, Default)]
pub struct ManifoldImpl {
    pub vert_pos: Vec<DVec3>,
    pub halfedge: Vec<Halfedge>,
    pub face_normal: Vec<DVec3>,
    pub vert_normal: Vec<DVec3>,
    pub properties: Vec<f64>,
    pub num_prop: usize,
    pub mesh_relation: MeshRelationD,
    pub collider: Option<Collider>,
    pub bbox: GeoBox,
    pub epsilon: f64,
    pub tolerance: f64,
}

#[derive(Debug, Clone, Copy)]
pub enum Shape {
    Tetrahedron,
    Cube,
    Octahedron,
}

impl ManifoldImpl {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_shape(shape: Shape, m: glam::DMat4) -> Self {
        let mut vert_pos = Vec::new();
        let mut tri_verts = Vec::new();
        match shape {
            Shape::Tetrahedron => {
                vert_pos = vec![
                    DVec3::new(-1.0, -1.0, 1.0),
                    DVec3::new(-1.0, 1.0, -1.0),
                    DVec3::new(1.0, -1.0, -1.0),
                    DVec3::new(1.0, 1.0, 1.0),
                ];
                tri_verts = vec![
                    IVec3::new(2, 0, 1),
                    IVec3::new(0, 3, 1),
                    IVec3::new(2, 3, 0),
                    IVec3::new(3, 2, 1),
                ];
            }
            Shape::Cube => {
                vert_pos = vec![
                    DVec3::new(0.0, 0.0, 0.0),
                    DVec3::new(0.0, 0.0, 1.0),
                    DVec3::new(0.0, 1.0, 0.0),
                    DVec3::new(0.0, 1.0, 1.0),
                    DVec3::new(1.0, 0.0, 0.0),
                    DVec3::new(1.0, 0.0, 1.0),
                    DVec3::new(1.0, 1.0, 0.0),
                    DVec3::new(1.0, 1.0, 1.0),
                ];
                tri_verts = vec![
                    IVec3::new(1, 0, 4),
                    IVec3::new(2, 4, 0),
                    IVec3::new(1, 3, 0),
                    IVec3::new(3, 1, 5),
                    IVec3::new(3, 2, 0),
                    IVec3::new(3, 7, 2),
                    IVec3::new(5, 4, 6),
                    IVec3::new(5, 1, 4),
                    IVec3::new(6, 4, 2),
                    IVec3::new(7, 6, 2),
                    IVec3::new(7, 3, 5),
                    IVec3::new(7, 5, 6),
                ];
            }
            Shape::Octahedron => {
                vert_pos = vec![
                    DVec3::new(1.0, 0.0, 0.0),
                    DVec3::new(-1.0, 0.0, 0.0),
                    DVec3::new(0.0, 1.0, 0.0),
                    DVec3::new(0.0, -1.0, 0.0),
                    DVec3::new(0.0, 0.0, 1.0),
                    DVec3::new(0.0, 0.0, -1.0),
                ];
                tri_verts = vec![
                    IVec3::new(0, 2, 4),
                    IVec3::new(1, 5, 3),
                    IVec3::new(2, 1, 4),
                    IVec3::new(3, 5, 0),
                    IVec3::new(1, 3, 4),
                    IVec3::new(0, 5, 2),
                    IVec3::new(3, 0, 4),
                    IVec3::new(2, 5, 1),
                ];
            }
        }

        let mut transformed_vert_pos = Vec::with_capacity(vert_pos.len());
        for v in vert_pos {
            let v_h = m * glam::DVec4::new(v.x, v.y, v.z, 1.0);
            transformed_vert_pos.push(DVec3::new(v_h.x, v_h.y, v_h.z));
        }

        let mut impl_ = Self {
            vert_pos: transformed_vert_pos,
            ..Self::default()
        };

        impl_.create_halfedges(&tri_verts);
        impl_.finish();
        impl_
    }

    pub fn subdivide<F>(&mut self, _edge_divisions: F)
    where
        F: Fn(DVec3, DVec4, DVec4) -> i32,
    {
    }

    pub fn create_halfedges(&mut self, tri_verts: &[IVec3]) {
        let num_tri = tri_verts.len();
        let num_halfedge = 3 * num_tri;
        self.halfedge = Vec::with_capacity(num_halfedge);

        for tri in tri_verts {
            self.halfedge.push(Halfedge {
                start_vert: tri.x,
                end_vert: tri.y,
                paired_halfedge: -1,
                prop_vert: tri.x,
            });
            self.halfedge.push(Halfedge {
                start_vert: tri.y,
                end_vert: tri.z,
                paired_halfedge: -1,
                prop_vert: tri.y,
            });
            self.halfedge.push(Halfedge {
                start_vert: tri.z,
                end_vert: tri.x,
                paired_halfedge: -1,
                prop_vert: tri.z,
            });
        }

        let mut edge_map: std::collections::HashMap<(i32, i32), usize> = Default::default();
        for i in 0..num_halfedge {
            let h = &self.halfedge[i];
            let key = if h.start_vert < h.end_vert {
                (h.start_vert, h.end_vert)
            } else {
                (h.end_vert, h.start_vert)
            };

            if let Some(&other) = edge_map.get(&key) {
                self.halfedge[i].paired_halfedge = other as i32;
                self.halfedge[other].paired_halfedge = i as i32;
                edge_map.remove(&key);
            } else {
                edge_map.insert(key, i);
            }
        }
    }

    pub fn get_face_box_morton(&self, face_box: &mut Vec<GeoBox>, face_morton: &mut Vec<u32>) {
        let num_tri = self.halfedge.len() / 3;
        face_box.resize(num_tri, GeoBox::default());
        face_morton.resize(num_tri, 0);

        for face in 0..num_tri {
            if self.halfedge[3 * face].paired_halfedge < 0 {
                face_morton[face] = K_NO_CODE;
                continue;
            }

            let mut center = DVec3::ZERO;
            let mut bbox = GeoBox::default();
            for i in 0..3 {
                let pos = self.vert_pos[self.halfedge[3 * face + i].start_vert as usize];
                center += pos;
                bbox.union_point(pos);
            }
            center /= 3.0;
            face_box[face] = bbox;
            face_morton[face] = morton_code(center, self.bbox);
        }
    }

    pub fn calculate_normals(&mut self) {
        let num_vert = self.vert_pos.len();
        let num_tri = self.halfedge.len() / 3;
        self.face_normal = vec![DVec3::new(0.0, 0.0, 1.0); num_tri];
        self.vert_normal = vec![DVec3::ZERO; num_vert];

        let mut vert_halfedge_map = vec![i32::MAX; num_vert];

        for face in 0..num_tri {
            let v0 = self.halfedge[3 * face].start_vert as usize;
            let v1 = self.halfedge[3 * face + 1].start_vert as usize;
            let v2 = self.halfedge[3 * face + 2].start_vert as usize;

            let edge0 = self.vert_pos[v1] - self.vert_pos[v0];
            let edge1 = self.vert_pos[v2] - self.vert_pos[v1];
            let tri_normal = edge0.cross(edge1).normalize_or_zero();
            self.face_normal[face] = if tri_normal.length_squared() > 0.0 {
                tri_normal
            } else {
                DVec3::new(0.0, 0.0, 1.0)
            };

            for i in 0..3 {
                let v = self.halfedge[3 * face + i].start_vert as usize;
                let edge_idx = (3 * face + i) as i32;
                if edge_idx < vert_halfedge_map[v] {
                    vert_halfedge_map[v] = edge_idx;
                }
            }
        }

        for v in 0..num_vert {
            let mut normal = DVec3::ZERO;
            let first_edge = vert_halfedge_map[v];
            if first_edge == i32::MAX {
                continue;
            }

            let mut curr = first_edge;
            loop {
                let face_idx = (curr / 3) as usize;
                normal += self.face_normal[face_idx];

                let paired = self.halfedge[curr as usize].paired_halfedge;
                if paired < 0 {
                    break;
                }
                curr = next_halfedge(paired);
                if curr == first_edge {
                    break;
                }
            }
            self.vert_normal[v] = normal.normalize_or_zero();
        }
    }

    pub fn finish(&mut self) {
        if self.halfedge.is_empty() {
            return;
        }

        self.calculate_normals();

        let mut bbox = GeoBox::default();
        for &v in &self.vert_pos {
            bbox.union_point(v);
        }
        self.bbox = bbox;

        let mut face_box = Vec::new();
        let mut face_morton = Vec::new();
        self.get_face_box_morton(&mut face_box, &mut face_morton);

        self.collider = Some(Collider::new(&face_box, &face_morton));
    }

    pub fn initialize_original(&mut self) {}
    pub fn mark_coplanar(&mut self) {}
}

pub use shadow01::shadow01_functions as shadow01_f;
