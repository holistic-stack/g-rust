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

use crate::kernel::ManifoldImpl;
use glam::{DMat3, DQuat, DVec3, DVec4};
use manifold_math::{K_PI, K_TWO_PI};
use manifold_types::{Barycentric, TriRef};
use std::collections::HashMap;

/// Smoothness information for an edge
#[derive(Debug, Clone, Copy)]
pub struct Smoothness {
    pub halfedge: i32,
    pub smoothness: f64,
}

fn safe_normalize(v: DVec3) -> DVec3 {
    v.normalize_or_zero()
}

fn k_precision() -> f64 {
    1e-5
}

/// Returns a normalized vector orthogonal to ref, in the plane of ref and in,
/// unless in and ref are colinear, in which case it falls back to the plane of
/// ref and altIn.
fn orthogonal_to(in_vec: DVec3, alt_in: DVec3, ref_vec: DVec3) -> DVec3 {
    let dot = in_vec.dot(ref_vec);
    let mut out = in_vec - dot * ref_vec;
    if out.dot(out) < k_precision() * in_vec.dot(in_vec) {
        let alt_dot = alt_in.dot(ref_vec);
        out = alt_in - alt_dot * ref_vec;
    }
    safe_normalize(out)
}

/// Wrap radians to [-π, π]
fn wrap(radians: f64) -> f64 {
    if radians < -K_PI {
        radians + K_TWO_PI
    } else if radians > K_PI {
        radians - K_TWO_PI
    } else {
        radians
    }
}

/// Get the angle between two unit-vectors.
fn angle_between(a: DVec3, b: DVec3) -> f64 {
    let dot = a.dot(b);
    if dot >= 1.0 {
        0.0
    } else if dot <= -1.0 {
        K_PI
    } else {
        dot.acos()
    }
}

/// Calculate a tangent vector in the form of a weighted cubic Bezier taking as
/// input the desired tangent direction (length doesn't matter) and the edge
/// vector to the neighboring vertex. In a symmetric situation where the tangents
/// at each end are mirror images of each other, this will result in a circular
/// arc.
fn circular_tangent(tangent: DVec3, edge_vec: DVec3) -> (DVec3, f64) {
    let dir = safe_normalize(tangent);

    let weight = 0.5f64.max(dir.dot(safe_normalize(edge_vec)));
    let bz2 = (dir * edge_vec.length() * 0.5, weight);
    let t = 2.0 / 3.0;
    let bz3_v = DVec3::ZERO.lerp(bz2.0, t);
    let bz3_w = 1.0 + (bz2.1 - 1.0) * t;

    if bz3_w == 0.0 {
        (bz3_v, 0.0)
    } else {
        (bz3_v / bz3_w, bz3_w)
    }
}

struct InterpTri<'a> {
    vert_pos: &'a mut [DVec3],
    vert_bary: &'a [Barycentric],
    impl_: &'a ManifoldImpl,
}

impl<'a> InterpTri<'a> {
    fn homogeneous_from_vec4(v: (DVec3, f64)) -> (DVec3, f64) {
        (v.0 * v.1, v.1)
    }

    fn homogeneous(v: DVec3) -> (DVec3, f64) {
        (v, 1.0)
    }

    fn h_normalize(v: (DVec3, f64)) -> DVec3 {
        if v.1 == 0.0 {
            v.0
        } else {
            v.0 / v.1
        }
    }

    fn bezier(point: DVec3, tangent: (DVec3, f64)) -> (DVec3, f64) {
        let sum_xyz = point + tangent.0;
        let sum_w = tangent.1;
        (sum_xyz * sum_w, sum_w)
    }

    fn cubic_bezier_2_linear(
        p0: (DVec3, f64),
        p1: (DVec3, f64),
        p2: (DVec3, f64),
        p3: (DVec3, f64),
        x: f64,
    ) -> [(DVec3, f64); 2] {
        let lerp = |a: (DVec3, f64), b: (DVec3, f64), t: f64| {
            (a.0.lerp(b.0, t), a.1 + (b.1 - a.1) * t)
        };
        let p12 = lerp(p1, p2, x);
        let out0 = lerp(lerp(p0, p1, x), p12, x);
        let out1 = lerp(p12, lerp(p2, p3, x), x);
        [out0, out1]
    }

    fn bezier_point(points: [(DVec3, f64); 2], x: f64) -> DVec3 {
        let lerp = |a: (DVec3, f64), b: (DVec3, f64), t: f64| {
            (a.0.lerp(b.0, t), a.1 + (b.1 - a.1) * t)
        };
        Self::h_normalize(lerp(points[0], points[1], x))
    }

    fn bezier_tangent(points: [(DVec3, f64); 2]) -> DVec3 {
        safe_normalize(Self::h_normalize(points[1]) - Self::h_normalize(points[0]))
    }

    fn rotate_from_to(v: DVec3, start: DQuat, end: DQuat) -> DVec3 {
        let v_prime = start.conjugate().mul_vec3(v);
        end.mul_vec3(v_prime)
    }

    fn slerp(x: DQuat, y: DQuat, a: f64, long_way: bool) -> DQuat {
        let mut z = y;
        let mut cos_theta = x.dot(y);

        if (cos_theta < 0.0) != long_way {
            z = -y;
            cos_theta = -cos_theta;
        }

        if cos_theta > 1.0 - f64::EPSILON {
            x.lerp(z, a)
        } else {
            let angle = cos_theta.acos();
            let sin_angle = angle.sin();
            let scale_x = ((1.0 - a) * angle).sin() / sin_angle;
            let scale_z = (a * angle).sin() / sin_angle;
            DQuat::from_xyzw(
                scale_x * x.x + scale_z * z.x,
                scale_x * x.y + scale_z * z.y,
                scale_x * x.z + scale_z * z.z,
                scale_x * x.w + scale_z * z.w,
            )
        }
    }

    fn bezier_2_bezier(
        corners: [DVec3; 2],
        tangents_x: [(DVec3, f64); 2],
        tangents_y: [(DVec3, f64); 2],
        x: f64,
        anchor: DVec3,
    ) -> [(DVec3, f64); 2] {
        let bez = Self::cubic_bezier_2_linear(
            Self::homogeneous(corners[0]),
            Self::bezier(corners[0], tangents_x[0]),
            Self::bezier(corners[1], tangents_x[1]),
            Self::homogeneous(corners[1]),
            x,
        );
        let end = Self::bezier_point(bez, x);
        let tangent = Self::bezier_tangent(bez);

        let n_tangents_x = [
            safe_normalize(tangents_x[0].0),
            -safe_normalize(tangents_x[1].0),
        ];
        let bi_tangents = [
            orthogonal_to(tangents_y[0].0, anchor - corners[0], n_tangents_x[0]),
            orthogonal_to(tangents_y[1].0, anchor - corners[1], n_tangents_x[1]),
        ];

        let q0 = DQuat::from_mat3(&DMat3::from_cols(
            n_tangents_x[0],
            bi_tangents[0],
            n_tangents_x[0].cross(bi_tangents[0]),
        ));
        let q1 = DQuat::from_mat3(&DMat3::from_cols(
            n_tangents_x[1],
            bi_tangents[1],
            n_tangents_x[1].cross(bi_tangents[1]),
        ));

        let edge = corners[1] - corners[0];
        let long_way = n_tangents_x[0].dot(edge) + n_tangents_x[1].dot(edge) < 0.0;
        let q_tmp = Self::slerp(q0, q1, x, long_way);

        let q_x_dir = q_tmp.mul_vec3(DVec3::X);
        let q = DQuat::from_rotation_arc(q_x_dir, tangent) * q_tmp;

        let delta = Self::rotate_from_to(tangents_y[0].0, q0, q)
            .lerp(Self::rotate_from_to(tangents_y[1].0, q1, q), x);
        let delta_w = tangents_y[0].1 + (tangents_y[1].1 - tangents_y[0].1) * x;

        [Self::homogeneous(end), (delta, delta_w)]
    }

    fn bezier_2d(
        corners: [DVec3; 4],
        tangents_x: [(DVec3, f64); 4],
        tangents_y: [(DVec3, f64); 4],
        x: f64,
        y: f64,
        centroid: DVec3,
    ) -> DVec3 {
        let bez0 = Self::bezier_2_bezier(
            [corners[0], corners[1]],
            [tangents_x[0], tangents_x[1]],
            [tangents_y[0], tangents_y[1]],
            x,
            centroid,
        );
        let bez1 = Self::bezier_2_bezier(
            [corners[2], corners[3]],
            [tangents_x[2], tangents_x[3]],
            [tangents_y[2], tangents_y[3]],
            1.0 - x,
            centroid,
        );

        let bez = Self::cubic_bezier_2_linear(
            bez0[0],
            Self::bezier(bez0[0].0, bez0[1]),
            Self::bezier(bez1[0].0, bez1[1]),
            bez1[0],
            y,
        );
        Self::bezier_point(bez, y)
    }

    fn call(&mut self, vert: usize) {
        let tri = self.vert_bary[vert].tri as usize;
        let uvw = self.vert_bary[vert].uvw;

        let mut halfedges = [0i32; 4];
        halfedges[0] = (tri * 3) as i32;
        halfedges[1] = (tri * 3 + 1) as i32;
        halfedges[2] = (tri * 3 + 2) as i32;

        let _pair = self.impl_.halfedge[tri * 3].paired_halfedge;
        halfedges[3] = -1;

        let corners = [
            self.impl_.vert_pos[self.impl_.halfedge[halfedges[0] as usize].start_vert as usize],
            self.impl_.vert_pos[self.impl_.halfedge[halfedges[1] as usize].start_vert as usize],
            self.impl_.vert_pos[self.impl_.halfedge[halfedges[2] as usize].start_vert as usize],
            if halfedges[3] >= 0 {
                self.impl_.vert_pos[self.impl_.halfedge[halfedges[3] as usize].start_vert as usize]
            } else {
                DVec3::ZERO
            },
        ];

        for i in 0..4 {
            let _idx = if i == 3 { 3 } else { i };
            let val = match i { 0 => uvw.x, 1 => uvw.y, 2 => uvw.z, _ => uvw.w };
            if val == 1.0 {
                self.vert_pos[vert] = corners[i];
                return;
            }
        }

        let mut pos_h = (DVec3::ZERO, 0.0);

        if halfedges[3] < 0 {
            // Tri
            let tangent_r = [
                self.impl_.halfedge_tangent[halfedges[0] as usize],
                self.impl_.halfedge_tangent[halfedges[1] as usize],
                self.impl_.halfedge_tangent[halfedges[2] as usize],
            ];
            let tangent_l = [
                self.impl_.halfedge_tangent[self.impl_.halfedge[halfedges[2] as usize].paired_halfedge as usize],
                self.impl_.halfedge_tangent[self.impl_.halfedge[halfedges[0] as usize].paired_halfedge as usize],
                self.impl_.halfedge_tangent[self.impl_.halfedge[halfedges[1] as usize].paired_halfedge as usize],
            ];
            let centroid = (corners[0] + corners[1] + corners[2]) / 3.0;

            for i in 0..3 {
                let j = (i + 1) % 3;
                let k = (i + 2) % 3;
                let uvw_arr = [uvw.x, uvw.y, uvw.z];
                let x = uvw_arr[k] / (1.0 - uvw_arr[i]);

                let bez = Self::bezier_2_bezier(
                    [corners[j], corners[k]],
                    [tangent_r[j], tangent_l[k]],
                    [tangent_l[j], tangent_r[k]],
                    x,
                    centroid,
                );

                let lerp = |a: (DVec3, f64), b: (DVec3, f64), t: f64| {
                    (a.0.lerp(b.0, t), a.1 + (b.1 - a.1) * t)
                };

                let bez1 = Self::cubic_bezier_2_linear(
                    bez[0],
                    Self::bezier(bez[0].0, bez[1]),
                    Self::bezier(corners[i], lerp(tangent_r[i], tangent_l[i], x)),
                    Self::homogeneous(corners[i]),
                    uvw_arr[i],
                );
                let p = Self::bezier_point(bez1, uvw_arr[i]);
                let add = Self::homogeneous_from_vec4((p, uvw_arr[j] * uvw_arr[k]));
                pos_h.0 += add.0;
                pos_h.1 += add.1;
            }
        }
        self.vert_pos[vert] = Self::h_normalize(pos_h);
    }
}

impl ManifoldImpl {
    pub fn get_normal(&self, halfedge: usize, normal_idx: usize) -> DVec3 {
        let prop = self.halfedge[halfedge].prop_vert as usize;
        let mut normal = DVec3::ZERO;
        for i in 0..3 {
            normal[i] = self.properties[prop * self.num_prop + normal_idx + i];
        }
        normal
    }

    pub fn tangent_from_normal(&self, normal: DVec3, halfedge: usize) -> (DVec3, f64) {
        let edge = self.halfedge[halfedge];
        let edge_vec =
            self.vert_pos[edge.end_vert as usize] - self.vert_pos[edge.start_vert as usize];
        let edge_normal =
            self.face_normal[halfedge / 3] + self.face_normal[edge.paired_halfedge as usize / 3];
        let dir = edge_normal.cross(edge_vec).cross(normal);
        circular_tangent(dir, edge_vec)
    }

    pub fn is_inside_quad(&self, halfedge: usize) -> bool {
        if !self.halfedge_tangent.is_empty() {
            return self.halfedge_tangent[halfedge].1 < 0.0;
        }
        let tri = halfedge / 3;
        let ref_tri = self.mesh_relation.tri_ref[tri];
        let pair = self.halfedge[halfedge].paired_halfedge as usize;
        let pair_tri = pair / 3;
        let pair_ref = self.mesh_relation.tri_ref[pair_tri];

        if !ref_tri.same_face(&pair_ref) {
            return false;
        }

        let next_halfedge = |h: usize| -> usize {
            if h % 3 == 2 {
                h - 2
            } else {
                h + 1
            }
        };

        let same_face = |h: usize, r: &TriRef| -> bool {
            let neighbor = self.halfedge[h].paired_halfedge as usize;
            r.same_face(&self.mesh_relation.tri_ref[neighbor / 3])
        };

        let mut neighbor = next_halfedge(halfedge);
        if same_face(neighbor, &ref_tri) {
            return false;
        }
        neighbor = next_halfedge(neighbor);
        if same_face(neighbor, &ref_tri) {
            return false;
        }
        neighbor = next_halfedge(pair);
        if same_face(neighbor, &pair_ref) {
            return false;
        }
        neighbor = next_halfedge(neighbor);
        if same_face(neighbor, &pair_ref) {
            return false;
        }

        true
    }

    pub fn is_marked_inside_quad(&self, halfedge: usize) -> bool {
        !self.halfedge_tangent.is_empty() && self.halfedge_tangent[halfedge].1 < 0.0
    }

    pub fn update_sharpened_edges(&self, sharpened_edges: &[Smoothness]) -> Vec<Smoothness> {
        let mut old_halfedge2new = HashMap::new();
        for tri in 0..self.halfedge.len() / 3 {
            let old_tri = self.mesh_relation.tri_ref[tri].face_id;
            for i in 0..3 {
                old_halfedge2new.insert(3 * old_tri + i as i32, 3 * tri + i);
            }
        }

        let mut new_sharp = sharpened_edges.to_vec();
        for edge in &mut new_sharp {
            if let Some(&new_he) = old_halfedge2new.get(&edge.halfedge) {
                edge.halfedge = new_he as i32;
            }
        }
        new_sharp
    }

    pub fn flat_faces(&self) -> Vec<bool> {
        let num_tri = self.halfedge.len() / 3;
        let mut tri_is_flat_face = vec![false; num_tri];

        for tri in 0..num_tri {
            let ref_tri = self.mesh_relation.tri_ref[tri];
            let mut face_neighbors = 0;
            let mut face_tris = [-1i32; 3];

            for j in 0..3 {
                let neighbor_tri = self.halfedge[3 * tri + j].paired_halfedge as usize / 3;
                let j_ref = self.mesh_relation.tri_ref[neighbor_tri];
                if j_ref.same_face(&ref_tri) {
                    face_neighbors += 1;
                    face_tris[j] = neighbor_tri as i32;
                }
            }

            if face_neighbors > 1 {
                tri_is_flat_face[tri] = true;
                for &j in &face_tris {
                    if j >= 0 {
                        tri_is_flat_face[j as usize] = true;
                    }
                }
            }
        }

        tri_is_flat_face
    }

    pub fn vert_flat_face(&self, flat_faces: &[bool]) -> Vec<i32> {
        let num_vert = self.vert_pos.len();
        let mut vert_flat_face = vec![-1i32; num_vert];
        let mut vert_ref = vec![TriRef::default(); num_vert];

        for tri in 0..flat_faces.len() {
            if flat_faces[tri] {
                for j in 0..3 {
                    let vert = self.halfedge[3 * tri + j].start_vert as usize;
                    if !vert_ref[vert].same_face(&self.mesh_relation.tri_ref[tri]) {
                        vert_ref[vert] = self.mesh_relation.tri_ref[tri];
                        vert_flat_face[vert] = if vert_flat_face[vert] == -1 {
                            tri as i32
                        } else {
                            -2
                        };
                    }
                }
            }
        }

        vert_flat_face
    }

    pub fn vert_halfedge(&self) -> Vec<usize> {
        let num_vert = self.vert_pos.len();
        let mut vert_halfedge = vec![0usize; num_vert];
        let mut counters = vec![0u8; num_vert];

        for idx in 0..self.halfedge.len() {
            let vert = self.halfedge[idx].start_vert as usize;
            let old = std::mem::replace(&mut counters[vert], 1);
            if old == 1 {
                vert_halfedge[vert] = idx;
            }
        }

        vert_halfedge
    }

    pub fn sharpen_edges(&self, min_sharp_angle: f64, min_smoothness: f64) -> Vec<Smoothness> {
        let mut sharpened_edges = Vec::new();
        let min_radians = min_sharp_angle.to_radians();

        for e in 0..self.halfedge.len() {
            if !self.halfedge[e].is_forward() {
                continue;
            }
            let pair = self.halfedge[e].paired_halfedge as usize;
            let dihedral = self.face_normal[e / 3]
                .dot(self.face_normal[pair / 3])
                .acos();
            if dihedral > min_radians {
                sharpened_edges.push(Smoothness {
                    halfedge: e as i32,
                    smoothness: min_smoothness,
                });
                sharpened_edges.push(Smoothness {
                    halfedge: pair as i32,
                    smoothness: min_smoothness,
                });
            }
        }

        sharpened_edges
    }

    pub fn sharpen_tangent(&mut self, halfedge: usize, smoothness: f64) {
        self.halfedge_tangent[halfedge].0 *= smoothness;
        if smoothness == 0.0 {
            self.halfedge_tangent[halfedge].1 = 0.0;
        }
    }

    pub fn set_normals(&mut self, normal_idx: i32, min_sharp_angle: f64) {
        if self.vert_pos.is_empty() || normal_idx < 0 {
            return;
        }

        let old_num_prop = self.num_prop;
        let tri_is_flat_face = self.flat_faces();
        let vert_flat_face = self.vert_flat_face(&tri_is_flat_face);
        let mut vert_num_sharp = vec![0i32; self.vert_pos.len()];

        for e in 0..self.halfedge.len() {
            if !self.halfedge[e].is_forward() {
                continue;
            }
            let pair = self.halfedge[e].paired_halfedge as usize;
            let tri1 = e / 3;
            let tri2 = pair / 3;
            let dihedral = self.face_normal[tri1]
                .dot(self.face_normal[tri2])
                .acos()
                .to_degrees();
            if dihedral > min_sharp_angle {
                vert_num_sharp[self.halfedge[e].start_vert as usize] += 1;
                vert_num_sharp[self.halfedge[e].end_vert as usize] += 1;
            } else {
                let face_split = tri_is_flat_face[tri1] != tri_is_flat_face[tri2]
                    || (tri_is_flat_face[tri1]
                        && tri_is_flat_face[tri2]
                        && !self.mesh_relation.tri_ref[tri1]
                            .same_face(&self.mesh_relation.tri_ref[tri2]));
                if vert_flat_face[self.halfedge[e].start_vert as usize] == -2 && face_split {
                    vert_num_sharp[self.halfedge[e].start_vert as usize] += 1;
                }
                if vert_flat_face[self.halfedge[e].end_vert as usize] == -2 && face_split {
                    vert_num_sharp[self.halfedge[e].end_vert as usize] += 1;
                }
            }
        }

        let num_prop = old_num_prop.max(normal_idx as usize + 3);
        let mut old_properties = vec![0.0f64; num_prop * self.num_prop_vert()];
        std::mem::swap(&mut self.properties, &mut old_properties);
        self.num_prop = num_prop;

        let mut old_halfedge_prop = vec![0i32; self.halfedge.len()];
        for i in 0..self.halfedge.len() {
            old_halfedge_prop[i] = self.halfedge[i].prop_vert;
            self.halfedge[i].prop_vert = -1;
        }

        let num_edge = self.halfedge.len();
        for start_edge in 0..num_edge {
            if self.halfedge[start_edge].prop_vert >= 0 {
                continue;
            }
            let vert = self.halfedge[start_edge].start_vert as usize;

            if vert_num_sharp[vert] < 2 {
                // vertex has single normal
                let normal = if vert_flat_face[vert] >= 0 {
                    self.face_normal[vert_flat_face[vert] as usize]
                } else {
                    self.vert_normal[vert]
                };
                let mut last_prop = -1;
                let vert_edges = self.collect_vert_edges(start_edge);
                for current in vert_edges {
                    let prop = old_halfedge_prop[current];
                    self.halfedge[current].prop_vert = prop;
                    if prop != last_prop {
                        last_prop = prop;
                        let start = prop as usize * old_num_prop;
                        let dest = prop as usize * num_prop;
                        for i in 0..old_num_prop {
                            self.properties[dest + i] = old_properties[start + i];
                        }
                        for i in 0..3 {
                            self.properties[dest + normal_idx as usize + i] = normal[i];
                        }
                    }
                }
            } else {
                // vertex has multiple normals
                let center_pos = self.vert_pos[vert];
                let mut group = Vec::new();
                let mut normals = Vec::new();
                let mut current = start_edge;
                let mut prev_face = current / 3;

                loop {
                    let next = crate::kernel::next_halfedge(self.halfedge[current].paired_halfedge as i32) as usize;
                    let face = next / 3;

                    let dihedral = self.face_normal[face].dot(self.face_normal[prev_face]).acos().to_degrees();
                    if dihedral > min_sharp_angle || tri_is_flat_face[face] != tri_is_flat_face[prev_face] || (tri_is_flat_face[face] && tri_is_flat_face[prev_face] && !self.mesh_relation.tri_ref[face].same_face(&self.mesh_relation.tri_ref[prev_face])) {
                        break;
                    }
                    current = next;
                    prev_face = face;
                    if current == start_edge {
                        break;
                    }
                }

                let end_edge = current;
                let vert_edges = self.collect_vert_edges(end_edge);

                struct FaceEdge {
                    face: usize,
                    edge_vec: DVec3,
                }

                let mut face_edges = Vec::new();
                for &curr in &vert_edges {
                    if self.is_inside_quad(curr) {
                        face_edges.push(FaceEdge { face: curr/3, edge_vec: DVec3::NAN });
                    } else {
                        let v = self.halfedge[curr].end_vert as usize;
                        let mut pos = self.vert_pos[v];
                        if vert_num_sharp[v] < 2 {
                            let normal = if vert_flat_face[v] >= 0 {
                                self.face_normal[vert_flat_face[v] as usize]
                            } else {
                                self.vert_normal[v]
                            };
                            let tangent = self.tangent_from_normal(normal, self.halfedge[curr].paired_halfedge as usize);
                            pos += tangent.0;
                        }
                        face_edges.push(FaceEdge { face: curr/3, edge_vec: safe_normalize(pos - center_pos) });
                    }
                }

                for i in 0..face_edges.len() {
                    let here = &face_edges[i];
                    let next = &face_edges[(i + 1) % face_edges.len()];

                    let dihedral = self.face_normal[here.face].dot(self.face_normal[next.face]).acos().to_degrees();
                    if dihedral > min_sharp_angle || tri_is_flat_face[here.face] != tri_is_flat_face[next.face] || (tri_is_flat_face[here.face] && tri_is_flat_face[next.face] && !self.mesh_relation.tri_ref[here.face].same_face(&self.mesh_relation.tri_ref[next.face])) {
                        normals.push(DVec3::ZERO);
                    }
                    group.push(normals.len() - 1);
                    if next.edge_vec.x.is_finite() {
                        let cross = next.edge_vec.cross(here.edge_vec);
                        normals.last_mut().unwrap().add_assign(safe_normalize(cross) * angle_between(here.edge_vec, next.edge_vec));
                    }
                }

                for n in &mut normals {
                    *n = safe_normalize(*n);
                }

                let mut last_group = 0;
                let mut last_prop = -1;

                for (idx, &current1) in vert_edges.iter().enumerate() {
                    let prop = old_halfedge_prop[current1];
                    let start = prop as usize * old_num_prop;

                    if group[idx] != last_group && group[idx] != 0 && prop == last_prop {
                        last_group = group[idx];
                        let new_prop = self.num_prop_vert();
                        let mut new_props = vec![0.0; num_prop];
                        for i in 0..old_num_prop {
                            new_props[i] = old_properties[start + i];
                        }
                        for i in 0..3 {
                            new_props[normal_idx as usize + i] = normals[group[idx]][i];
                        }
                        self.properties.extend(new_props);
                        self.halfedge[current1].prop_vert = new_prop as i32;
                    } else if prop != last_prop {
                        last_prop = prop;
                        let dest = prop as usize * num_prop;
                        for i in 0..old_num_prop {
                            self.properties[dest + i] = old_properties[start + i];
                        }
                        for i in 0..3 {
                            self.properties[dest + normal_idx as usize + i] = normals[group[idx]][i];
                        }
                        self.halfedge[current1].prop_vert = prop;
                    } else {
                        self.halfedge[current1].prop_vert = prop;
                    }
                }
            }
        }
    }

    pub fn linearize_flat_tangents(&mut self) {
        let n = self.halfedge_tangent.len();
        for halfedge in 0..n {
            let paired = self.halfedge[halfedge].paired_halfedge;
            if paired < 0 {
                continue;
            }
            let paired_idx = paired as usize;

            let mut t1 = self.halfedge_tangent[halfedge];
            let mut t2 = self.halfedge_tangent[paired_idx];

            let flat = [t1.1 == 0.0, t2.1 == 0.0];
            if !self.halfedge[halfedge].is_forward() || (!flat[0] && !flat[1]) {
                continue;
            }

            let edge_vec = self.vert_pos[self.halfedge[halfedge].end_vert as usize]
                - self.vert_pos[self.halfedge[halfedge].start_vert as usize];

            if flat[0] && flat[1] {
                t1.0 = edge_vec / 3.0;
                t1.1 = 1.0;
                t2.0 = -edge_vec / 3.0;
                t2.1 = 1.0;
            } else if flat[0] {
                t1.0 = (edge_vec + t2.0) * 0.5;
                t1.1 = 1.0;
            } else {
                t2.0 = (-edge_vec + t1.0) * 0.5;
                t2.1 = 1.0;
            }

            self.halfedge_tangent[halfedge] = t1;
            self.halfedge_tangent[paired_idx] = t2;
        }
    }

    pub fn distribute_tangents(&mut self, fixed_halfedges: &[bool]) {
        let num_halfedge = fixed_halfedges.len();
        for halfedge in 0..num_halfedge {
            if !fixed_halfedges[halfedge] {
                continue;
            }

            let mut current = halfedge;
            if self.is_marked_inside_quad(halfedge) {
                current = self.halfedge[current].paired_halfedge as usize + 1;
            }

            let approx_normal = self.vert_normal[self.halfedge[current].start_vert as usize];
            let center = self.vert_pos[self.halfedge[current].start_vert as usize];
            let mut last_edge_vec =
                safe_normalize(self.vert_pos[self.halfedge[current].end_vert as usize] - center);
            let first_tangent = safe_normalize(self.halfedge_tangent[current].0);
            let mut last_tangent = first_tangent;
            let mut normal = DVec3::ZERO;
            let mut current_angle = Vec::new();
            let mut desired_angle = Vec::new();

            loop {
                current = self.halfedge[current].paired_halfedge as usize;
                if current % 3 == 2 {
                    current -= 2;
                } else {
                    current += 1;
                }

                if self.is_marked_inside_quad(current) {
                    continue;
                }

                let this_edge_vec = safe_normalize(
                    self.vert_pos[self.halfedge[current].end_vert as usize] - center,
                );
                let this_tangent = safe_normalize(self.halfedge_tangent[current].0);
                normal += this_tangent.cross(last_tangent);

                let prev_angle = desired_angle.last().copied().unwrap_or(0.0);
                desired_angle.push(angle_between(this_edge_vec, last_edge_vec) + prev_angle);

                if current == halfedge {
                    current_angle.push(K_TWO_PI);
                } else {
                    current_angle.push(angle_between(this_tangent, first_tangent));
                    if approx_normal.dot(this_tangent.cross(first_tangent)) < 0.0 {
                        if let Some(v) = current_angle.last_mut() {
                            *v = K_TWO_PI - *v;
                        }
                    }
                }

                last_edge_vec = this_edge_vec;
                last_tangent = this_tangent;

                if fixed_halfedges[current] {
                    break;
                }
            }

            if current_angle.len() == 1 || normal.dot(normal) == 0.0 {
                continue;
            }

            let scale = *current_angle.last().unwrap() / *desired_angle.last().unwrap();
            let mut offset = 0.0;

            if current == halfedge {
                for i in 0..current_angle.len() {
                    offset += wrap(current_angle[i] - scale * desired_angle[i]);
                }
                offset /= current_angle.len() as f64;
            }

            let mut i = 0;
            current = halfedge;
            loop {
                current = self.halfedge[current].paired_halfedge as usize;
                if current % 3 == 2 {
                    current -= 2;
                } else {
                    current += 1;
                }

                if self.is_marked_inside_quad(current) {
                    continue;
                }

                desired_angle[i] *= scale;
                let last_angle = if i > 0 { desired_angle[i - 1] } else { 0.0 };

                if desired_angle[i] - last_angle > K_PI {
                    desired_angle[i] = last_angle + K_PI;
                } else if i + 1 < desired_angle.len()
                    && scale * desired_angle[i + 1] - desired_angle[i] > K_PI
                {
                    desired_angle[i] = scale * desired_angle[i + 1] - K_PI;
                }

                let angle = current_angle[i] - desired_angle[i] - offset;
                let mut tangent = self.halfedge_tangent[current];
                let q = DQuat::from_axis_angle(normal.normalize(), angle);
                tangent.0 = q.mul_vec3(tangent.0);
                self.halfedge_tangent[current] = tangent;
                i += 1;

                if fixed_halfedges[current] {
                    break;
                }
            }
        }
    }

    pub fn create_tangents(&mut self, normal_idx: i32) {
        let num_vert = self.vert_pos.len();
        let num_halfedge = self.halfedge.len();
        self.halfedge_tangent.clear();
        self.halfedge_tangent.resize(num_halfedge, (DVec3::ZERO, 0.0));
        let mut fixed_halfedge = vec![false; num_halfedge];

        let vert_halfedge = self.vert_halfedge();

        for e in 0..num_vert {
            struct FlatNormal {
                is_flat_face: bool,
                normal: DVec3,
            }
            let mut face_edges = [-1i32; 2];
            let mut flat_normals = Vec::new();

            let vert_edges = self.collect_vert_edges(vert_halfedge[e]);
            for &halfedge in &vert_edges {
                let normal = self.get_normal(halfedge, normal_idx as usize);
                let diff = self.face_normal[halfedge / 3] - normal;
                flat_normals.push(FlatNormal {
                    is_flat_face: diff.dot(diff) < k_precision() * k_precision(),
                    normal,
                });
            }

            for (idx, &halfedge) in vert_edges.iter().enumerate() {
                let here = &flat_normals[idx];
                let next = &flat_normals[(idx + 1) % flat_normals.len()];

                if self.is_inside_quad(halfedge) {
                    self.halfedge_tangent[halfedge] = (DVec3::ZERO, -1.0);
                    continue;
                }

                let diff = next.normal - here.normal;
                let different_normals = diff.dot(diff) > k_precision() * k_precision();

                if different_normals || here.is_flat_face != next.is_flat_face {
                    fixed_halfedge[halfedge] = true;
                    if face_edges[0] == -1 {
                        face_edges[0] = halfedge as i32;
                    } else if face_edges[1] == -1 {
                        face_edges[1] = halfedge as i32;
                    } else {
                        face_edges[0] = -2;
                    }
                }

                if different_normals {
                    let edge_vec = self.vert_pos[self.halfedge[halfedge].end_vert as usize] - self.vert_pos[self.halfedge[halfedge].start_vert as usize];
                    let dir = here.normal.cross(next.normal);
                    let sign = if dir.dot(edge_vec) < 0.0 { -1.0 } else { 1.0 };
                    self.halfedge_tangent[halfedge] = circular_tangent(dir * sign, edge_vec);
                } else {
                    self.halfedge_tangent[halfedge] = self.tangent_from_normal(here.normal, halfedge);
                }
            }

            if face_edges[0] >= 0 && face_edges[1] >= 0 {
                let edge0 = self.vert_pos[self.halfedge[face_edges[0] as usize].end_vert as usize] - self.vert_pos[self.halfedge[face_edges[0] as usize].start_vert as usize];
                let edge1 = self.vert_pos[self.halfedge[face_edges[1] as usize].end_vert as usize] - self.vert_pos[self.halfedge[face_edges[1] as usize].start_vert as usize];
                let new_tangent = safe_normalize(edge0) - safe_normalize(edge1);
                self.halfedge_tangent[face_edges[0] as usize] = circular_tangent(new_tangent, edge0);
                self.halfedge_tangent[face_edges[1] as usize] = circular_tangent(-new_tangent, edge1);
            } else if face_edges[0] == -1 && face_edges[1] == -1 {
                fixed_halfedge[vert_halfedge[e]] = true;
            }
        }

        self.distribute_tangents(&fixed_halfedge);
    }

    pub fn refine<F>(&mut self, edge_divisions: F, _keep_interior: bool)
    where
        F: Fn(DVec3, DVec4, DVec4) -> i32,
    {
        if self.vert_pos.is_empty() {
            return;
        }

        let old = self.clone();

        let vert_bary = self.subdivide(edge_divisions);

        if vert_bary.is_empty() {
            return;
        }

        if old.halfedge_tangent.len() == old.halfedge.len() {
            let len = self.vert_pos.len();
            let mut interp = InterpTri {
                vert_pos: &mut self.vert_pos,
                vert_bary: &vert_bary,
                impl_: &old,
            };

            for i in 0..len {
                interp.call(i);
            }
        }

        self.halfedge_tangent.clear();
        self.finish();
        if old.halfedge_tangent.len() == old.halfedge.len() {
            self.mark_coplanar();
        }
        self.mesh_relation.original_id = -1;
    }

    // Helper to collect edges around a vertex
    fn collect_vert_edges(&self, start_edge: usize) -> Vec<usize> {
        let mut edges = Vec::new();
        let mut current = start_edge;
        loop {
            edges.push(current);
            let paired = self.halfedge[current].paired_halfedge;
            if paired < 0 {
                break;
            }
            let paired = paired as usize;
            if paired % 3 == 2 {
                current = paired - 2;
            } else {
                current = paired + 1;
            }
            if current == start_edge {
                break;
            }
        }
        edges
    }
}

trait AddAssign {
    fn add_assign(&mut self, other: Self);
}

impl AddAssign for DVec3 {
    fn add_assign(&mut self, other: DVec3) {
        *self = *self + other;
    }
}
