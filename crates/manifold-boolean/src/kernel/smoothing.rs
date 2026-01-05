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

use glam::{DQuat, DVec3};
use manifold_math::{degrees, radians, K_PI, K_TWO_PI};
use manifold_types::{Barycentric, Halfedge, TriRef};

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
    1e-6
}

use super::ManifoldImpl;

/// Returns a normalized vector orthogonal to ref, in the plane of ref and in,
/// unless in and ref are colinear, in which case it falls back to the plane of
/// ref and altIn.
fn orthogonal_to(in_vec: DVec3, alt_in: DVec3, ref_vec: DVec3) -> DVec3 {
    let dot = in_vec.dot(ref_vec);
    let mut out = in_vec - dot * ref_vec;
    let k_precision = 1e-6;
    if out.dot(out) < k_precision * in_vec.dot(in_vec) {
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
fn circular_tangent(tangent: DVec3, edge_vec: DVec3) -> DVec3 {
    let dir = safe_normalize(tangent);

    let weight = 0.5f64.max(dir.normalize_or_zero().dot(safe_normalize(edge_vec)));
    // Quadratic weighted bezier for circular interpolation
    let bz2 = (dir * edge_vec.length() * 0.5, weight);
    // Equivalent cubic weighted bezier
    let bz3 = (
        glam::DVec3::ZERO.lerp(DVec3::from(bz2.0), 2.0 / 3.0),
        bz2.1 * (2.0 / 3.0),
    );
    // Convert from homogeneous form to geometric form
    DVec3::new(bz3.0.x / bz3.1, bz3.0.y / bz3.1, bz3.0.z / bz3.1)
}

struct InterpTri<'a> {
    vert_pos: &'a mut [DVec3],
    vert_bary: &'a [Barycentric],
    impl_: &'a ManifoldImpl,
}

impl<'a> InterpTri<'a> {
    fn homogeneous(v4: (DVec3, f64)) -> (DVec3, f64) {
        let v3 = v4.0 * v4.1;
        (v3, v4.1)
    }

    fn homogeneous_from_vec3(v: DVec3) -> (DVec3, f64) {
        (v, 1.0)
    }

    fn h_normalize(v4: (DVec3, f64)) -> DVec3 {
        if v4.1 == 0.0 {
            v4.0
        } else {
            v4.0 / v4.1
        }
    }

    fn scale(v4: (DVec3, f64), scale: f64) -> (DVec3, f64) {
        (v4.0 * scale, v4.1)
    }

    fn bezier(point: DVec3, tangent: (DVec3, f64)) -> (DVec3, f64) {
        let point_v4 = (point, 0.0);
        (point_v4.0 + tangent.0, tangent.1)
    }

    fn cubic_bezier_2_linear(
        p0: (DVec3, f64),
        p1: (DVec3, f64),
        p2: (DVec3, f64),
        p3: (DVec3, f64),
        x: f64,
    ) -> [(DVec3, f64); 2] {
        let p12 = lerp_v4(&p1, &p2, x);
        let out_0 = lerp_v4(&lerp_v4(&p0, &p1, x), &p12, x);
        let out_1 = lerp_v4(&p12, &lerp_v4(&p2, &p3, x), x);
        [out_0, out_1]
    }

    fn bezier_point(points: [(DVec3, f64); 2], x: f64) -> DVec3 {
        Self::h_normalize(lerp_v4(&points[0], &points[1], x))
    }

    fn bezier_tangent(points: [(DVec3, f64); 2]) -> DVec3 {
        let p0 = Self::h_normalize(points[0]);
        let p1 = Self::h_normalize(points[1]);
        safe_normalize(p1 - p0)
    }

    fn rotate_from_to(v: DVec3, start: glam::DQuat, end: glam::DQuat) -> DVec3 {
        let q = start.inverse() * end;
        q.mul_vec3(v)
    }

    fn slerp(x: glam::DQuat, y: glam::DQuat, a: f64, long_way: bool) -> glam::DQuat {
        let mut z = y;
        let mut cos_theta = x.dot(y);

        // Take the long way around the sphere only when requested
        if (cos_theta < 0.0) != long_way {
            z = -y;
            cos_theta = -cos_theta;
        }

        if cos_theta > 1.0 - f64::EPSILON {
            // for numerical stability use standard slerp
            x.slerp(z, a as f32 as f64)
        } else {
            let angle = cos_theta.acos();
            let sin_angle = angle.sin();
            let scale_x = ((1.0 - a) * angle).sin() / sin_angle;
            let scale_z = (a * angle).sin() / sin_angle;
            // Manual quaternion scaling and addition
            glam::DQuat::from_xyzw(
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
            Self::homogeneous_from_vec3(corners[0]),
            Self::bezier(corners[0], tangents_x[0]),
            Self::bezier(corners[1], tangents_x[1]),
            Self::homogeneous_from_vec3(corners[1]),
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

        let q0 = glam::DQuat::from_mat3(&glam::DMat3::from_cols(
            n_tangents_x[0],
            bi_tangents[0],
            n_tangents_x[0].cross(bi_tangents[0]),
        ));
        let q1 = glam::DQuat::from_mat3(&glam::DMat3::from_cols(
            n_tangents_x[1],
            bi_tangents[1],
            n_tangents_x[1].cross(bi_tangents[1]),
        ));
        let edge = corners[1] - corners[0];
        let long_way = n_tangents_x[0].dot(edge) + n_tangents_x[1].dot(edge) < 0.0;
        let q_tmp = Self::slerp(q0, q1, x, long_way);

        // Create rotation quaternion that aligns q_tmp's x-axis with the tangent
        let x_axis = q_tmp.mul_vec3(glam::DVec3::X);
        let q = glam::DQuat::from_rotation_arc(x_axis, tangent) * q_tmp;

        let delta = Self::rotate_from_to(tangents_y[0].0, q0, q)
            .lerp(Self::rotate_from_to(tangents_y[1].0, q1, q), x);
        let delta_w = tangents_y[0].1 + (tangents_y[1].1 - tangents_y[0].1) * x;

        [Self::homogeneous_from_vec3(end), (delta, delta_w)]
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

    fn interpolate(&mut self, vert: usize) {
        let pos = &mut self.vert_pos[vert];
        let tri = self.vert_bary[vert].tri as usize;
        let uvw = self.vert_bary[vert].uvw;

        let halfedges = [
            tri * 3,
            tri * 3 + 1,
            tri * 3 + 2,
            if self.impl_.halfedge[tri * 3].paired_halfedge >= 0 {
                self.impl_.halfedge[tri * 3].paired_halfedge as usize
            } else {
                usize::MAX
            },
        ];

        let corners = [
            self.impl_.vert_pos[self.impl_.halfedge[halfedges[0]].start_vert as usize],
            self.impl_.vert_pos[self.impl_.halfedge[halfedges[1]].start_vert as usize],
            self.impl_.vert_pos[self.impl_.halfedge[halfedges[2]].start_vert as usize],
            if halfedges[3] != usize::MAX {
                self.impl_.vert_pos[self.impl_.halfedge[halfedges[3]].start_vert as usize]
            } else {
                DVec3::ZERO
            },
        ];

        for i in 0..4 {
            if uvw[i] == 1.0 {
                *pos = corners[i];
                return;
            }
        }

        let mut pos_h = (DVec3::ZERO, 0.0);

        if halfedges[3] == usize::MAX {
            // tri
            let tangent_r = [
                (
                    self.impl_.halfedge_tangent[halfedges[0]].0,
                    self.impl_.halfedge_tangent[halfedges[0]].1,
                ),
                (
                    self.impl_.halfedge_tangent[halfedges[1]].0,
                    self.impl_.halfedge_tangent[halfedges[1]].1,
                ),
                (
                    self.impl_.halfedge_tangent[halfedges[2]].0,
                    self.impl_.halfedge_tangent[halfedges[2]].1,
                ),
            ];
            let tangent_l = [
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[2]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[2]].paired_halfedge as usize]
                        .1,
                ),
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[0]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[0]].paired_halfedge as usize]
                        .1,
                ),
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[1]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[1]].paired_halfedge as usize]
                        .1,
                ),
            ];
            let centroid = (corners[0] + corners[1] + corners[2]) / 3.0;

            for i in 0..3 {
                let j = (i + 1) % 3;
                let k = (i + 2) % 3;
                let x = uvw[k] / (1.0 - uvw[i]);

                let bez = Self::bezier_2_bezier(
                    [corners[j], corners[k]],
                    [tangent_r[j], tangent_l[k]],
                    [tangent_l[j], tangent_r[k]],
                    x,
                    centroid,
                );

                let bez1 = Self::cubic_bezier_2_linear(
                    bez[0],
                    Self::bezier(bez[0].0, bez[1]),
                    Self::bezier(
                        corners[i],
                        (
                            (tangent_r[i].0 * (1.0 - x) + tangent_l[i].0 * x),
                            tangent_r[i].1 * (1.0 - x) + tangent_l[i].1 * x,
                        ),
                    ),
                    Self::homogeneous_from_vec3(corners[i]),
                    uvw[i],
                );
                let p = Self::bezier_point(bez1, uvw[i]);
                pos_h.0 += p * (uvw[j] * uvw[k]);
            }
        } else {
            // quad
            let tangents_x = [
                (
                    self.impl_.halfedge_tangent[halfedges[0]].0,
                    self.impl_.halfedge_tangent[halfedges[0]].1,
                ),
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[0]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[0]].paired_halfedge as usize]
                        .1,
                ),
                (
                    self.impl_.halfedge_tangent[halfedges[2]].0,
                    self.impl_.halfedge_tangent[halfedges[2]].1,
                ),
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[2]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[2]].paired_halfedge as usize]
                        .1,
                ),
            ];
            let tangents_y = [
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[3]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[3]].paired_halfedge as usize]
                        .1,
                ),
                (
                    self.impl_.halfedge_tangent[halfedges[1]].0,
                    self.impl_.halfedge_tangent[halfedges[1]].1,
                ),
                (
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[1]].paired_halfedge as usize]
                        .0,
                    self.impl_.halfedge_tangent
                        [self.impl_.halfedge[halfedges[1]].paired_halfedge as usize]
                        .1,
                ),
                (
                    self.impl_.halfedge_tangent[halfedges[3]].0,
                    self.impl_.halfedge_tangent[halfedges[3]].1,
                ),
            ];
            let centroid = (corners[0] + corners[1] + corners[2] + corners[3]) * 0.25;
            let x = uvw[1] + uvw[2];
            let y = uvw[2] + uvw[3];
            let p_x = Self::bezier_2d(corners, tangents_x, tangents_y, x, y, centroid);
            let p_y = Self::bezier_2d(
                [corners[1], corners[2], corners[3], corners[0]],
                [tangents_y[1], tangents_y[2], tangents_y[3], tangents_y[0]],
                [tangents_x[1], tangents_x[2], tangents_x[3], tangents_x[0]],
                y,
                1.0 - x,
                centroid,
            );
            pos_h.0 += p_x * (x * (1.0 - x));
            pos_h.0 += p_y * (y * (1.0 - y));
        }
        *pos = Self::h_normalize(pos_h);
    }
}

fn lerp_v4(a: &(DVec3, f64), b: &(DVec3, f64), t: f64) -> (DVec3, f64) {
    (a.0.lerp(b.0, t), a.1 + (b.1 - a.1) * t)
}

impl ManifoldImpl {
    /// Get the property normal associated with the startVert of this halfedge, where
    /// normalIdx shows the beginning of where normals are stored in the properties.
    pub fn get_normal(&self, halfedge: usize, normal_idx: usize) -> DVec3 {
        let prop = self.halfedge[halfedge].prop_vert as usize;
        let mut normal = DVec3::ZERO;
        for i in 0..3 {
            normal[i] = self.properties[prop * self.num_prop + normal_idx + i];
        }
        normal
    }

    /// Returns a circular tangent for the requested halfedge, orthogonal to the
    /// given normal vector, and avoiding folding.
    pub fn tangent_from_normal(&self, normal: DVec3, halfedge: usize) -> (DVec3, f64) {
        let edge = self.halfedge[halfedge];
        let edge_vec =
            self.vert_pos[edge.end_vert as usize] - self.vert_pos[edge.start_vert as usize];
        let edge_normal =
            self.face_normal[halfedge / 3] + self.face_normal[edge.paired_halfedge as usize / 3];
        let dir = edge_normal.cross(edge_vec).cross(normal);
        let tangent = circular_tangent(dir, edge_vec);
        (tangent, 1.0)
    }

    /// Returns true if this halfedge should be marked as the interior of a quad, as
    /// defined by its two triangles referring to the same face, and those triangles
    /// having no further face neighbors beyond.
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

        // Check neighbors
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

    /// Returns true if this halfedge is an interior of a quad, as defined by its
    /// halfedge tangent having negative weight.
    pub fn is_marked_inside_quad(&self, halfedge: usize) -> bool {
        !self.halfedge_tangent.is_empty() && self.halfedge_tangent[halfedge].1 < 0.0
    }

    /// Sharpened edges are referenced to the input Mesh, but the triangles have
    /// been sorted in creating the Manifold, so the indices are converted using
    /// meshRelation_.faceID, which temporarily holds the mapping.
    pub fn update_sharpened_edges(&self, sharpened_edges: &[Smoothness]) -> Vec<Smoothness> {
        use std::collections::HashMap;

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

    /// Find faces containing at least 3 triangles - these will not have
    /// interpolated normals - all their vert normals must match their face normal.
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

    /// Returns a vector of length numVert that has a tri that is part of a
    /// neighboring flat face if there is only one flat face. If there are none it
    /// gets -1, and if there are more than one it gets -2.
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
                // arbitrary, last one wins
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

    /// Sharpen tangents that intersect an edge to sharpen that edge. The weight is
    /// unchanged, as this has a squared effect on radius of curvature, except
    /// in the case of zero radius, which is marked with weight = 0.
    pub fn sharpen_tangent(&mut self, halfedge: usize, smoothness: f64) {
        self.halfedge_tangent[halfedge].0 *= smoothness;
        if smoothness == 0.0 {
            self.halfedge_tangent[halfedge].1 = 0.0;
        }
    }

    /// Instead of calculating the internal shared normals like CalculateNormals
    /// does, this method fills in vertex properties, unshared across edges that
    /// are bent more than minSharpAngle.
    pub fn set_normals(&mut self, normal_idx: i32, min_sharp_angle: f64) {
        if self.vert_pos.is_empty() {
            return;
        }
        if normal_idx < 0 {
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
                let mut current = start_edge;
                loop {
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

                    let paired = self.halfedge[current].paired_halfedge;
                    if paired < 0 {
                        break;
                    }
                    let paired_idx = paired as usize;
                    if paired_idx % 3 == 2 {
                        current = paired_idx - 2;
                    } else {
                        current = paired_idx + 1;
                    }
                    if current == start_edge {
                        break;
                    }
                }
            } else {
                // vertex has multiple normals - simplified stub
                // Full implementation requires complex pseudo-normal calculation
            }
        }
    }

    /// Tangents get flattened to create sharp edges by setting their weight to zero.
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

    /// Redistribute the tangents around each vertex so that the angles between them
    /// have the same ratios as the angles of the triangles between the corresponding
    /// edges.
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

                // cumulative sum
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
                // only one - find average offset
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

                // shrink obtuse angles
                if desired_angle[i] - last_angle > K_PI {
                    desired_angle[i] = last_angle + K_PI;
                } else if i + 1 < desired_angle.len()
                    && scale * desired_angle[i + 1] - desired_angle[i] > K_PI
                {
                    desired_angle[i] = scale * desired_angle[i + 1] - K_PI;
                }

                let angle = current_angle[i] - desired_angle[i] - offset;
                let mut tangent = self.halfedge_tangent[current];
                let q = glam::DQuat::from_axis_angle(normal.normalize(), angle);
                tangent.0 = q.mul_vec3(tangent.0);
                self.halfedge_tangent[current] = tangent;
                i += 1;

                if fixed_halfedges[current] {
                    break;
                }
            }
        }
    }

    /// Helper function to iterate around a vertex
    fn for_vert<F>(&self, start_edge: usize, mut callback: F)
    where
        F: FnMut(usize),
    {
        let mut current = start_edge;
        loop {
            callback(current);
            let paired = self.halfedge[current].paired_halfedge as usize;
            if paired < 0 {
                break;
            }
            if paired % 3 == 2 {
                current = paired - 2;
            } else {
                current = paired + 1;
            }
            if current == start_edge {
                break;
            }
        }
    }
}
