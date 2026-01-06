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

use crate::kernel::{Intersections, ManifoldImpl};
use glam::DVec3;
use manifold_math::Box as GeoBox;
use manifold_types::{Halfedge, OpType, TriRef};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy)]
pub struct EdgePos {
    pub edge_pos: f64,
    pub vert: i32,
    pub collision_id: i32,
    pub is_start: bool,
}

impl PartialEq for EdgePos {
    fn eq(&self, other: &Self) -> bool {
        (self.edge_pos - other.edge_pos).abs() < f64::EPSILON
            && self.collision_id == other.collision_id
    }
}

impl Eq for EdgePos {}

impl PartialOrd for EdgePos {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.edge_pos.partial_cmp(&other.edge_pos) {
            Some(std::cmp::Ordering::Equal) => self.collision_id.partial_cmp(&other.collision_id),
            ord => ord,
        }
    }
}

impl Ord for EdgePos {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}

pub struct DuplicateVerts<'a> {
    pub vert_pos_r: Vec<DVec3>,
    pub inclusion: &'a Vec<i32>,
    pub vert_r: &'a Vec<i32>,
    pub vert_pos_p: &'a Vec<DVec3>,
}

impl<'a> DuplicateVerts<'a> {
    pub fn new(
        n: usize,
        inclusion: &'a Vec<i32>,
        vert_r: &'a Vec<i32>,
        vert_pos_p: &'a Vec<DVec3>,
    ) -> Self {
        Self {
            vert_pos_r: Vec::with_capacity(n),
            inclusion,
            vert_r,
            vert_pos_p,
        }
    }

    pub fn populate(&mut self) {
        for i in 0..self.vert_r.len() {
            let vert = self.vert_r[i];
            if vert >= 0 {
                let abs_inclusion = self.inclusion[i].abs() as usize;
                for _ in 0..abs_inclusion {
                    self.vert_pos_r.push(self.vert_pos_p[i]);
                }
            }
        }
    }
}

#[derive(Debug)]
pub struct BooleanResult {
    pub mesh_r: ManifoldImpl,
}

pub fn pair_up(edge_pos: &mut [EdgePos]) -> Vec<Halfedge> {
    debug_assert!(edge_pos.len() % 2 == 0);
    let n_edges = edge_pos.len() / 2;
    let (mut starts, mut ends): (Vec<_>, Vec<_>) =
        edge_pos.iter().copied().partition(|x| x.is_start);
    starts.sort_by(|a: &EdgePos, b: &EdgePos| a.partial_cmp(b).unwrap());
    ends.sort_by(|a: &EdgePos, b: &EdgePos| a.partial_cmp(b).unwrap());

    let mut edges = Vec::with_capacity(n_edges);
    for i in 0..n_edges {
        edges.push(Halfedge {
            start_vert: starts[i].vert,
            end_vert: ends[i].vert,
            paired_halfedge: -1,
            prop_vert: starts[i].vert,
        });
    }
    edges
}

pub fn add_new_edge_verts(
    edges_p: &mut BTreeMap<i32, Vec<EdgePos>>,
    _edges_new: &mut BTreeMap<(i32, i32), Vec<EdgePos>>,
    p1q2: &[[i32; 2]],
    i12: &[i32],
    v12_r: &[i32],
    halfedge_p: &[Halfedge],
    forward: bool,
    offset: usize,
) {
    for i in 0..p1q2.len() {
        let edge_p = p1q2[i][if forward { 0 } else { 1 }];
        let face_q = p1q2[i][if forward { 1 } else { 0 }];
        let vert = v12_r[i];
        let inclusion = i12[i];

        let halfedge = halfedge_p[edge_p as usize];
        let mut key_right = (halfedge.paired_halfedge / 3, face_q);
        let mut key_left = (edge_p / 3, face_q);
        if !forward {
            std::mem::swap(&mut key_right, &mut key_left);
        }

        let direction = inclusion < 0;

        let edge_pos = EdgePos {
            edge_pos: 0.0,
            vert: vert,
            collision_id: (i + offset) as i32,
            is_start: direction,
        };
        edges_p.entry(edge_p).or_default().push(edge_pos);
    }
}

pub fn size_output(
    out_r: &mut ManifoldImpl,
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    i03: &[i32],
    i30: &[i32],
    i12: &[i32],
    i21: &[i32],
    p1q2: &[[i32; 2]],
    p2q1: &[[i32; 2]],
    invert_q: bool,
) -> (Vec<i32>, Vec<i32>) {
    let num_tri_p = in_p.halfedge.len() / 3;
    let num_tri_q = in_q.halfedge.len() / 3;
    let mut sides_per_face_pq = vec![0i32; num_tri_p + num_tri_q];

    for i in 0..in_p.halfedge.len() {
        sides_per_face_pq[i / 3] += i03[in_p.halfedge[i].start_vert as usize].abs();
    }
    for i in 0..in_q.halfedge.len() {
        sides_per_face_pq[num_tri_p + i / 3] += i30[in_q.halfedge[i].start_vert as usize].abs();
    }

    for i in 0..i12.len() {
        let edge_p = p1q2[i][0];
        let face_q = p1q2[i][1];
        let inclusion = i12[i].abs();
        sides_per_face_pq[face_q as usize + num_tri_p] += inclusion;
        sides_per_face_pq[edge_p as usize / 3] += inclusion;
        let paired = in_p.halfedge[edge_p as usize].paired_halfedge;
        if paired >= 0 {
            sides_per_face_pq[paired as usize / 3] += inclusion;
        }
    }

    for i in 0..i21.len() {
        let edge_q = p2q1[i][0];
        let face_p = p2q1[i][1];
        let inclusion = i21[i].abs();
        sides_per_face_pq[face_p as usize] += inclusion;
        sides_per_face_pq[num_tri_p + edge_q as usize / 3] += inclusion;
        let paired = in_q.halfedge[edge_q as usize].paired_halfedge;
        if paired >= 0 {
            sides_per_face_pq[num_tri_p + paired as usize / 3] += inclusion;
        }
    }

    let mut face_pq2r = vec![0i32; num_tri_p + num_tri_q + 1];
    let mut curr_face_r = 0;
    for i in 0..sides_per_face_pq.len() {
        face_pq2r[i + 1] = curr_face_r;
        if sides_per_face_pq[i] > 0 {
            curr_face_r += 1;
        }
    }
    let num_face_r = curr_face_r;
    face_pq2r.resize(num_tri_p + num_tri_q, 0);

    out_r.face_normal = Vec::with_capacity(num_face_r as usize);
    for i in 0..num_tri_p {
        if sides_per_face_pq[i] > 0 {
            out_r.face_normal.push(in_p.face_normal[i]);
        }
    }
    for i in 0..num_tri_q {
        if sides_per_face_pq[i + num_tri_p] > 0 {
            let normal = in_q.face_normal[i];
            out_r
                .face_normal
                .push(if invert_q { -normal } else { normal });
        }
    }

    let mut face_edge = vec![0i32; num_face_r as usize + 1];
    let mut curr_edge = 0;
    let mut j = 0;
    for i in 0..sides_per_face_pq.len() {
        if sides_per_face_pq[i] > 0 {
            face_edge[j + 1] = curr_edge;
            curr_edge += sides_per_face_pq[i];
            j += 1;
        }
    }
    out_r.halfedge = vec![Halfedge::default(); curr_edge as usize];

    (face_edge, face_pq2r)
}

pub fn append_whole_edges(
    mesh_r: &mut ManifoldImpl,
    face_ptr: &mut [i32],
    halfedge_ref: &mut [TriRef],
    in_p: &ManifoldImpl,
    whole_halfedge_p: &[bool],
    i03: &[i32],
    v_p2r: &[i32],
    face_p2r: &[i32],
    forward: bool,
) {
    for idx in 0..in_p.halfedge.len() {
        if !whole_halfedge_p[idx] {
            continue;
        }
        let mut halfedge = in_p.halfedge[idx];
        if !halfedge.is_forward() {
            continue;
        }

        let inclusion = i03[halfedge.start_vert as usize];
        if inclusion == 0 {
            continue;
        }
        if inclusion < 0 {
            std::mem::swap(&mut halfedge.start_vert, &mut halfedge.end_vert);
        }
        halfedge.start_vert = v_p2r[halfedge.start_vert as usize];
        halfedge.end_vert = v_p2r[halfedge.end_vert as usize];

        let face_left_p = idx / 3;
        let new_face = face_p2r[face_left_p];
        let face_right_p = halfedge.paired_halfedge / 3;
        let face_right = face_p2r[face_right_p as usize];

        let forward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_left_p as i32,
            coplanar_id: -1,
        };
        let backward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_right_p,
            coplanar_id: -1,
        };

        for _ in 0..inclusion.abs() {
            let forward_edge = face_ptr[new_face as usize];
            face_ptr[new_face as usize] += 1;
            let backward_edge = face_ptr[face_right as usize];
            face_ptr[face_right as usize] += 1;

            halfedge.paired_halfedge = backward_edge;
            mesh_r.halfedge[forward_edge as usize] = halfedge;
            mesh_r.halfedge[backward_edge as usize] = Halfedge {
                start_vert: halfedge.end_vert,
                end_vert: halfedge.start_vert,
                paired_halfedge: forward_edge,
                prop_vert: halfedge.end_vert,
            };

            halfedge_ref[forward_edge as usize] = forward_ref;
            halfedge_ref[backward_edge as usize] = backward_ref;

            halfedge.start_vert += 1;
            halfedge.end_vert += 1;
        }
    }
}

pub fn append_partial_edges(
    mesh_r: &mut ManifoldImpl,
    whole_halfedge_p: &mut [bool],
    face_ptr_r: &mut [i32],
    edges_p: &mut BTreeMap<i32, Vec<EdgePos>>,
    halfedge_ref: &mut [TriRef],
    in_p: &ManifoldImpl,
    i03: &[i32],
    v_p2r: &[i32],
    face_p2r: &[i32],
    forward: bool,
) {
    let vert_pos_p = &in_p.vert_pos;
    let halfedge_p = &in_p.halfedge;

    for (&edge_p_idx, edge_pos_p_vec) in edges_p.iter_mut() {
        edge_pos_p_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let halfedge = halfedge_p[edge_p_idx as usize];
        whole_halfedge_p[edge_p_idx as usize] = false;
        whole_halfedge_p[halfedge.paired_halfedge as usize] = false;

        let v_start = halfedge.start_vert;
        let v_end = halfedge.end_vert;
        let edge_vec = vert_pos_p[v_end as usize] - vert_pos_p[v_start as usize];

        for edge in edge_pos_p_vec.iter_mut() {
            edge.edge_pos = mesh_r.vert_pos[edge.vert as usize].dot(edge_vec);
        }

        let inclusion_start = i03[v_start as usize];
        let mut edge_pos_obj = EdgePos {
            edge_pos: mesh_r.vert_pos[v_p2r[v_start as usize] as usize].dot(edge_vec),
            vert: v_p2r[v_start as usize],
            collision_id: i32::MAX,
            is_start: inclusion_start > 0,
        };
        for _ in 0..inclusion_start.abs() {
            edge_pos_p_vec.push(edge_pos_obj);
            edge_pos_obj.vert += 1;
        }

        let inclusion_end = i03[v_end as usize];
        edge_pos_obj = EdgePos {
            edge_pos: mesh_r.vert_pos[v_p2r[v_end as usize] as usize].dot(edge_vec),
            vert: v_p2r[v_end as usize],
            collision_id: i32::MAX,
            is_start: inclusion_end < 0,
        };
        for _ in 0..inclusion_end.abs() {
            edge_pos_p_vec.push(edge_pos_obj);
            edge_pos_obj.vert += 1;
        }

        let edges = pair_up(edge_pos_p_vec);

        let face_left_p = edge_p_idx / 3;
        let face_left = face_p2r[face_left_p as usize];
        let face_right_p = halfedge.paired_halfedge / 3;
        let face_right = face_p2r[face_right_p as usize];

        let forward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_left_p as i32,
            coplanar_id: -1,
        };
        let backward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_right_p as i32,
            coplanar_id: -1,
        };

        for mut e in edges {
            let forward_edge = face_ptr_r[face_left as usize];
            face_ptr_r[face_left as usize] += 1;
            let backward_edge = face_ptr_r[face_right as usize];
            face_ptr_r[face_right as usize] += 1;

            e.paired_halfedge = backward_edge;
            mesh_r.halfedge[forward_edge as usize] = e;
            halfedge_ref[forward_edge as usize] = forward_ref;

            std::mem::swap(&mut e.start_vert, &mut e.end_vert);
            e.paired_halfedge = forward_edge;
            mesh_r.halfedge[backward_edge as usize] = e;
            halfedge_ref[backward_edge as usize] = backward_ref;
        }
    }
}

pub fn append_new_edges(
    mesh_r: &mut ManifoldImpl,
    face_ptr_r: &mut [i32],
    edges_new: &mut BTreeMap<(i32, i32), Vec<EdgePos>>,
    halfedge_ref: &mut [TriRef],
    face_pq2r: &[i32],
    num_face_p: usize,
) {
    for (&(face_p, face_q), edge_pos_vec) in edges_new.iter_mut() {
        edge_pos_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut bbox = GeoBox::default();
        for edge in edge_pos_vec.iter() {
            bbox.union_point(mesh_r.vert_pos[edge.vert as usize]);
        }
        let size = bbox.max - bbox.min;
        let axis = if size.x > size.y && size.x > size.z {
            0
        } else if size.y > size.z {
            1
        } else {
            2
        };

        for edge in edge_pos_vec.iter_mut() {
            edge.edge_pos = mesh_r.vert_pos[edge.vert as usize].to_array()[axis];
        }

        let edges = pair_up(edge_pos_vec);

        let face_left = face_pq2r[face_p as usize];
        let face_right = face_pq2r[num_face_p + face_q as usize];
        let forward_ref = TriRef {
            mesh_id: 0,
            original_id: -1,
            face_id: face_p,
            coplanar_id: -1,
        };
        let backward_ref = TriRef {
            mesh_id: 1,
            original_id: -1,
            face_id: face_q,
            coplanar_id: -1,
        };

        for mut e in edges {
            let forward_edge = face_ptr_r[face_left as usize];
            face_ptr_r[face_left as usize] += 1;
            let backward_edge = face_ptr_r[face_right as usize];
            face_ptr_r[face_right as usize] += 1;

            e.paired_halfedge = backward_edge;
            mesh_r.halfedge[forward_edge as usize] = e;
            halfedge_ref[forward_edge as usize] = forward_ref;

            std::mem::swap(&mut e.start_vert, &mut e.end_vert);
            e.paired_halfedge = forward_edge;
            mesh_r.halfedge[backward_edge as usize] = e;
            halfedge_ref[backward_edge as usize] = backward_ref;
        }
    }
}

pub fn create_boolean_result(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    xv12: &Intersections,
    xv21: &Intersections,
    w03: &Vec<i32>,
    w30: &Vec<i32>,
    op: OpType,
) -> BooleanResult {
    let c1 = if matches!(op, OpType::Intersect) {
        0
    } else {
        1
    };
    let c2 = if matches!(op, OpType::Add) { 1 } else { 0 };
    let c3 = if matches!(op, OpType::Intersect) {
        1
    } else {
        -1
    };

    let i12: Vec<i32> = xv12.x12.iter().map(|&v| c3 * v).collect();
    let i21: Vec<i32> = xv21.x12.iter().map(|&v| c3 * v).collect();
    let i03: Vec<i32> = w03.iter().map(|&v| c1 + c3 * v).collect();
    let i30: Vec<i32> = w30.iter().map(|&v| c2 + c3 * v).collect();

    let mut v_p2r = vec![0i32; in_p.vert_pos.len()];
    let mut curr_offset = 0;
    for i in 0..in_p.vert_pos.len() {
        v_p2r[i] = curr_offset;
        curr_offset += i03[i].abs();
    }
    let n_pv = curr_offset;

    let mut v_q2r = vec![0i32; in_q.vert_pos.len()];
    for i in 0..in_q.vert_pos.len() {
        v_q2r[i] = curr_offset;
        curr_offset += i30[i].abs();
    }
    let n_qv = curr_offset - n_pv;

    let mut v12_r = vec![0i32; xv12.v12.len()];
    for i in 0..xv12.v12.len() {
        v12_r[i] = curr_offset;
        curr_offset += i12[i].abs();
    }
    let n12 = curr_offset - n_pv - n_qv;

    let mut v21_r = vec![0i32; xv21.v12.len()];
    for i in 0..xv21.v12.len() {
        v21_r[i] = curr_offset;
        curr_offset += i21[i].abs();
    }

    let mut mesh_r = ManifoldImpl::new();
    let num_vert_r = curr_offset as usize;
    if num_vert_r == 0 {
        return BooleanResult { mesh_r };
    }

    mesh_r.vert_pos.resize(num_vert_r, DVec3::ZERO);
    for (i, &v) in in_p.vert_pos.iter().enumerate() {
        for j in 0..i03[i].abs() {
            mesh_r.vert_pos[(v_p2r[i] + j) as usize] = v;
        }
    }
    for (i, &v) in in_q.vert_pos.iter().enumerate() {
        for j in 0..i30[i].abs() {
            mesh_r.vert_pos[(v_q2r[i] + j) as usize] = v;
        }
    }
    for (i, &v) in xv12.v12.iter().enumerate() {
        for j in 0..i12[i].abs() {
            mesh_r.vert_pos[(v12_r[i] + j) as usize] = v;
        }
    }
    for (i, &v) in xv21.v12.iter().enumerate() {
        for j in 0..i21[i].abs() {
            mesh_r.vert_pos[(v21_r[i] + j) as usize] = v;
        }
    }

    let mut edges_p: BTreeMap<i32, Vec<EdgePos>> = BTreeMap::new();
    let mut edges_q: BTreeMap<i32, Vec<EdgePos>> = BTreeMap::new();
    let mut edges_new: BTreeMap<(i32, i32), Vec<EdgePos>> = BTreeMap::new();

    add_new_edge_verts(
        &mut edges_p,
        &mut edges_new,
        &xv12.p1q2,
        &i12,
        &v12_r,
        &in_p.halfedge,
        true,
        0,
    );
    add_new_edge_verts(
        &mut edges_q,
        &mut edges_new,
        &xv21.p1q2,
        &i21,
        &v21_r,
        &in_q.halfedge,
        false,
        xv12.v12.len(),
    );

    let (face_edge, face_pq2r) = size_output(
        &mut mesh_r,
        in_p,
        in_q,
        &i03,
        &i30,
        &i12,
        &i21,
        &xv12.p1q2,
        &xv21.p1q2,
        matches!(op, OpType::Subtract),
    );

    let num_tri_p = in_p.halfedge.len() / 3;
    let mut face_ptr_r = face_edge.clone();
    let mut whole_halfedge_p = vec![true; in_p.halfedge.len()];
    let mut whole_halfedge_q = vec![true; in_q.halfedge.len()];
    let mut halfedge_ref = vec![TriRef::default(); mesh_r.halfedge.len()];

    append_partial_edges(
        &mut mesh_r,
        &mut whole_halfedge_p,
        &mut face_ptr_r,
        &mut edges_p,
        &mut halfedge_ref,
        in_p,
        &i03,
        &v_p2r,
        &face_pq2r,
        true,
    );
    append_partial_edges(
        &mut mesh_r,
        &mut whole_halfedge_q,
        &mut face_ptr_r,
        &mut edges_q,
        &mut halfedge_ref,
        in_q,
        &i30,
        &v_q2r,
        &face_pq2r[num_tri_p..],
        false,
    );

    append_new_edges(
        &mut mesh_r,
        &mut face_ptr_r,
        &mut edges_new,
        &mut halfedge_ref,
        &face_pq2r,
        num_tri_p,
    );

    append_whole_edges(
        &mut mesh_r,
        &mut face_ptr_r,
        &mut halfedge_ref,
        in_p,
        &whole_halfedge_p,
        &i03,
        &v_p2r,
        &face_pq2r,
        true,
    );
    append_whole_edges(
        &mut mesh_r,
        &mut face_ptr_r,
        &mut halfedge_ref,
        in_q,
        &whole_halfedge_q,
        &i30,
        &v_q2r,
        &face_pq2r[num_tri_p..],
        false,
    );

    BooleanResult { mesh_r }
}
