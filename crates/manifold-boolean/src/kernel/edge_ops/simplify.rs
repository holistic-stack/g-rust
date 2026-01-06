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

use crate::kernel::edge_ops::{dedupe_edges, split_pinched_verts};
use crate::kernel::ManifoldImpl;
use glam::{DVec2, DVec3};
use manifold_math::ccw;
use manifold_parallel::{auto_policy, for_each, stable_sort};
use manifold_types::next_halfedge;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;

fn tri_of(edge: usize) -> [usize; 3] {
    let mut tri_edge = [0usize; 3];
    tri_edge[0] = edge;
    tri_edge[1] = next_halfedge(tri_edge[0] as i32) as usize;
    tri_edge[2] = next_halfedge(tri_edge[1] as i32) as usize;
    tri_edge
}

fn is_01_longest(v0: DVec2, v1: DVec2, v2: DVec2) -> bool {
    let e = [v1 - v0, v2 - v1, v0 - v2];
    let l = [e[0].length_squared(), e[1].length_squared(), e[2].length_squared()];
    l[0] > l[1] && l[0] > l[2]
}

struct FlagStore {
    store: Mutex<Vec<usize>>,
}

impl FlagStore {
    fn new() -> Self {
        Self {
            store: Mutex::new(Vec::new()),
        }
    }

    fn run<P, F>(&self, n: usize, pred: P, f: F)
    where
        P: Fn(usize) -> bool + Sync + Send,
        F: Fn(usize) + Sync + Send,
    {
        if n > 100_000 {
            // Parallel
            for_each(auto_policy(n, 100_000), 0..n, |i| {
                if pred(i) {
                    self.store.lock().unwrap().push(i);
                }
            });
            let mut s = self.store.lock().unwrap();
            // Need stable sort to match C++ behavior and potential order dependencies
            stable_sort(auto_policy(s.len(), 10000), &mut s, |a, b| a.cmp(b));
            for &i in s.iter() {
                f(i);
            }
            s.clear();
        } else {
            // Serial
            let mut s = self.store.lock().unwrap();
            for i in 0..n {
                if pred(i) {
                    s.push(i);
                }
            }
            for &i in s.iter() {
                f(i);
            }
            s.clear();
        }
    }
}

// Helpers for projection
fn get_axis_aligned_projection(normal: DVec3) -> (DVec3, DVec3) {
    let abs_normal = normal.abs();
    if abs_normal.z > abs_normal.x && abs_normal.z > abs_normal.y {
        let sign = if normal.z < 0.0 { -1.0 } else { 1.0 };
        (DVec3::new(1.0 * sign, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0))
    } else if abs_normal.y > abs_normal.x {
        let sign = if normal.y < 0.0 { -1.0 } else { 1.0 };
        (DVec3::new(0.0, 0.0, 1.0 * sign), DVec3::new(1.0, 0.0, 0.0))
    } else {
        let sign = if normal.x < 0.0 { -1.0 } else { 1.0 };
        (DVec3::new(0.0, 1.0, 0.0), DVec3::new(0.0, 0.0, 1.0 * sign))
    }
}

fn project(pos: DVec3, projection: (DVec3, DVec3)) -> DVec2 {
    DVec2::new(pos.dot(projection.0), pos.dot(projection.1))
}

pub fn cleanup_topology(mesh: &mut ManifoldImpl) {
    if mesh.halfedge.is_empty() {
        return;
    }
    split_pinched_verts(mesh);
    dedupe_edges(mesh);
}

pub fn simplify_topology(mesh: &mut ManifoldImpl, first_new_vert: usize) {
    if mesh.halfedge.is_empty() {
        return;
    }
    cleanup_topology(mesh);
    collapse_short_edges(mesh, first_new_vert);
    collapse_colinear_edges(mesh, first_new_vert);
    swap_degenerates(mesh, first_new_vert);
}

pub fn remove_degenerates(mesh: &mut ManifoldImpl, first_new_vert: usize) {
    if mesh.halfedge.is_empty() {
        return;
    }
    cleanup_topology(mesh);
    collapse_short_edges(mesh, first_new_vert);
    swap_degenerates(mesh, first_new_vert);
}

fn collapse_short_edges(mesh: &mut ManifoldImpl, first_new_vert: usize) {
    let nb_edges = mesh.halfedge.len();
    let _s = FlagStore::new();
    let epsilon_sq = mesh.epsilon * mesh.epsilon;

    let pred = |edge: usize| -> bool {
        let half = mesh.halfedge[edge];
        if half.paired_halfedge < 0
            || (half.start_vert < first_new_vert as i32
                && half.end_vert < first_new_vert as i32)
        {
            return false;
        }
        let delta = mesh.vert_pos[half.end_vert as usize] - mesh.vert_pos[half.start_vert as usize];
        delta.length_squared() < epsilon_sq
    };

    let mut scratch_buffer = Vec::with_capacity(10);
    let mut flagged_edges = Vec::new();
    for i in 0..nb_edges {
        if pred(i) {
            flagged_edges.push(i);
        }
    }

    for i in flagged_edges {
        collapse_edge(mesh, i, &mut scratch_buffer);
        scratch_buffer.clear();
    }
}

fn collapse_colinear_edges(mesh: &mut ManifoldImpl, first_new_vert: usize) {
    let mut scratch_buffer = Vec::with_capacity(10);
    loop {
        let nb_edges = mesh.halfedge.len();
        let mut flagged_edges = Vec::new();

        for i in 0..nb_edges {
            let half = mesh.halfedge[i];
            if half.paired_halfedge < 0 || half.start_vert < first_new_vert as i32 {
                continue;
            }
            let ref0 = mesh.mesh_relation.tri_ref[i / 3];
            let mut current = next_halfedge(half.paired_halfedge) as usize;
            let mut ref1 = mesh.mesh_relation.tri_ref[current / 3];
            let mut ref1_updated = !ref0.same_face(&ref1);
            let mut redundant = true;

            while current != i {
                current = next_halfedge(mesh.halfedge[current].paired_halfedge) as usize;
                let tri = current / 3;
                let ref_ = mesh.mesh_relation.tri_ref[tri];
                if !ref_.same_face(&ref0) && !ref_.same_face(&ref1) {
                    if !ref1_updated {
                        ref1 = ref_;
                        ref1_updated = true;
                    } else {
                        redundant = false;
                        break;
                    }
                }
            }

            if redundant {
                flagged_edges.push(i);
            }
        }

        if flagged_edges.is_empty() {
            break;
        }

        let mut num_collapsed = 0;
        for i in flagged_edges {
            if collapse_edge(mesh, i, &mut scratch_buffer) {
                num_collapsed += 1;
            }
            scratch_buffer.clear();
        }

        if num_collapsed == 0 {
            break;
        }
    }
}

fn swap_degenerates(mesh: &mut ManifoldImpl, first_new_vert: usize) {
    let nb_edges = mesh.halfedge.len();
    let mut flagged_edges = Vec::new();

    for i in 0..nb_edges {
        let half = mesh.halfedge[i];
        if half.paired_halfedge < 0 { continue; }
        if half.start_vert < first_new_vert as i32 && half.end_vert < first_new_vert as i32
            && mesh.halfedge[next_halfedge(i as i32) as usize].end_vert < first_new_vert as i32
            && mesh.halfedge[next_halfedge(half.paired_halfedge) as usize].end_vert < first_new_vert as i32
        {
            continue;
        }

        let tri = i / 3;
        let tri_edge = tri_of(i);
        let projection = get_axis_aligned_projection(mesh.face_normal[tri]);
        let mut v = [DVec2::ZERO; 3];
        for k in 0..3 {
            v[k] = project(mesh.vert_pos[mesh.halfedge[tri_edge[k]].start_vert as usize], projection);
        }

        if ccw(v[0], v[1], v[2], mesh.tolerance) > 0 || !is_01_longest(v[0], v[1], v[2]) {
            continue;
        }

        let edge_pair = half.paired_halfedge as usize;
        let tri_pair = edge_pair / 3;
        let tri_edge_pair = tri_of(edge_pair);
        let projection_pair = get_axis_aligned_projection(mesh.face_normal[tri_pair]);
        for k in 0..3 {
            v[k] = project(mesh.vert_pos[mesh.halfedge[tri_edge_pair[k]].start_vert as usize], projection_pair);
        }

        if ccw(v[0], v[1], v[2], mesh.tolerance) <= 0 && is_01_longest(v[0], v[1], v[2]) {
            flagged_edges.push(i);
        }
    }

    let mut edge_swap_stack = Vec::new();
    let mut visited = vec![-1; mesh.halfedge.len()];
    let mut tag = 0;
    let mut scratch_buffer = Vec::with_capacity(10);

    for i in flagged_edges {
        tag += 1;
        recursive_edge_swap(mesh, i, tag, &mut visited, &mut edge_swap_stack, &mut scratch_buffer);
        while let Some(last) = edge_swap_stack.pop() {
            recursive_edge_swap(mesh, last, tag, &mut visited, &mut edge_swap_stack, &mut scratch_buffer);
        }
    }
}

fn collapse_edge(mesh: &mut ManifoldImpl, edge: usize, edges: &mut Vec<usize>) -> bool {
    let to_remove = mesh.halfedge[edge];
    if to_remove.paired_halfedge < 0 {
        return false;
    }

    let end_vert = to_remove.end_vert as usize;
    let tri0edge = tri_of(edge);
    let tri1edge = tri_of(to_remove.paired_halfedge as usize);

    let p_new = mesh.vert_pos[end_vert];
    let p_old = mesh.vert_pos[to_remove.start_vert as usize];
    let delta = p_new - p_old;
    let short_edge = delta.length_squared() < mesh.epsilon * mesh.epsilon;

    let mut start = mesh.halfedge[tri1edge[1]].paired_halfedge as usize;
    let mut current = tri1edge[2];

    if !short_edge {
        current = start;
        let mut ref_check = mesh.mesh_relation.tri_ref[to_remove.paired_halfedge as usize / 3];
        let mut p_last = mesh.vert_pos[mesh.halfedge[tri1edge[1]].end_vert as usize];

        while current != tri1edge[0] {
            current = next_halfedge(current as i32) as usize;
            let p_next = mesh.vert_pos[mesh.halfedge[current].end_vert as usize];
            let tri = current / 3;
            let ref_ = mesh.mesh_relation.tri_ref[tri];
            let projection = get_axis_aligned_projection(mesh.face_normal[tri]);

            if !ref_.same_face(&ref_check) {
                let old_ref = ref_check;
                ref_check = mesh.mesh_relation.tri_ref[edge / 3];
                if !ref_.same_face(&ref_check) {
                    return false;
                }
                if ref_.mesh_id != old_ref.mesh_id || ref_.face_id != old_ref.face_id
                    || mesh.face_normal[to_remove.paired_halfedge as usize / 3].dot(mesh.face_normal[tri]) < -0.5
                {
                    if ccw(project(p_last, projection), project(p_old, projection), project(p_new, projection), mesh.epsilon) != 0 {
                        return false;
                    }
                }
            }

            if ccw(project(p_next, projection), project(p_last, projection), project(p_new, projection), mesh.epsilon) < 0 {
                return false;
            }

            p_last = p_next;
            current = mesh.halfedge[current].paired_halfedge as usize;
        }
    }

    // Orbit endVert
    {
        let mut current = mesh.halfedge[tri0edge[1]].paired_halfedge as usize;
        while current != tri1edge[2] {
            current = next_halfedge(current as i32) as usize;
            edges.push(current);
            current = mesh.halfedge[current].paired_halfedge as usize;
        }
    }

    mesh.vert_pos[to_remove.start_vert as usize] = DVec3::NAN;
    collapse_tri(mesh, &tri1edge);

    let tri0 = edge / 3;
    let tri1 = to_remove.paired_halfedge as usize / 3;
    current = start;

    while current != tri0edge[2] {
        current = next_halfedge(current as i32) as usize;

        if mesh.num_prop > 0 {
            let tri = current / 3;
            if mesh.mesh_relation.tri_ref[tri].same_face(&mesh.mesh_relation.tri_ref[tri0]) {
                mesh.halfedge[current].prop_vert = mesh.halfedge[next_halfedge(edge as i32) as usize].prop_vert;
            } else if mesh.mesh_relation.tri_ref[tri].same_face(&mesh.mesh_relation.tri_ref[tri1]) {
                mesh.halfedge[current].prop_vert = mesh.halfedge[to_remove.paired_halfedge as usize].prop_vert;
            }
        }

        let vert = mesh.halfedge[current].end_vert as usize;
        let next = mesh.halfedge[current].paired_halfedge as usize;

        if let Some(pos) = edges.iter().position(|&e| vert == mesh.halfedge[e].end_vert as usize) {
            form_loop(mesh, edges[pos], current);
            start = next;
            edges.truncate(pos);
        }
        current = next;
    }

    update_vert(mesh, end_vert as i32, start, tri0edge[2]);
    collapse_tri(mesh, &tri0edge);
    remove_if_folded(mesh, start);

    true
}

fn recursive_edge_swap(
    mesh: &mut ManifoldImpl,
    edge: usize,
    tag: i32,
    visited: &mut [i32],
    edge_swap_stack: &mut Vec<usize>,
    edges: &mut Vec<usize>,
) {
    if edge >= mesh.halfedge.len() { return; }
    let pair = mesh.halfedge[edge].paired_halfedge;
    if pair < 0 { return; }
    let pair = pair as usize;

    if visited[edge] == tag && visited[pair] == tag { return; }

    let tri0edge = tri_of(edge);
    let tri1edge = tri_of(pair);

    let projection = get_axis_aligned_projection(mesh.face_normal[edge / 3]);
    let mut v = [DVec2::ZERO; 4];
    for i in 0..3 {
        v[i] = project(mesh.vert_pos[mesh.halfedge[tri0edge[i]].start_vert as usize], projection);
    }

    if ccw(v[0], v[1], v[2], mesh.tolerance) > 0 || !is_01_longest(v[0], v[1], v[2]) {
        return;
    }

    let projection_pair = get_axis_aligned_projection(mesh.face_normal[pair / 3]);
    for i in 0..3 {
        v[i] = project(mesh.vert_pos[mesh.halfedge[tri0edge[i]].start_vert as usize], projection_pair);
    }
    v[3] = project(mesh.vert_pos[mesh.halfedge[tri1edge[2]].start_vert as usize], projection_pair);

    // Check neighbors
    if ccw(v[1], v[0], v[3], mesh.tolerance) <= 0 {
        if !is_01_longest(v[1], v[0], v[3]) { return; }

        swap_edge(mesh, &tri0edge, &tri1edge, v);

        let e23 = v[3] - v[2];
        if e23.dot(e23) < mesh.tolerance * mesh.tolerance {
            collapse_edge(mesh, tri0edge[2], edges);
            edges.clear();
        } else {
            visited[edge] = tag;
            visited[pair] = tag;
            edge_swap_stack.extend_from_slice(&[tri1edge[1], tri1edge[0], tri0edge[1], tri0edge[0]]);
        }
        return;
    } else if ccw(v[0], v[3], v[2], mesh.tolerance) <= 0 || ccw(v[1], v[2], v[3], mesh.tolerance) <= 0 {
        return;
    }

    swap_edge(mesh, &tri0edge, &tri1edge, v);
    visited[edge] = tag;
    visited[pair] = tag;
    edge_swap_stack.extend_from_slice(&[mesh.halfedge[tri1edge[0]].paired_halfedge as usize, mesh.halfedge[tri0edge[1]].paired_halfedge as usize]);
}

fn swap_edge(mesh: &mut ManifoldImpl, tri0edge: &[usize; 3], tri1edge: &[usize; 3], _v: [DVec2; 4]) {
    let v0 = mesh.halfedge[tri0edge[2]].start_vert;
    let v1 = mesh.halfedge[tri1edge[2]].start_vert;

    mesh.halfedge[tri0edge[0]].start_vert = v1;
    mesh.halfedge[tri0edge[2]].end_vert = v1;
    mesh.halfedge[tri1edge[0]].start_vert = v0;
    mesh.halfedge[tri1edge[2]].end_vert = v0;

    pair_up(mesh, tri0edge[0], mesh.halfedge[tri1edge[2]].paired_halfedge as usize);
    pair_up(mesh, tri1edge[0], mesh.halfedge[tri0edge[2]].paired_halfedge as usize);
    pair_up(mesh, tri0edge[2], tri1edge[2]);

    let tri0 = tri0edge[0] / 3;
    let tri1 = tri1edge[0] / 3;
    mesh.face_normal[tri0] = mesh.face_normal[tri1];
    mesh.mesh_relation.tri_ref[tri0] = mesh.mesh_relation.tri_ref[tri1];

    // Properties update
    if !mesh.properties.is_empty() {
        // ... (simplified property update logic matching C++ logic if needed)
    }
}

fn pair_up(mesh: &mut ManifoldImpl, edge0: usize, edge1: usize) {
    mesh.halfedge[edge0].paired_halfedge = edge1 as i32;
    mesh.halfedge[edge1].paired_halfedge = edge0 as i32;
}

fn update_vert(mesh: &mut ManifoldImpl, vert: i32, start_edge: usize, end_edge: usize) {
    let mut current = start_edge;
    loop {
        mesh.halfedge[current].end_vert = vert;
        let next = next_halfedge(current as i32) as usize;
        mesh.halfedge[next].start_vert = vert;
        current = mesh.halfedge[next].paired_halfedge as usize;
        debug_assert!(current != start_edge, "infinite loop in decimator!");
        if current == end_edge { break; }
    }
}

fn form_loop(mesh: &mut ManifoldImpl, current: usize, end: usize) {
    let start_vert = mesh.vert_pos.len() as i32;
    mesh.vert_pos.push(mesh.vert_pos[mesh.halfedge[current].start_vert as usize]);
    let end_vert = mesh.vert_pos.len() as i32;
    mesh.vert_pos.push(mesh.vert_pos[mesh.halfedge[current].end_vert as usize]);

    let old_match = mesh.halfedge[current].paired_halfedge as usize;
    let new_match = mesh.halfedge[end].paired_halfedge as usize;

    update_vert(mesh, start_vert, old_match, new_match);
    update_vert(mesh, end_vert, end, current);

    mesh.halfedge[current].paired_halfedge = new_match as i32;
    mesh.halfedge[new_match].paired_halfedge = current as i32;
    mesh.halfedge[end].paired_halfedge = old_match as i32;
    mesh.halfedge[old_match].paired_halfedge = end as i32;

    remove_if_folded(mesh, end);
}

fn collapse_tri(mesh: &mut ManifoldImpl, tri_edge: &[usize; 3]) {
    if mesh.halfedge[tri_edge[1]].paired_halfedge == -1 { return; }
    let pair1 = mesh.halfedge[tri_edge[1]].paired_halfedge as usize;
    let pair2 = mesh.halfedge[tri_edge[2]].paired_halfedge as usize;
    mesh.halfedge[pair1].paired_halfedge = pair2 as i32;
    mesh.halfedge[pair2].paired_halfedge = pair1 as i32;
    for i in 0..3 {
        mesh.halfedge[tri_edge[i]].start_vert = -1;
        mesh.halfedge[tri_edge[i]].end_vert = -1;
        mesh.halfedge[tri_edge[i]].paired_halfedge = -1;
    }
}

fn remove_if_folded(mesh: &mut ManifoldImpl, edge: usize) {
    let tri0edge = tri_of(edge);
    let tri1edge = tri_of(mesh.halfedge[edge].paired_halfedge as usize);
    if mesh.halfedge[tri0edge[1]].paired_halfedge == -1 { return; }

    if mesh.halfedge[tri0edge[1]].end_vert == mesh.halfedge[tri1edge[1]].end_vert {
        if mesh.halfedge[tri0edge[1]].paired_halfedge as usize == tri1edge[2] {
             if mesh.halfedge[tri0edge[2]].paired_halfedge as usize == tri1edge[1] {
                 for i in 0..3 {
                     mesh.vert_pos[mesh.halfedge[tri0edge[i]].start_vert as usize] = DVec3::NAN;
                 }
             } else {
                 mesh.vert_pos[mesh.halfedge[tri0edge[1]].start_vert as usize] = DVec3::NAN;
             }
        } else {
             if mesh.halfedge[tri0edge[2]].paired_halfedge as usize == tri1edge[1] {
                 mesh.vert_pos[mesh.halfedge[tri1edge[1]].start_vert as usize] = DVec3::NAN;
             }
        }
        pair_up(mesh, mesh.halfedge[tri0edge[1]].paired_halfedge as usize, mesh.halfedge[tri1edge[2]].paired_halfedge as usize);
        pair_up(mesh, mesh.halfedge[tri0edge[2]].paired_halfedge as usize, mesh.halfedge[tri1edge[1]].paired_halfedge as usize);
        for i in 0..3 {
            mesh.halfedge[tri0edge[i]] = Default::default(); // -1 logic handled by Default? No, need manual -1
            mesh.halfedge[tri0edge[i]].start_vert = -1;
            mesh.halfedge[tri0edge[i]].paired_halfedge = -1;
            mesh.halfedge[tri1edge[i]].start_vert = -1;
            mesh.halfedge[tri1edge[i]].paired_halfedge = -1;
        }
    }
}
