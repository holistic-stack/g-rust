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

// Winding03 - Vertex winding number computation
// C++ Reference: `submodules/manifold/src/boolean3.cpp:427-528`

use crate::kernel::kernel02::{kernel02_false_true, kernel02_true_true};
use crate::kernel::{Intersections, ManifoldImpl};
use manifold_collider::{Query, Recorder};
use manifold_parallel::{auto_policy, for_each, for_each_n, ExecutionPolicy, K_SEQ_THRESHOLD};
use manifold_types::DisjointSets;
use std::collections::HashSet;
use std::sync::Mutex;

struct WindingRecorder<'a, F>
where
    F: Fn(i32, i32),
{
    callback: F,
    _phantom: std::marker::PhantomData<&'a ()>,
}

impl<'a, F> WindingRecorder<'a, F>
where
    F: Fn(i32, i32) + Sync + Send,
{
    fn new(callback: F) -> Self {
        Self {
            callback,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<'a, F> Recorder for WindingRecorder<'a, F>
where
    F: Fn(i32, i32) + Sync + Send,
{
    fn record(&self, query_idx: i32, leaf_idx: i32) {
        (self.callback)(query_idx, leaf_idx);
    }
}

struct VertQuery<'a> {
    vert_pos: &'a [glam::DVec3],
    verts: &'a [i32],
}

impl<'a> Query for VertQuery<'a> {
    fn get_bbox(&self, i: i32) -> manifold_math::Box {
        let p = self.vert_pos[self.verts[i as usize] as usize];
        manifold_math::Box::new(p, p)
    }
}

fn winding03_inner<const EXPAND_P: bool, const FORWARD: bool>(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    intersections: &Intersections,
) -> Vec<i32> {
    // a: 0 (vert), b: 2 (face)
    let a = if FORWARD { in_p } else { in_q };
    let b = if FORWARD { in_q } else { in_p };
    let p1q2 = &intersections.p1q2;
    let index = if FORWARD { 0 } else { 1 };

    let num_verts = a.vert_pos.len();
    let u_a = DisjointSets::new(num_verts as u32);

    let halfedge_policy = auto_policy(a.halfedge.len(), K_SEQ_THRESHOLD);
    for_each_n(halfedge_policy, 0, a.halfedge.len(), |edge| {
        let he = &a.halfedge[edge];
        if !he.is_forward() {
            return;
        }

        // check if the edge is broken
        let edge_i32 = edge as i32;
        let is_broken = p1q2
            .binary_search_by(|probe| probe[index].cmp(&edge_i32))
            .is_ok();

        if !is_broken {
            u_a.unite(he.start_vert as u32, he.end_vert as u32);
        }
    });

    // find components
    let mut components = HashSet::new();
    for v in 0..num_verts {
        components.insert(u_a.find(v as u32));
    }

    let mut verts: Vec<i32> = components.into_iter().map(|c| c as i32).collect();
    verts.sort();

    let mut w03 = vec![0i32; num_verts];
    let w03_mutexes: Vec<Mutex<i32>> = (0..num_verts).map(|_| Mutex::new(0)).collect();

    let recorder_callback = |i: i32, b_idx: i32| {
        let vert_idx = verts[i as usize];
        let (s02, z02) = if EXPAND_P {
            kernel02_true_true(vert_idx, b_idx, a, b)
        } else {
            kernel02_false_true(vert_idx, b_idx, a, b)
        };
        if z02.is_finite() {
            let mut val = w03_mutexes[vert_idx as usize].lock().unwrap();
            *val += s02 * if FORWARD { 1 } else { -1 };
        }
    };

    let recorder = WindingRecorder::new(recorder_callback);
    let query = VertQuery {
        vert_pos: &a.vert_pos,
        verts: &verts,
    };

    if let Some(ref collider) = b.collider {
        collider.collisions::<false, _, _>(&query, verts.len() as i32, &recorder, true);
    }

    // Copy results from mutexes to w03
    for i in 0..num_verts {
        w03[i] = *w03_mutexes[i].lock().unwrap();
    }

    // flood fill
    let w03_policy = auto_policy(num_verts, K_SEQ_THRESHOLD);
    let w03_ptr = w03.as_mut_ptr() as usize;
    for_each_n(w03_policy, 0, num_verts, |i| {
        let root = u_a.find(i as u32) as usize;
        if root != i {
            unsafe {
                let w03_mut = std::slice::from_raw_parts_mut(w03_ptr as *mut i32, num_verts);
                w03_mut[i] = w03_mut[root];
            }
        }
    });

    w03
}

pub fn winding03_true_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    intersections: &Intersections,
    expand_p: bool,
) -> Vec<i32> {
    if expand_p {
        winding03_inner::<true, true>(in_p, in_q, intersections)
    } else {
        winding03_inner::<false, true>(in_p, in_q, intersections)
    }
}

pub fn winding03_false_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    intersections: &Intersections,
    expand_p: bool,
) -> Vec<i32> {
    if expand_p {
        winding03_inner::<true, false>(in_p, in_q, intersections)
    } else {
        winding03_inner::<false, false>(in_p, in_q, intersections)
    }
}

pub fn winding03(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    intersections: &Intersections,
    expand_p: bool,
) -> Vec<i32> {
    winding03_true_true(in_p, in_q, intersections, expand_p)
}

pub fn winding03_false_true(
    in_p: &ManifoldImpl,
    _in_q: &ManifoldImpl,
    _intersections: &Intersections,
    _expand_p: bool,
) -> Vec<i32> {
    vec![0i32; in_p.vert_pos.len()]
}

pub fn winding03(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    intersections: &Intersections,
    expand_p: bool,
    _forward: bool,
) -> Vec<i32> {
    if expand_p {
        winding03_true_true(in_p, in_q, intersections, true)
    } else {
        winding03_true_true(in_p, in_q, intersections, false)
    }
}
