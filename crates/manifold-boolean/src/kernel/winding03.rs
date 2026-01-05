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

use crate::kernel::kernel02;
use crate::kernel::{Intersections, ManifoldImpl};
use manifold_collider::{Query, Recorder};
use manifold_math::Box as GeoBox;
use manifold_types::DisjointSets;
use std::collections::HashSet;
use std::sync::atomic::{AtomicI32, Ordering};

pub struct Winding03Local {
    pub intersections: Intersections,
}

impl Winding03Local {
    pub fn new() -> Self {
        Self {
            intersections: Intersections::new(),
        }
    }
}

struct WindingQuery<'a> {
    pub verts: &'a [u32],
    pub vert_pos: &'a [glam::DVec3],
}

impl<'a> Query for WindingQuery<'a> {
    fn get_bbox(&self, i: i32) -> GeoBox {
        let p = self.vert_pos[self.verts[i as usize] as usize];
        GeoBox::new(p, p)
    }
}

struct WindingRecorder<'a> {
    pub w03: &'a [AtomicI32],
    pub verts: &'a [u32],
    pub in_p: &'a ManifoldImpl,
    pub in_q: &'a ManifoldImpl,
    pub forward: bool,
    pub expand_p: bool,
}

impl<'a> Recorder for WindingRecorder<'a> {
    fn record(&self, query_idx: i32, leaf_idx: i32) {
        let vert_idx = self.verts[query_idx as usize];
        let face_idx = leaf_idx;

        let (s02, z02) = if self.forward {
            if self.expand_p {
                kernel02::kernel02_true_true(vert_idx as i32, face_idx, self.in_p, self.in_q)
            } else {
                kernel02::kernel02_false_true(vert_idx as i32, face_idx, self.in_p, self.in_q)
            }
        } else {
            if self.expand_p {
                kernel02::kernel02_true_false(vert_idx as i32, face_idx, self.in_q, self.in_p)
            } else {
                kernel02::kernel02_false_false(vert_idx as i32, face_idx, self.in_q, self.in_p)
            }
        };

        if z02.is_finite() {
            let val = s02 * if self.forward { 1 } else { -1 };
            self.w03[vert_idx as usize].fetch_add(val, Ordering::SeqCst);
        }
    }
}

pub fn winding03_true_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    xv12: &Intersections,
    expand_p: bool,
) -> Vec<i32> {
    let n_p = in_p.vert_pos.len();
    let dset = DisjointSets::new(n_p as u32);

    let p1q2 = &xv12.p1q2;
    for (i, he) in in_p.halfedge.iter().enumerate() {
        if !he.is_forward() {
            continue;
        }

        let mut broken = false;
        let edge_idx = i as i32;
        let mut low = 0;
        let mut high = p1q2.len();
        while low < high {
            let mid = low + (high - low) / 2;
            if p1q2[mid][0] < edge_idx {
                low = mid + 1;
            } else {
                high = mid;
            }
        }

        if low < p1q2.len() && p1q2[low][0] == edge_idx {
            broken = true;
        }

        if !broken {
            dset.unite(he.start_vert as u32, he.end_vert as u32);
        }
    }

    let mut components = HashSet::new();
    for v in 0..n_p {
        components.insert(dset.find(v as u32));
    }
    let verts: Vec<u32> = components.into_iter().collect();

    let w03_atomic: Vec<AtomicI32> = (0..n_p).map(|_| AtomicI32::new(0)).collect();

    if let Some(collider) = &in_q.collider {
        let query = WindingQuery {
            verts: &verts,
            vert_pos: &in_p.vert_pos,
        };
        let recorder = WindingRecorder {
            w03: &w03_atomic,
            verts: &verts,
            in_p,
            in_q,
            forward: true,
            expand_p,
        };
        collider.collisions::<false, _, _>(&query, verts.len() as i32, &recorder, true);
    }

    let mut w03 = vec![0i32; n_p];
    for i in 0..n_p {
        let root = dset.find(i as u32);
        w03[i] = w03_atomic[root as usize].load(Ordering::SeqCst);
    }

    w03
}

pub fn winding03_false_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    xv21: &Intersections,
    expand_p: bool,
) -> Vec<i32> {
    let n_q = in_q.vert_pos.len();
    let dset = DisjointSets::new(n_q as u32);

    let p2q1 = &xv21.p1q2;
    for (i, he) in in_q.halfedge.iter().enumerate() {
        if !he.is_forward() {
            continue;
        }

        let mut broken = false;
        let edge_idx = i as i32;
        let mut low = 0;
        let mut high = p2q1.len();
        while low < high {
            let mid = low + (high - low) / 2;
            if p2q1[mid][0] < edge_idx {
                low = mid + 1;
            } else {
                high = mid;
            }
        }

        if low < p2q1.len() && p2q1[low][0] == edge_idx {
            broken = true;
        }

        if !broken {
            dset.unite(he.start_vert as u32, he.end_vert as u32);
        }
    }

    let mut components = HashSet::new();
    for v in 0..n_q {
        components.insert(dset.find(v as u32));
    }
    let verts: Vec<u32> = components.into_iter().collect();

    let w30_atomic: Vec<AtomicI32> = (0..n_q).map(|_| AtomicI32::new(0)).collect();

    if let Some(collider) = &in_p.collider {
        let query = WindingQuery {
            verts: &verts,
            vert_pos: &in_q.vert_pos,
        };
        let recorder = WindingRecorder {
            w03: &w30_atomic,
            verts: &verts,
            in_p: in_p,
            in_q: in_q,
            forward: false,
            expand_p,
        };
        collider.collisions::<false, _, _>(&query, verts.len() as i32, &recorder, true);
    }

    let mut w30 = vec![0i32; n_q];
    for i in 0..n_q {
        let root = dset.find(i as u32);
        w30[i] = w30_atomic[root as usize].load(Ordering::SeqCst);
    }

    w30
}
