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

//! Intersect12 - Edge-Face intersection computation
//!
//! C++ Reference: `submodules/manifold/src/boolean3.cpp:378-424`

use crate::kernel::kernel12::{
    kernel12_false_false, kernel12_false_true, kernel12_true_false, kernel12_true_true,
};
use crate::kernel::{Intersections, ManifoldImpl};
use manifold_collider::{Query, Recorder};
use manifold_math::Box as GeoBox;
use std::sync::Mutex;

/// Query wrapper for edge bounding box queries
struct EdgeBBoxQuery {
    vert_pos: Vec<glam::DVec3>,
    halfedge: Vec<manifold_types::Halfedge>,
}

impl EdgeBBoxQuery {
    fn new(mesh: &ManifoldImpl) -> Self {
        Self {
            vert_pos: mesh.vert_pos.clone(),
            halfedge: mesh.halfedge.clone(),
        }
    }
}

impl Query for EdgeBBoxQuery {
    fn get_bbox(&self, halfedge_idx: i32) -> GeoBox {
        let he = &self.halfedge[halfedge_idx as usize];
        if he.is_forward() {
            let start = self.vert_pos[he.start_vert as usize];
            let end = self.vert_pos[he.end_vert as usize];
            GeoBox::new(start, end)
        } else {
            GeoBox::default()
        }
    }
}

/// Thread-safe intersection storage
///
/// DETERMINISM: Intersections may be added in nondeterministic order due to
/// multi-threaded collision detection. However, determinism is restored
/// by stable sorting in intersect12_inner before returning results.
struct IntersectionStore {
    p1q2: Vec<[i32; 2]>,
    x12: Vec<i32>,
    v12: Vec<glam::DVec3>,
}

impl IntersectionStore {
    fn new() -> Self {
        Self {
            p1q2: Vec::new(),
            x12: Vec::new(),
            v12: Vec::new(),
        }
    }

    fn add(&mut self, p1q2: [i32; 2], x12: i32, v12: glam::DVec3) {
        self.p1q2.push(p1q2);
        self.x12.push(x12);
        self.v12.push(v12);
    }
}

enum Kernel12Fn<'a> {
    TrueTrue(&'a ManifoldImpl, &'a ManifoldImpl),
    TrueFalse(&'a ManifoldImpl, &'a ManifoldImpl),
    FalseTrue(&'a ManifoldImpl, &'a ManifoldImpl),
    FalseFalse(&'a ManifoldImpl, &'a ManifoldImpl),
}

impl<'a> Kernel12Fn<'a> {
    fn call(&self, a1: i32, b2: i32) -> (i32, glam::DVec3) {
        match self {
            Kernel12Fn::TrueTrue(in_p, in_q) => kernel12_true_true(a1, b2, in_p, in_q),
            Kernel12Fn::TrueFalse(in_p, in_q) => kernel12_true_false(a1, b2, in_p, in_q),
            Kernel12Fn::FalseTrue(in_p, in_q) => kernel12_false_true(a1, b2, in_p, in_q),
            Kernel12Fn::FalseFalse(in_p, in_q) => kernel12_false_false(a1, b2, in_p, in_q),
        }
    }
}

/// Recorder for collecting edge-face intersections
/// Uses Mutex for thread-safe interior mutability
struct Kernel12RecorderImpl<'a> {
    store: Mutex<IntersectionStore>,
    k12: Kernel12Fn<'a>,
    forward: bool,
}

impl<'a> Kernel12RecorderImpl<'a> {
    fn new(k12: Kernel12Fn<'a>, forward: bool) -> Self {
        Self {
            store: Mutex::new(IntersectionStore::new()),
            k12,
            forward,
        }
    }

    fn into_intersections(self) -> Intersections {
        let store = self.store.into_inner().unwrap();
        Intersections {
            p1q2: store.p1q2,
            x12: store.x12,
            v12: store.v12,
        }
    }
}

impl<'a> Recorder for Kernel12RecorderImpl<'a> {
    fn record(&self, query_idx: i32, leaf_idx: i32) {
        let (x12, v12) = self.k12.call(query_idx, leaf_idx);
        if v12.x.is_finite() {
            let p1q2 = if self.forward {
                [query_idx, leaf_idx]
            } else {
                [leaf_idx, query_idx]
            };
            let mut store = self.store.lock().unwrap();
            store.add(p1q2, x12, v12);
        }
    }
}

/// Internal Intersect12 implementation for specific expandP and forward combination
/// C++ Reference: `submodules/manifold/src/boolean3.cpp:378-415`
fn intersect12_inner<const EXPAND_P: bool, const FORWARD: bool>(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
) -> Intersections {
    // a: 1 (edge), b: 2 (face)
    let a = if FORWARD { in_p } else { in_q };
    let b = if FORWARD { in_q } else { in_p };

    let k12: Kernel12Fn<'_> = if EXPAND_P {
        if FORWARD {
            Kernel12Fn::TrueTrue(in_p, in_q)
        } else {
            Kernel12Fn::TrueFalse(in_p, in_q)
        }
    } else {
        if FORWARD {
            Kernel12Fn::FalseTrue(in_p, in_q)
        } else {
            Kernel12Fn::FalseFalse(in_p, in_q)
        }
    };

    let recorder = Kernel12RecorderImpl::new(k12, FORWARD);

    // Get the collider from mesh b and run collision detection
    if let Some(ref collider) = b.collider {
        let num_halfedges = a.halfedge.len() as i32;
        let query = EdgeBBoxQuery::new(a);
        collider.collisions::<false, _, _>(&query, num_halfedges, &recorder, true);
    }

    // Get the intersections from the recorder
    let intersections = recorder.into_intersections();

    // Sort p1q2 according to edges
    // DETERMINISM: Rust's sort_by is stable, so equal elements preserve relative order
    // This restores determinism after potentially nondeterministic Mutex-based collection
    let mut i12: Vec<usize> = (0..intersections.p1q2.len()).collect();
    let index = if FORWARD { 0 } else { 1 };
    i12.sort_by(|&a_idx, &b_idx| {
        let a_pair = &intersections.p1q2[a_idx];
        let b_pair = &intersections.p1q2[b_idx];
        a_pair[index]
            .cmp(&b_pair[index])
            .then_with(|| a_pair[1 - index].cmp(&b_pair[1 - index]))
    });

    // Permute the arrays according to sorted indices
    let mut sorted_p1q2 = Vec::with_capacity(intersections.p1q2.len());
    let mut sorted_x12 = Vec::with_capacity(intersections.x12.len());
    let mut sorted_v12 = Vec::with_capacity(intersections.v12.len());

    for &idx in &i12 {
        sorted_p1q2.push(intersections.p1q2[idx]);
        sorted_x12.push(intersections.x12[idx]);
        sorted_v12.push(intersections.v12[idx]);
    }

    Intersections {
        p1q2: sorted_p1q2,
        x12: sorted_x12,
        v12: sorted_v12,
    }
}

/// Intersect12 with forward=true
/// C++ Reference: `submodules/manifold/src/boolean3.cpp:417-424`
pub fn intersect12_true_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    expand_p: bool,
) -> Intersections {
    if expand_p {
        intersect12_inner::<true, true>(in_p, in_q)
    } else {
        intersect12_inner::<false, true>(in_p, in_q)
    }
}

/// Intersect12 with forward=false
/// C++ Reference: `submodules/manifold/src/boolean3.cpp:417-424`
pub fn intersect12_false_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    expand_p: bool,
) -> Intersections {
    if expand_p {
        intersect12_inner::<true, false>(in_p, in_q)
    } else {
        intersect12_inner::<false, false>(in_p, in_q)
    }
}
