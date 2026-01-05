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

use crate::helpers::intersect;
use crate::kernel::kernel02;
use crate::kernel::kernel11;
use crate::kernel::{Intersections, ManifoldImpl};
use glam::{DVec3, IVec3};

pub struct Kernel12Recorder {
    pub intersections: Intersections,
}

impl Kernel12Recorder {
    pub fn new() -> Self {
        Self {
            intersections: Intersections::new(),
        }
    }

    pub fn record(&mut self, p1q2: [i32; 2], x12: i32, v12: DVec3) {
        if v12.x.is_finite() {
            self.intersections.p1q2.push(p1q2);
            self.intersections.x12.push(x12);
            self.intersections.v12.push(v12);
        }
    }
}

pub fn make_simple_recorder() -> Kernel12Recorder {
    Kernel12Recorder::new()
}

pub fn kernel12_false_true(
    a1: i32,
    b2: i32,
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
) -> (i32, DVec3) {
    let mut x12 = 0;
    let mut k = 0;
    let mut xzy_lr0 = [DVec3::ZERO; 2];
    let mut xzy_lr1 = [DVec3::ZERO; 2];
    let mut shadows_state = false;

    let edge_a = in_p.halfedge[a1 as usize];
    for vert_a in [edge_a.start_vert, edge_a.end_vert] {
        let (s, z): (i32, f64) = kernel02::kernel02_false_true(vert_a, b2, in_p, in_q);
        if z.is_finite() {
            x12 += s * if vert_a == edge_a.start_vert { 1 } else { -1 };
            if k < 2 && (k == 0 || (s != 0) != shadows_state) {
                shadows_state = s != 0;
                let mut pos = in_p.vert_pos[vert_a as usize];
                std::mem::swap(&mut pos.y, &mut pos.z);
                xzy_lr0[k] = pos;
                xzy_lr1[k] = pos;
                xzy_lr1[k].y = z;
                k += 1;
            }
        }
    }

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_q.halfedge[b1 as usize];
        let b1_f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge
        };
        let (s, xyzz) = kernel11::kernel11_false(in_p, in_q, a1, b1_f);
        if xyzz.x.is_finite() {
            x12 -= s * if edge_b.is_forward() { 1 } else { -1 };
            if k < 2 && (k == 0 || (s != 0) != shadows_state) {
                shadows_state = s != 0;
                xzy_lr0[k] = DVec3::new(xyzz.x, xyzz.z, xyzz.y);
                xzy_lr1[k] = DVec3::new(xyzz.x, xyzz.w, xyzz.y);
                k += 1;
            }
        }
    }

    if x12 == 0 {
        (0, DVec3::splat(f64::NAN))
    } else {
        debug_assert!(k == 2, "Boolean manifold error: v12");
        let xzyy = intersect(xzy_lr0[0], xzy_lr0[1], xzy_lr1[0], xzy_lr1[1]);
        (x12, DVec3::new(xzyy.x, xzyy.z, xzyy.y))
    }
}

pub fn kernel12_true_true(
    a1: i32,
    b2: i32,
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
) -> (i32, DVec3) {
    let mut x12 = 0;
    let mut k = 0;
    let mut xzy_lr0 = [DVec3::ZERO; 2];
    let mut xzy_lr1 = [DVec3::ZERO; 2];
    let mut shadows_state = false;

    let edge_a = in_p.halfedge[a1 as usize];
    for vert_a in [edge_a.start_vert, edge_a.end_vert] {
        let (s, z): (i32, f64) = kernel02::kernel02_true_true(vert_a, b2, in_p, in_q);
        if z.is_finite() {
            x12 += s * if vert_a == edge_a.start_vert { 1 } else { -1 };
            if k < 2 && (k == 0 || (s != 0) != shadows_state) {
                shadows_state = s != 0;
                let mut pos = in_p.vert_pos[vert_a as usize];
                std::mem::swap(&mut pos.y, &mut pos.z);
                xzy_lr0[k] = pos;
                xzy_lr1[k] = pos;
                xzy_lr1[k].y = z;
                k += 1;
            }
        }
    }

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_q.halfedge[b1 as usize];
        let b1_f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge
        };
        let (s, xyzz) = kernel11::kernel11_true(in_p, in_q, a1, b1_f);
        if xyzz.x.is_finite() {
            x12 -= s * if edge_b.is_forward() { 1 } else { -1 };
            if k < 2 && (k == 0 || (s != 0) != shadows_state) {
                shadows_state = s != 0;
                xzy_lr0[k] = DVec3::new(xyzz.x, xyzz.z, xyzz.y);
                xzy_lr1[k] = DVec3::new(xyzz.x, xyzz.w, xyzz.y);
                k += 1;
            }
        }
    }

    if x12 == 0 {
        (0, DVec3::splat(f64::NAN))
    } else {
        debug_assert!(k == 2, "Boolean manifold error: v12");
        let xzyy = intersect(xzy_lr0[0], xzy_lr0[1], xzy_lr1[0], xzy_lr1[1]);
        (x12, DVec3::new(xzyy.x, xzyy.z, xzyy.y))
    }
}
