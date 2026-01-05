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

use super::ManifoldImpl;
use crate::helpers::{intersect, shadows, with_sign};
use crate::kernel::shadow01::shadow01_functions as shadow01;

/// Kernel11 template function for expandP=false
/// From boolean3.cpp lines 125-189
pub fn kernel11_false(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    p1: i32,
    q1: i32,
) -> (i32, glam::DVec4) {
    let mut xyzz11 = glam::DVec4::new(f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    let mut s11;

    // For pRL[k], qRL[k], k==0 is left and k==1 is right.
    let mut k = 0i32;
    let mut p_rl = [glam::DVec3::ZERO, glam::DVec3::ZERO];
    let mut q_rl = [glam::DVec3::ZERO, glam::DVec3::ZERO];

    // Either left or right must shadow, but not both. This ensures
    // intersection is between left and right.
    let mut shadows_val = false;
    s11 = 0;

    let p0 = [
        in_p.halfedge[p1 as usize].start_vert,
        in_p.halfedge[p1 as usize].end_vert,
    ];

    for i in 0..=1 {
        let (s01, yz01) = shadow01::shadow01_false_true(p0[i], q1, in_p, in_q);

        // If the value is NaN, then these do not overlap.
        if yz01.x.is_finite() {
            s11 += s01 * if i == 0 { -1 } else { 1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadows_val) {
                shadows_val = s01 != 0;
                p_rl[k as usize] = in_p.vert_pos[p0[i] as usize];
                q_rl[k as usize] = glam::DVec3::new(p_rl[k as usize].x, yz01.x, yz01.y);
                k += 1;
            }
        }
    }

    let q0 = [
        in_q.halfedge[q1 as usize].start_vert,
        in_q.halfedge[q1 as usize].end_vert,
    ];

    for i in 0..=1 {
        let (s10, yz10) = shadow01::shadow01_false_false(q0[i], p1, in_q, in_p);

        // If the value is NaN, then these do not overlap.
        if yz10.x.is_finite() {
            s11 += s10 * if i == 0 { -1 } else { 1 };
            if k < 2 && (k == 0 || (s10 != 0) != shadows_val) {
                shadows_val = s10 != 0;
                q_rl[k as usize] = in_q.vert_pos[q0[i] as usize];
                p_rl[k as usize] = glam::DVec3::new(q_rl[k as usize].x, yz10.x, yz10.y);
                k += 1;
            }
        }
    }

    if s11 != 0 {
        debug_assert!(k == 2, "Boolean manifold error: s11");
        xyzz11 = intersect(p_rl[0], p_rl[1], q_rl[0], q_rl[1]);

        let p1pair = in_p.halfedge[p1 as usize].paired_halfedge;
        let dir_p =
            in_p.face_normal[(p1 / 3) as usize].z + in_p.face_normal[(p1pair / 3) as usize].z;

        let q1pair = in_q.halfedge[q1 as usize].paired_halfedge;
        let dir_q =
            in_q.face_normal[(q1 / 3) as usize].z + in_q.face_normal[(q1pair / 3) as usize].z;

        if !shadows(xyzz11.z, xyzz11.w, with_sign(false, dir_p) - dir_q) {
            s11 = 0;
        }
    }

    (s11, xyzz11)
}

/// Kernel11 template function for expandP=true
/// From boolean3.cpp lines 125-189
pub fn kernel11_true(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    p1: i32,
    q1: i32,
) -> (i32, glam::DVec4) {
    let mut xyzz11 = glam::DVec4::new(f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    let mut s11;

    // For pRL[k], qRL[k], k==0 is left and k==1 is right.
    let mut k = 0i32;
    let mut p_rl = [glam::DVec3::ZERO, glam::DVec3::ZERO];
    let mut q_rl = [glam::DVec3::ZERO, glam::DVec3::ZERO];

    // Either left or right must shadow, but not both. This ensures
    // intersection is between left and right.
    let mut shadows_val = false;
    s11 = 0;

    let p0 = [
        in_p.halfedge[p1 as usize].start_vert,
        in_p.halfedge[p1 as usize].end_vert,
    ];

    for i in 0..=1 {
        let (s01, yz01) = shadow01::shadow01_true_true(p0[i], q1, in_p, in_q);

        // If the value is NaN, then these do not overlap.
        if yz01.x.is_finite() {
            s11 += s01 * if i == 0 { -1 } else { 1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadows_val) {
                shadows_val = s01 != 0;
                p_rl[k as usize] = in_p.vert_pos[p0[i] as usize];
                q_rl[k as usize] = glam::DVec3::new(p_rl[k as usize].x, yz01.x, yz01.y);
                k += 1;
            }
        }
    }

    let q0 = [
        in_q.halfedge[q1 as usize].start_vert,
        in_q.halfedge[q1 as usize].end_vert,
    ];

    for i in 0..=1 {
        let (s10, yz10) = shadow01::shadow01_true_false(q0[i], p1, in_q, in_p);

        // If the value is NaN, then these do not overlap.
        if yz10.x.is_finite() {
            s11 += s10 * if i == 0 { -1 } else { 1 };
            if k < 2 && (k == 0 || (s10 != 0) != shadows_val) {
                shadows_val = s10 != 0;
                q_rl[k as usize] = in_q.vert_pos[q0[i] as usize];
                p_rl[k as usize] = glam::DVec3::new(q_rl[k as usize].x, yz10.x, yz10.y);
                k += 1;
            }
        }
    }

    if s11 != 0 {
        debug_assert!(k == 2, "Boolean manifold error: s11");
        xyzz11 = intersect(p_rl[0], p_rl[1], q_rl[0], q_rl[1]);

        let p1pair = in_p.halfedge[p1 as usize].paired_halfedge;
        let dir_p =
            in_p.face_normal[(p1 / 3) as usize].z + in_p.face_normal[(p1pair / 3) as usize].z;

        let q1pair = in_q.halfedge[q1 as usize].paired_halfedge;
        let dir_q =
            in_q.face_normal[(q1 / 3) as usize].z + in_q.face_normal[(q1pair / 3) as usize].z;

        if !shadows(xyzz11.z, xyzz11.w, with_sign(true, dir_p) - dir_q) {
            s11 = 0;
        }
    }

    (s11, xyzz11)
}
