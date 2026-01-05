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
use crate::helpers::{interpolate, shadows, with_sign};
use crate::kernel::shadow01::shadow01_functions as shadow01;
use glam::DVec3;

pub fn kernel02_false_true(
    a0: i32,
    b2: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, f64) {
    let mut s02 = 0;
    let mut k = 0;
    let mut yzz_rl = [DVec3::ZERO; 2];
    let mut shadows_state = false;

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_b.halfedge[b1 as usize];
        let b1_f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge
        };

        let (s01, yz01) = shadow01::shadow01_false_true(a0, b1_f, in_a, in_b);
        if yz01.x.is_finite() {
            s02 += s01 * if edge_b.is_forward() { -1 } else { 1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadows_state) {
                shadows_state = s01 != 0;
                yzz_rl[k as usize] = DVec3::new(yz01.x, yz01.y, yz01.y);
                k += 1;
            }
        }
    }

    if s02 == 0 {
        (0, f64::NAN)
    } else {
        debug_assert!(k == 2, "Boolean manifold error: s02");
        let vert_pos_a = in_a.vert_pos[a0 as usize];
        let z02 = interpolate(yzz_rl[0], yzz_rl[1], vert_pos_a.y).y;
        if !shadows(vert_pos_a.z, z02, -in_b.face_normal[b2 as usize].z) {
            return (0, z02);
        }
        (s02, z02)
    }
}

pub fn kernel02_false_false(
    a0: i32,
    b2: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, f64) {
    let mut s02 = 0;
    let mut k = 0;
    let mut yzz_rl = [DVec3::ZERO; 2];
    let mut shadows_state = false;

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_b.halfedge[b1 as usize];
        let b1_f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge
        };

        let (s01, yz01) = shadow01::shadow01_false_false(a0, b1_f, in_a, in_b);
        if yz01.x.is_finite() {
            s02 += s01 * if edge_b.is_forward() { 1 } else { -1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadows_state) {
                shadows_state = s01 != 0;
                yzz_rl[k as usize] = DVec3::new(yz01.x, yz01.y, yz01.y);
                k += 1;
            }
        }
    }

    if s02 == 0 {
        (0, f64::NAN)
    } else {
        debug_assert!(k == 2, "Boolean manifold error: s02");
        let vert_pos_a = in_a.vert_pos[a0 as usize];
        let z02 = interpolate(yzz_rl[0], yzz_rl[1], vert_pos_a.y).y;
        if !shadows(
            z02,
            vert_pos_a.z,
            with_sign(false, in_b.face_normal[b2 as usize].z),
        ) {
            return (0, z02);
        }
        (s02, z02)
    }
}

pub fn kernel02_true_true(
    a0: i32,
    b2: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, f64) {
    let mut s02 = 0;
    let mut k = 0;
    let mut yzz_rl = [DVec3::ZERO; 2];
    let mut shadows_state = false;

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_b.halfedge[b1 as usize];
        let b1_f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge
        };

        let (s01, yz01) = shadow01::shadow01_true_true(a0, b1_f, in_a, in_b);
        if yz01.x.is_finite() {
            s02 += s01 * if edge_b.is_forward() { -1 } else { 1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadows_state) {
                shadows_state = s01 != 0;
                yzz_rl[k as usize] = DVec3::new(yz01.x, yz01.y, yz01.y);
                k += 1;
            }
        }
    }

    if s02 == 0 {
        (0, f64::NAN)
    } else {
        debug_assert!(k == 2, "Boolean manifold error: s02");
        let vert_pos_a = in_a.vert_pos[a0 as usize];
        let z02 = interpolate(yzz_rl[0], yzz_rl[1], vert_pos_a.y).y;
        if !shadows(vert_pos_a.z, z02, -in_b.face_normal[b2 as usize].z) {
            return (0, z02);
        }
        (s02, z02)
    }
}

pub fn kernel02_true_false(
    a0: i32,
    b2: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, f64) {
    let mut s02 = 0;
    let mut k = 0;
    let mut yzz_rl = [DVec3::ZERO; 2];
    let mut shadows_state = false;

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_b.halfedge[b1 as usize];
        let b1_f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge
        };

        let (s01, yz01) = shadow01::shadow01_true_false(a0, b1_f, in_a, in_b);
        if yz01.x.is_finite() {
            s02 += s01 * if edge_b.is_forward() { 1 } else { -1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadows_state) {
                shadows_state = s01 != 0;
                yzz_rl[k as usize] = DVec3::new(yz01.x, yz01.y, yz01.y);
                k += 1;
            }
        }
    }

    if s02 == 0 {
        (0, f64::NAN)
    } else {
        debug_assert!(k == 2, "Boolean manifold error: s02");
        let vert_pos_a = in_a.vert_pos[a0 as usize];
        let z02 = interpolate(yzz_rl[0], yzz_rl[1], vert_pos_a.y).y;
        if !shadows(
            z02,
            vert_pos_a.z,
            with_sign(true, in_b.face_normal[b2 as usize].z),
        ) {
            return (0, z02);
        }
        (s02, z02)
    }
}
