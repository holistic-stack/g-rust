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

use crate::helpers::{interpolate, shadows, with_sign};
use crate::kernel::ManifoldImpl;
use glam::DVec2;

/// Shadow01 template function for expandP=false, forward=true
/// From boolean3.cpp lines 83-123
pub fn shadow01_false_true(
    a0: i32,
    b1: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, DVec2) {
    let b1s = in_b.halfedge[b1 as usize].start_vert;
    let b1e = in_b.halfedge[b1 as usize].end_vert;
    let a0x = in_a.vert_pos[a0 as usize].x;
    let b1sx = in_b.vert_pos[b1s as usize].x;
    let b1ex = in_b.vert_pos[b1e as usize].x;
    let a0xp = in_a.vert_normal[a0 as usize].x;
    let b1sxp = in_b.vert_normal[b1s as usize].x;
    let b1exp = in_b.vert_normal[b1e as usize].x;

    let s01 = (if shadows(a0x, b1ex, with_sign(false, a0xp) - b1exp) {
        1
    } else {
        0
    }) - (if shadows(a0x, b1sx, with_sign(false, a0xp) - b1sxp) {
        1
    } else {
        0
    });

    let mut yz01 = DVec2::new(f64::NAN, f64::NAN);

    if s01 != 0 {
        yz01 = interpolate(
            in_b.vert_pos[b1s as usize],
            in_b.vert_pos[b1e as usize],
            in_a.vert_pos[a0 as usize].x,
        );

        let b1pair = in_b.halfedge[b1 as usize].paired_halfedge;
        let dir = in_b.face_normal[(b1 / 3) as usize].y + in_b.face_normal[(b1pair / 3) as usize].y;

        if !shadows(in_a.vert_pos[a0 as usize].y, yz01.x, -dir) {
            return (0, yz01);
        }
    }

    (s01, yz01)
}

/// Shadow01 template function for expandP=false, forward=false
/// From boolean3.cpp lines 83-123
pub fn shadow01_false_false(
    a0: i32,
    b1: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, DVec2) {
    let b1s = in_b.halfedge[b1 as usize].start_vert;
    let b1e = in_b.halfedge[b1 as usize].end_vert;
    let a0x = in_a.vert_pos[a0 as usize].x;
    let b1sx = in_b.vert_pos[b1s as usize].x;
    let b1ex = in_b.vert_pos[b1e as usize].x;
    let a0xp = in_a.vert_normal[a0 as usize].x;
    let b1sxp = in_b.vert_normal[b1s as usize].x;
    let b1exp = in_b.vert_normal[b1e as usize].x;

    let s01 = (if shadows(b1sx, a0x, with_sign(false, b1sxp) - a0xp) {
        1
    } else {
        0
    }) - (if shadows(b1ex, a0x, with_sign(false, b1exp) - a0xp) {
        1
    } else {
        0
    });

    let mut yz01 = DVec2::new(f64::NAN, f64::NAN);

    if s01 != 0 {
        yz01 = interpolate(
            in_b.vert_pos[b1s as usize],
            in_b.vert_pos[b1e as usize],
            in_a.vert_pos[a0 as usize].x,
        );

        let b1pair = in_b.halfedge[b1 as usize].paired_halfedge;
        let dir = in_b.face_normal[(b1 / 3) as usize].y + in_b.face_normal[(b1pair / 3) as usize].y;

        if !shadows(yz01.x, in_a.vert_pos[a0 as usize].y, with_sign(false, dir)) {
            return (0, yz01);
        }
    }

    (s01, yz01)
}

/// Shadow01 template function for expandP=true, forward=true
/// From boolean3.cpp lines 83-123
pub fn shadow01_true_true(
    a0: i32,
    b1: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, DVec2) {
    let b1s = in_b.halfedge[b1 as usize].start_vert;
    let b1e = in_b.halfedge[b1 as usize].end_vert;
    let a0x = in_a.vert_pos[a0 as usize].x;
    let b1sx = in_b.vert_pos[b1s as usize].x;
    let b1ex = in_b.vert_pos[b1e as usize].x;
    let a0xp = in_a.vert_normal[a0 as usize].x;
    let b1sxp = in_b.vert_normal[b1s as usize].x;
    let b1exp = in_b.vert_normal[b1e as usize].x;

    let s01 = (if shadows(a0x, b1ex, with_sign(true, a0xp) - b1exp) {
        1
    } else {
        0
    }) - (if shadows(a0x, b1sx, with_sign(true, a0xp) - b1sxp) {
        1
    } else {
        0
    });

    let mut yz01 = DVec2::new(f64::NAN, f64::NAN);

    if s01 != 0 {
        yz01 = interpolate(
            in_b.vert_pos[b1s as usize],
            in_b.vert_pos[b1e as usize],
            in_a.vert_pos[a0 as usize].x,
        );

        let b1pair = in_b.halfedge[b1 as usize].paired_halfedge;
        let dir = in_b.face_normal[(b1 / 3) as usize].y + in_b.face_normal[(b1pair / 3) as usize].y;

        if !shadows(in_a.vert_pos[a0 as usize].y, yz01.x, -dir) {
            return (0, yz01);
        }
    }

    (s01, yz01)
}

/// Shadow01 template function for expandP=true, forward=false
/// From boolean3.cpp lines 83-123
pub fn shadow01_true_false(
    a0: i32,
    b1: i32,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
) -> (i32, DVec2) {
    let b1s = in_b.halfedge[b1 as usize].start_vert;
    let b1e = in_b.halfedge[b1 as usize].end_vert;
    let a0x = in_a.vert_pos[a0 as usize].x;
    let b1sx = in_b.vert_pos[b1s as usize].x;
    let b1ex = in_b.vert_pos[b1e as usize].x;
    let a0xp = in_a.vert_normal[a0 as usize].x;
    let b1sxp = in_b.vert_normal[b1s as usize].x;
    let b1exp = in_b.vert_normal[b1e as usize].x;

    let s01 = (if shadows(b1sx, a0x, with_sign(true, b1sxp) - a0xp) {
        1
    } else {
        0
    }) - (if shadows(b1ex, a0x, with_sign(true, b1exp) - a0xp) {
        1
    } else {
        0
    });

    let mut yz01 = DVec2::new(f64::NAN, f64::NAN);

    if s01 != 0 {
        yz01 = interpolate(
            in_b.vert_pos[b1s as usize],
            in_b.vert_pos[b1e as usize],
            in_a.vert_pos[a0 as usize].x,
        );

        let b1pair = in_b.halfedge[b1 as usize].paired_halfedge;
        let dir = in_b.face_normal[(b1 / 3) as usize].y + in_b.face_normal[(b1pair / 3) as usize].y;

        if !shadows(yz01.x, in_a.vert_pos[a0 as usize].y, with_sign(true, dir)) {
            return (0, yz01);
        }
    }

    (s01, yz01)
}
