// Copyright 2020 The Manifold Authors.
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

use glam::{DVec2, DVec3, DVec4};

/// Apply sign based on boolean condition.
/// From boolean3.cpp line 40
#[inline]
pub fn with_sign(pos: bool, v: f64) -> f64 {
    if pos {
        v
    } else {
        -v
    }
}

/// Interpolate between two 3D points at given x coordinate.
/// Carefully designed to minimize rounding error.
/// From boolean3.cpp lines 42-56
#[inline]
pub fn interpolate(a_l: DVec3, a_r: DVec3, x: f64) -> DVec2 {
    let dx_l = x - a_l.x;
    let dx_r = x - a_r.x;

    debug_assert!(dx_l * dx_r <= 0.0, "Boolean manifold error: not in domain");

    let use_l = dx_l.abs() < dx_r.abs();
    let d_lr = a_r - a_l;
    let lambda = if use_l { dx_l } else { dx_r } / d_lr.x;

    if !lambda.is_finite() || !d_lr.y.is_finite() || !d_lr.z.is_finite() {
        return DVec2::new(a_l.y, a_l.z);
    }

    DVec2::new(
        lambda * d_lr.y + if use_l { a_l.y } else { a_r.y },
        lambda * d_lr.z + if use_l { a_l.z } else { a_r.z },
    )
}

/// Find intersection of two line segments in 3D.
/// Only place where floating-point operations happen in Boolean function.
/// From boolean3.cpp lines 58-77
pub fn intersect(a_l: DVec3, a_r: DVec3, b_l: DVec3, b_r: DVec3) -> DVec4 {
    let dy_l = b_l.y - a_l.y;
    let dy_r = b_r.y - a_r.y;

    debug_assert!(
        dy_l * dy_r <= 0.0,
        "Boolean manifold error: no intersection"
    );

    let use_l = dy_l.abs() < dy_r.abs();
    let dx = a_r.x - a_l.x;
    let mut lambda = if use_l { dy_l } else { dy_r } / (dy_l - dy_r);
    if !lambda.is_finite() {
        lambda = 0.0;
    }

    let a_dy = a_r.y - a_l.y;
    let b_dy = b_r.y - b_l.y;
    let use_a = a_dy.abs() < b_dy.abs();

    let x = lambda * dx + if use_l { a_l.x } else { a_r.x };
    let y = lambda * (if use_a { a_dy } else { b_dy })
        + if use_l {
            if use_a {
                a_l.y
            } else {
                b_l.y
            }
        } else {
            if use_a {
                a_r.y
            } else {
                b_r.y
            }
        };
    let z = lambda * (a_r.z - a_l.z) + if use_l { a_l.z } else { a_r.z };
    let w = lambda * (b_r.z - b_l.z) + if use_l { b_l.z } else { b_r.z };

    DVec4::new(x, y, z, w)
}

/// Check if one value shadows another along given direction.
/// From boolean3.cpp lines 79-81
#[inline]
pub fn shadows(p: f64, q: f64, dir: f64) -> bool {
    if p == q {
        return dir < 0.0;
    }
    p < q
}
