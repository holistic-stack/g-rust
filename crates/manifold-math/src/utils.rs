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

use glam::{DMat3, DVec3};

#[inline]
pub fn next3(i: usize) -> usize {
    match i {
        0 => 1,
        1 => 2,
        2 => 0,
        _ => i + 1,
    }
}

#[inline]
pub fn prev3(i: usize) -> usize {
    match i {
        0 => 2,
        1 => 0,
        2 => 1,
        _ => i + 2,
    }
}

/// Robust 2D orientation test. Returns 1 (CCW), -1 (CW), or 0 (colinear).
pub fn ccw(p0: glam::DVec2, p1: glam::DVec2, p2: glam::DVec2, tol: f64) -> i32 {
    let v1 = p1 - p0;
    let v2 = p2 - p0;
    let area = v1.x * v2.y - v1.y * v2.x;
    let base2 = v1.length_squared().max(v2.length_squared());
    if area * area * 4.0 <= base2 * tol * tol {
        0
    } else if area > 0.0 {
        1
    } else {
        -1
    }
}

/// Robust calculation of barycentric coordinates for a point relative to a triangle.
pub fn get_barycentric(v: DVec3, tri_pos: &DMat3, tolerance: f64) -> DVec3 {
    let p = [tri_pos.col(0), tri_pos.col(1), tri_pos.col(2)];

    for i in 0..3 {
        if (v - p[i]).length_squared() < tolerance * tolerance {
            let mut uvw = DVec3::ZERO;
            uvw[i] = 1.0;
            return uvw;
        }
    }

    let v0 = p[1] - p[0];
    let v1 = p[2] - p[0];
    let v2 = v - p[0];

    let n = v0.cross(v1);
    let area2 = n.length_squared();

    if area2 > tolerance * tolerance * v0.length_squared().max(v1.length_squared()) {
        let u = v2.cross(v1).dot(n) / area2;
        let v_coord = v0.cross(v2).dot(n) / area2;
        DVec3::new(1.0 - u - v_coord, u, v_coord)
    } else {
        // Line or point case
        let e12 = p[2] - p[1];
        let l2 = [
            v0.length_squared(),
            v1.length_squared(),
            e12.length_squared(),
        ];
        let max_l2 = l2[0].max(l2[1]).max(l2[2]);

        if max_l2 < tolerance * tolerance {
            return DVec3::new(1.0, 0.0, 0.0);
        }

        if l2[0] >= l2[1] && l2[0] >= l2[2] {
            let t = v2.dot(v0) / l2[0];
            DVec3::new(1.0 - t, t, 0.0)
        } else if l2[1] >= l2[0] && l2[1] >= l2[2] {
            let t = v2.dot(v1) / l2[1];
            DVec3::new(1.0 - t, 0.0, t)
        } else {
            let t = (v - p[1]).dot(e12) / l2[2];
            DVec3::new(0.0, 1.0 - t, t)
        }
    }
}
