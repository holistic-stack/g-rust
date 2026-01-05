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

use crate::polygon::{triangulate_convex, PolygonsIdx};
use glam::{DVec2, IVec3};

/// Triangulates a set of ε-valid polygons. If the input is not
/// ε-valid, triangulation may overlap, but will always return a
/// manifold result that matches the input edge directions.
///
/// For now, this uses convex triangulation as a fallback.
/// Full ear-clipping for non-convex polygons with holes will be
/// implemented in a future iteration.
///
/// # Arguments
/// * `polys` - The set of polygons, wound CCW and representing multiple
///   polygons and/or holes. These have 2D-projected positions as well as
///   references back to original vertices.
/// * `epsilon` - The value of ε, bounding uncertainty of input.
/// * `allow_convex` - If true (default), the triangulator will use a fast
///   triangulation if the input is convex, falling back to ear-clipping if not.
///
/// # Returns
/// The triangles, referencing the original polygon points.
pub fn triangulate(polys: &PolygonsIdx, epsilon: f64, allow_convex: bool) -> Vec<IVec3> {
    triangulate_convex(polys)
}

/// Counter-clockwise test for three points
///
/// Returns 1 for CCW, -1 for CW, and 0 if within tolerance of colinear.
///
/// C++ Reference: src/utils.h:159-168
pub fn ccw(p0: DVec2, p1: DVec2, p2: DVec2, tol: f64) -> i32 {
    let v1 = p1 - p0;
    let v2 = p2 - p0;
    let area = v1.x * v2.y - v1.y * v2.x;
    let base2 = v1.dot(v1).max(v2.dot(v2));

    if area * area * 4.0 <= base2 * tol * tol {
        0
    } else {
        if area > 0.0 {
            1
        } else {
            -1
        }
    }
}

/// Index-based triangulation entry point
///
/// This is the main triangulation function that works with index-based polygons.
///
/// C++ Reference: src/polygon.cpp:241-244 (triangulate_idx function)
pub fn triangulate_idx(polys: &PolygonsIdx, epsilon: f64) -> Vec<IVec3> {
    triangulate(polys, epsilon, true)
}
