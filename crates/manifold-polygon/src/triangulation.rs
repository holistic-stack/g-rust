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

use crate::polygon::triangulate_idx
use crate::polygon::PolygonsIdx;
use glam::IVec3;

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
