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

use glam::DVec2;

/// A vertex in a polygon, containing 2D position and index
///
/// C++ Reference: include/manifold/polygon.h:27-30
#[derive(Debug, Clone, Copy, Default)]
pub struct PolyVert {
    /// 2D position of the vertex
    pub pos: DVec2,
    /// Index (e.g. original vertex index)
    pub idx: i32,
}

/// A simple polygon (contiguous, no holes)
///
/// C++ Reference: include/manifold/polygon.h:33-36
pub type SimplePolygon = Vec<PolyVert>;

/// Index-based version of SimplePolygon for triangulation
///
/// C++ Reference: include/manifold/polygon.h:39-42
pub type SimplePolygonIdx = Vec<PolyVert>; // Changed from Vec<i32> to Vec<PolyVert> to match C++ SimplePolygonIdx usage in some contexts?
// Wait, C++ PolygonsIdx is vector of SimplePolygonIdx.
// include/manifold/polygon.h: using SimplePolygonIdx = std::vector<PolyVert>;
// So SimplePolygonIdx IS vector of PolyVert.
// My previous definition was Vec<i32>, which was wrong.

/// A set of simple polygons (may represent multiple polygons and/or holes)
///
/// C++ Reference: include/manifold/polygon.h:45-48
pub type Polygons = Vec<SimplePolygon>;

/// Index-based version of Polygons for triangulation
///
/// C++ Reference: include/manifold/polygon.h:51-54
pub type PolygonsIdx = Vec<SimplePolygonIdx>;

/// Fast convex polygon triangulation
///
/// Triangulates a convex polygon using a fan approach.
/// Assumes the polygon is convex and vertices are ordered counter-clockwise.
///
/// C++ Reference: src/polygon.cpp:209-232
pub fn triangulate_convex(polys: &PolygonsIdx) -> Vec<glam::IVec3> {
    let mut triangles = Vec::new();
    let num_tri = polys.iter().map(|p| if p.len() >= 3 { p.len() - 2 } else { 0 }).sum();
    triangles.reserve(num_tri);

    for poly in polys {
        if poly.len() < 3 {
            continue;
        }

        // Fan triangulation logic from C++ TriangulateConvex (alternating fan)
        let mut i = 0;
        let mut k = poly.len() - 1;
        let mut right = true;

        while i + 1 < k {
            let j = if right { i + 1 } else { k - 1 };
            triangles.push(glam::IVec3::new(
                poly[i].idx,
                poly[j].idx,
                poly[k].idx,
            ));
            if right {
                i = j;
            } else {
                k = j;
            }
            right = !right;
        }
    }

    triangles
}
