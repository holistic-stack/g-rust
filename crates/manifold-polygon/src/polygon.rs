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
use manifold_types::TriRef;

/// A vertex in a polygon, containing 2D position and reference to original vertex
///
/// C++ Reference: include/manifold/polygon.h:27-30
#[derive(Debug, Clone, Copy)]
pub struct PolyVert {
    /// 2D position of the vertex
    pub pos: DVec2,
    /// Reference to original vertex in 3D mesh
    pub ref_vert: TriRef,
}

/// A simple polygon (contiguous, no holes)
///
/// C++ Reference: include/manifold/polygon.h:33-36
pub type SimplePolygon = Vec<PolyVert>;

/// Index-based version of SimplePolygon for triangulation
///
/// C++ Reference: include/manifold/polygon.h:39-42
pub type SimplePolygonIdx = Vec<i32>;

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

    for poly in polys {
        if poly.len() < 3 {
            continue;
        }

        // Fan triangulation from first vertex
        for i in 1..poly.len() - 1 {
            triangles.push(glam::IVec3::new(
                poly[0],
                poly[i as usize],
                poly[(i + 1) as usize],
            ));
        }
    }

    triangles
}
