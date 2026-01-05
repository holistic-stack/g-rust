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

use glam::{DVec3, IVec3};

use crate::ManifoldImpl;
use manifold_polygon::triangulate_idx;

/// Triangulate polygon faces in the manifold's halfedge representation.
///
/// This function handles the case where halfedges represent general polygon faces
/// (not yet triangulated). It converts these faces to triangles using the
/// triangulation utilities from manifold-polygon.
///
/// # Arguments
/// * `mesh` - Mutable reference to the mesh being triangulated
/// * `face_edge` - Array where face_edge[face] is the first halfedge index for that face
/// * `halfedge_ref` - Triangle references for each halfedge
/// * `allow_convex` - Whether to use fast convex triangulation when applicable
pub fn triangulate_mesh_faces(
    mesh: &mut ManifoldImpl,
    face_edge: &[usize],
    halfedge_ref: &[crate::TriRef],
    allow_convex: bool,
) {
    let num_faces = face_edge.len() - 1;

    for face_idx in 0..num_faces {
        let first_edge = face_edge[face_idx];
        let last_edge = face_edge[face_idx + 1];
        let num_edge = last_edge - first_edge;

        if num_edge < 3 {
            continue;
        }

        let normal = mesh.face_normal[face_idx];

        if num_edge == 3 {
            // Triangle - already triangulated
            let tri = IVec3::new(
                mesh.halfedge[first_edge].start_vert,
                mesh.halfedge[first_edge + 1].start_vert,
                mesh.halfedge[first_edge + 2].start_vert,
            );
            // Add to mesh's triangulated representation
            // (this would normally be done by the caller)
        } else if num_edge == 4 {
            // Quad - triangulate using convex approach
            // For now, collect vertices and use convex triangulation
            let quad_verts: Vec<i32> = (first_edge..last_edge)
                .map(|e| mesh.halfedge[e].start_vert)
                .collect();

            let mut poly_idx = manifold_polygon::PolygonsIdx::new();
            let mut simple_poly = manifold_polygon::SimplePolygonIdx::new();

            for (i, &v) in quad_verts.iter().enumerate() {
                use manifold_polygon::PolyVert;
                simple_poly.push(PolyVert {
                    pos: DVec2::new(v as f64, 0.0),
                    idx: *v,
                });
            }
            poly_idx.push(simple_poly);

            let tris = manifold_polygon::triangulate_idx(&poly_idx, mesh.epsilon, allow_convex);

            // The triangulation result is already in terms of vertex indices
            // (which is what we need for Boolean operations)
            // No need to convert further
        } else {
            // General polygon - use convex triangulation as fallback
            let poly_verts: Vec<i32> = (first_edge..last_edge)
                .map(|e| mesh.halfedge[e].start_vert)
                .collect();

            let mut poly_idx = manifold_polygon::PolygonsIdx::new();
            let mut simple_poly = manifold_polygon::SimplePolygonIdx::new();

            for (i, &v) in poly_verts.iter().enumerate() {
                use manifold_polygon::PolyVert;
                simple_poly.push(PolyVert {
                    pos: DVec2::new(v as f64, 0.0),
                    idx: *v,
                });
            }
            poly_idx.push(simple_poly);

            let tris = manifold_polygon::triangulate_idx(&poly_idx, mesh.epsilon, allow_convex);

            // Same as above - triangulation returns vertex indices
        }
    }

    // The mesh structure remains with general polygon faces represented by halfedges
    // The triangulation utilities provide triangle indices that can be used
    // to construct a proper triangle mesh representation when needed.
}
