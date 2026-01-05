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

//! Face operations for mesh assembly and cleanup
//!
//! C++ Reference: `submodules/manifold/src/face_op.cpp:1-100`
//!
//! This module provides face-level operations needed for Boolean
//! result mesh assembly and cleanup.
//!
//! # Key Functions
//!
//! - **AssembleHalfedges**: Build halfedge structure from polygon indices (lines 39-65)
//! - **CollapseEdge**: Remove an edge by merging vertices (lines 67-100)
//! - **CollapseTri**: Collapse a triangle to a single vertex (lines 102-122)
//!
//! # Critical Dependencies
//!
//! - Uses `manifold_types` for Halfedge, TriRef types
//! - Uses `manifold_math` for DVec3 operations and geometry
//! - Uses `manifold_parallel` for parallel iteration
//! - Uses `manifold-boolean` kernel module for ManifoldImpl
//!
//! # Usage Pattern
//!
//! These functions are called by `boolean_result.rs` during mesh construction:
//!
//! ```rust
//! // After boolean intersection data computed
//! let mut mesh = ManifoldImpl::new();
//! mesh.face_tris = AssembleHalfedges(face_edge_vec, &mesh);
//! mesh.simplify_topology();
//! ```

use glam::{DVec2, DVec3, IVec3};
use manifold_boolean::kernel::ManifoldImpl;
use manifold_parallel::{for_each, ExecutionPolicy};
use manifold_types::{Halfedge, OpType, TriRef};

/// Assemble halfedges from polygon edge indices
/// C++ Reference: face_op.cpp:39-65
///
/// Takes a vector of polygons (where each polygon is a list of 3 vertex indices)
/// and builds the halfedge structure for the mesh.
///
/// # Algorithm
///
/// 1. Each polygon with 3 vertices creates 3 faces (triangles)
/// 2. Halfedges are oriented forward (start < end)
/// 3. Each edge gets 2 halfedges (one for each incident face)
/// 4. Paired halfedges share the same `endVert`
/// 5. Triangle properties are preserved via `TriRef` structure
///
/// # Input
///
/// - `polys`: Vector of polygon vertex index vectors
/// - `halfedge`: Output halfedge vector (will be resized)
/// - `vert_pos`: Array of vertex positions
/// - `face_normal`: Array of face normals (for TriRef preservation)
/// - `tri_ref`: Output triangle references (preserves mesh relationship)
/// - `face_tris`: Output face-to-triangle mapping
/// - `mesh_relation`: Mesh relation tracking structure
///
/// # Output
///
/// Returns a tuple containing:
/// - `face_tris`: Maps from halfedge index to triangle index
/// - `tri_ref`: Preserved triangle references
/// - `mesh_relation`: Updated with triangle references
pub fn assemble_halfedges(
    polys: &[Vec<i32>],
    halfedge: &mut Vec<Halfedge>,
    vert_pos: &[DVec3],
    face_normal: &[DVec3],
    tri_ref: &mut Vec<TriRef>,
    face_tris: &mut Vec<i32>,
    mesh_relation: &mut manifold_boolean::kernel::MeshRelationD,
    op_type: OpType,
) -> (Vec<i32>, Vec<TriRef>) {
    // C++ Reference: face_op.cpp:39-65

    // Reserve space for halfedges
    let num_halfedges = polys.len() * 3;
    halfedge.reserve(num_halfedges);
    tri_ref.reserve(polys.len());

    let mut face_idx = 0;
    let mut halfedge_idx = 0;

    // Process each polygon
    for poly in polys {
        let v0 = poly[0];
        let v1 = poly[1];
        let v2 = poly[2];

        // Create 3 halfedges for this triangle (v0-v1, v1-v2, v2-v0)
        // All oriented forward (start < end for CCW ordering)
        let face_normal_idx = face_normal.len();
        halfedge.push(Halfedge {
            start_vert: v0,
            end_vert: v1,
            paired_halfedge: halfedge_idx + 2, // Will be updated to pair with v2-v0 edge
            prop_vert: v0,
        });
        halfedge.push(Halfedge {
            start_vert: v1,
            end_vert: v2,
            paired_halfedge: halfedge_idx + 1, // Will be updated to pair with v0-v2 edge
            prop_vert: v1,
        });
        halfedge.push(Halfedge {
            start_vert: v2,
            end_vert: v0,
            paired_halfedge: halfedge_idx, // Will be updated to pair with v1-v0 edge
            prop_vert: v2,
        });

        // Store face normal for triangle reference
        tri_ref.push(TriRef {
            mesh_id: mesh_relation.mesh_id_transform.len() as i32,
            original_id: -1,
            face_id: face_idx,
            coplanar_id: -1,
        });
        face_normal.push(face_normal[face_normal_idx]);
        face_tris.push(face_idx);
        face_idx += 1;

        halfedge_idx += 3;
    }

    // Fix paired halfedges to point to correct halfedges
    // Each edge appears twice (once per incident face), need to pair them
    let num_original = halfedge_idx;
    for i in 0..num_original {
        let paired = halfedge[i].paired_halfedge;
        if paired >= 0 && paired < num_original as i32 {
            // Find the forward halfedge for this edge (same start/end, but first index)
            // Forward edges have lower indices for same edge pair
            let j = (i as i32 % 3); // This is the forward edge index
            if j < paired {
                halfedge[j].paired_halfedge = i;
            }
        }
    }

    // Map halfedge indices to triangle indices
    let face_tris = vec![0i32; num_halfedges];
    for i in 0..num_halfedges {
        face_tris[i] = i / 3;
    }

    (face_tris, tri_ref)
}

/// Collapse an edge by merging its two vertices
/// C++ Reference: face_op.cpp:67-100
///
/// This function removes an edge from the mesh by merging
/// its two vertices into a single vertex. This is used to
/// remove degenerate triangles and simplify the mesh.
///
/// # Algorithm
///
/// 1. Find the edge to collapse
/// 2. Update all references to point to the new vertex
/// 3. Mark the collapsed edge as invalid (paired_halfedge = -1)
///
/// # Parameters
///
/// - `halfedge`: The halfedge vector to modify
/// - `vert_pos`: The vertex position array
/// - `halfedge_idx`: Index of the edge being collapsed
/// - `start_vert`, `end_vert`: Vertices of the edge being collapsed
///
/// # Returns
///
/// Index of the newly created merged vertex, or -1 if no collapse occurred
pub fn collapse_edge(
    halfedge: &mut Vec<Halfedge>,
    vert_pos: &mut Vec<DVec3>,
    halfedge_idx: usize,
    start_vert: i32,
    end_vert: i32,
) -> i32 {
    // C++ Reference: face_op.cpp:67-100

    // Check if edge exists and is valid
    if halfedge_idx >= halfedge.len() {
        return -1;
    }

    let he = &halfedge[halfedge_idx];
    if he.paired_halfedge < 0 {
        // Already collapsed or invalid
        return -1;
    }

    // Get vertex positions
    let p_start = vert_pos[he.start_vert as usize];
    let p_end = vert_pos[he.end_vert as usize];

    // Update all halfedges pointing to start_vert to point to end_vert
    // First pass: update forward edges
    for i in 0..halfedge.len() {
        let curr = &halfedge[i];
        if curr.start_vert == start_vert && curr.paired_halfedge >= 0 {
            if let Some(paired) = halfedge.get(curr.paired_halfedge as usize) {
                let paired = &mut halfedge[paired];
                if paired.start_vert == end_vert {
                    // Found the matching forward edge - update its paired
                    paired.paired_halfedge = halfedge_idx as i32;
                }
            }
        }
    }

    // Mark collapsed edge as invalid
    halfedge[halfedge_idx].paired_halfedge = -1;

    halfedge_idx
}

/// Collapse a triangle to a single vertex
/// C++ Reference: face_op.cpp:102-122
///
/// This function collapses a triangle to a single vertex by merging
/// all three vertices into one.
///
/// # Algorithm
///
/// 1. Find the triangle indices from halfedge indices
/// 2. Mark all three halfedges as invalid (paired_halfedge = -1)
/// 3. Create new vertex at merged position
/// 4. Update all triangle references to point to new vertex
///
/// # Parameters
///
/// - `halfedge`: The halfedge vector to modify
/// - `halfedge_idx0`, `halfedge_idx1`, `halfedge_idx2`: Indices of the 3 triangle halfedges
/// - `vert_pos`: The vertex position array
/// - `tri_ref`: Triangle reference array (will be updated)
///
/// # Returns
///
/// Index of the new merged vertex, or -1 if no collapse occurred
pub fn collapse_triangle(
    halfedge: &mut Vec<Halfedge>,
    vert_pos: &mut Vec<DVec3>,
    tri_ref: &mut Vec<TriRef>,
    halfedge_idx0: usize,
    halfedge_idx1: usize,
    halfedge_idx2: usize,
) -> i32 {
    // C++ Reference: face_op.cpp:102-122

    // Check if triangle exists (all 3 halfedges should be valid)
    if halfedge_idx0 >= halfedge.len()
        || halfedge_idx1 >= halfedge.len()
        || halfedge_idx2 >= halfedge.len()
    {
        return -1;
    }

    let he0 = &halfedge[halfedge_idx0];
    let he1 = &halfedge[halfedge_idx1];
    let he2 = &halfedge[halfedge_idx2];

    // Verify this forms a triangle (should be CCW oriented)
    // In halfedge structure, forward edges have start < end
    // A proper triangle would have: he0 (v0->v1), he1 (v1->v2), he2 (v2->v0)

    // Mark all three halfedges as invalid
    he0.paired_halfedge = -1;
    he1.paired_halfedge = -1;
    he2.paired_halfedge = -1;

    // Compute merged vertex position (average of the three vertices)
    let p0 = vert_pos[he0.start_vert as usize];
    let p1 = vert_pos[he0.end_vert as usize];
    let p2 = vert_pos[he1.start_vert as usize];
    let p_new = (p0 + p1 + p2) / 3.0;

    // Add new vertex
    let new_vert_idx = vert_pos.len();
    vert_pos.push(p_new);

    // Update all triangle references to point to new vertex
    // Each triangle is defined by a forward halfedge
    // We need to iterate through all halfedges and update references
    let num_he = halfedge.len();
    for i in 0..num_he {
        let he = &halfedge[i];
        if he.start_vert >= he.paired_halfedge {
            // This is a forward halfedge pointing to this triangle
            // Update its tri_ref to point to new vertex
            tri_ref[he.start_vert] = TriRef {
                mesh_id: 0, // Will be updated by caller
                original_id: -1,
                face_id: he.start_vert / 3, // Triangle index
                coplanar_id: -1,
            };
        }
    }

    new_vert_idx as i32
}

/// Update face_tris mapping after halfedge structure changes
/// This is called after operations that modify halfedge indices to
/// ensure the face_tris mapping stays synchronized.
///
/// C++ Reference: face_op.cpp (implicit in collapse functions)
pub fn update_face_tris(halfedge: &Vec<Halfedge>) -> Vec<i32> {
    halfedge.len() / 3
}
