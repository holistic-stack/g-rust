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

use manifold::ManifoldImpl;
use manifold_boolean::triangulate_mesh_faces;

/// Integration test for Boolean operations
///
/// This test verifies that Boolean operations produce valid triangulated meshes.
#[test]
fn test_boolean_triangulation() {
    // Create a simple mesh with polygon faces (not triangles)
    let mut mesh = ManifoldImpl::new();

    // Create vertices for a cube (8 corners)
    mesh.vert_pos = vec![
        glam::DVec3::new(0.0, 0.0, 0.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
    ];

    // Create halfedges for a cube (6 faces, each with 4 edges)
    // Face 0 (bottom): 0-1-2-3
    // Face 1 (front): 0-4-5-1
    // Face 2 (top): 2-3-7-6
    // Face 3 (back): 2-6-7-5
    // Face 4 (left): 2-5-6-1
    // Face 5 (right): 1-6-7-4
    // Face 6 (front): 0-5-7-4
    // Face 7 (back): 1-4-7-5

    mesh.halfedge = vec![
        // Face 0 edges
        manifold_boolean::Halfedge { start_vert: 0, end_vert: 1, paired_halfedge: -1, prop_vert: 0 },
        manifold_boolean::Halfedge { start_vert: 1, end_vert: 2, paired_halfedge: -1, prop_vert: 1 },
        manifold_boolean::Halfedge { start_vert: 2, end_vert: 3, paired_halfedge: -1, prop_vert: 2 },
        manifold_boolean::Halfedge { start_vert: 3, end_vert: 0, paired_halfedge: -1, prop_vert: 3 },

        // Face 1 edges
        manifold_boolean::Halfedge { start_vert: 1, end_vert: 4, paired_halfedge: -1, prop_vert: 1 },
        manifold_boolean::Halfedge { start_vert: 4, end_vert: 5, paired_halfedge: -1, prop_vert: 4 },
        manifold_boolean::Halfedge { start_vert: 5, end_vert: 1, paired_halfedge: -1, prop_vert: 5 },
        manifold_boolean::Halfedge { start_vert: 1, end_vert: 7, paired_halfedge: -1, prop_vert: 1 },
        manifold_boolean::Halfedge { start_vert: 7, end_vert: 6, paired_halfedge: -1, prop_vert: 6 },

        // Face 2 edges
        manifold_boolean::Halfedge { start_vert: 4, end_vert: 7, paired_halfedge: -1, prop_vert: 4 },
        manifold_boolean::Halfedge { start_vert: 7, end_vert: 6, paired_halfedge: -1, prop_vert: 7 },
        manifold_boolean::Halfedge { start_vert: 6, end_vert: 2, paired_halfedge: -1, prop_vert: 6 },
        manifold_boolean::Halfedge { start_vert: 2, end_vert: 6, paired_halfedge: -1, prop_vert: 2 },

        // Face 3 edges
        manifold_boolean::Halfedge { start_vert: 6, end_vert: 7, paired_halfedge: -1, prop_vert: 6 },
        manifold_boolean::Halfedge { start_vert: 7, end_vert: 5, paired_halfedge: -1, prop_vert: 7 },
        manifold_boolean::Halfedge { start_vert: 5, end_vert: 2, paired_halfedge: -1, prop_vert: 5 },
        manifold_boolean::Halfedge { start_vert: 2, end_vert: 5, paired_halfedge: -1, prop_vert: 2 },

        // Face 4 edges
        manifold_boolean::Halfedge { start_vert: 5, end_vert: 6, paired_halfedge: -1, prop_vert: 5 },
        manifold_boolean::Halfedge { start_vert: 6, end_vert: 1, paired_halfedge: -1, prop_vert: 6 },
        manifold_boolean::Halfedge { start_vert: 1, end_vert: 7, paired_halfedge: -1, prop_vert: 1 },
        manifold_boolean::Halfedge { start_vert: 7, end_vert: 4, paired_halfedge: -1, prop_vert: 7 },

        // Face 5 edges
        manifold_boolean::Halfedge { start_vert: 5, end_vert: 7, paired_halfedge: -1, prop_vert: 5 },
        manifold_boolean::Halfedge { start_vert: 7, end_vert: 4, paired_halfedge: -1, prop_vert: 4 },

        // Face 6 edges
        manifold_boolean::Halfedge { start_vert: 6, end_vert: 1, paired_halfedge: -1, prop_vert: 6 },
        manifold_boolean::Halfedge { start_vert: 1, end_vert: 4, paired_halfedge: -1, prop_vert: 4 },
        manifold_boolean::Halfedge { start_vert: 4, end_vert: 0, paired_halfedge: -1, prop_vert: 4 },
    ];

    // Create face edge array - 8 faces, so face_edge has 9 entries (0..8)
    let face_edge: Vec<usize> = (0..=8).collect();
    let halfedge_ref = vec![manifold_boolean::TriRef {
        mesh_id: 0,
        original_id: 0,
        face_id: 0,
        coplanar_id: 0,
    }; 24];

    // Face normals
    mesh.face_normal = vec![
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(0.0, 1.0, -1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        glam::DVec3::new(0.0, -1.0, 0.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        glam::DVec3::new(0.0, 1.0, 1.0),
    ];

    // Set epsilon
    mesh.epsilon = 1e-12;

    // Triangulate faces
    triangulate_mesh_faces(&mut mesh, &face_edge, &halfedge_ref, true);

    // Verify result
    assert!(!mesh.halfedge.is_empty(), "Triangulation should produce halfedges");
    assert!(mesh.halfedge.len() % 3 == 0, "Should have complete triangles");
    assert!(mesh.face_normal.len() == mesh.halfedge.len() / 3, "Each triangle should have a normal");

    println!("Boolean triangulation test passed!");
}
