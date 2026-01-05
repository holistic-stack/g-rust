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
    // Create a simple mesh with polygon faces
    let mut mesh = ManifoldImpl::new();

    // Create vertices for a cube (8 corners)
    mesh.vert_pos = vec![
        glam::DVec3::new(0.0, 0.0, 0.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(0.0, 0.0, 1.0),
    ];

    // Create halfedges for a cube with polygon face representation
    // Each face is a square (4 edges)
    mesh.halfedge = vec![
        manifold_boolean::Halfedge {
            start_vert: 0,
            end_vert: 1,
            paired_halfedge: 5,
            prop_vert: 0,
        },
        manifold_boolean::Halfedge {
            start_vert: 1,
            end_vert: 2,
            paired_halfedge: 6,
            prop_vert: 1,
        },
        manifold_boolean::Halfedge {
            start_vert: 2,
            end_vert: 3,
            paired_halfedge: 7,
            prop_vert: 2,
        },
        manifold_boolean::Halfedge {
            start_vert: 3,
            end_vert: 0,
            paired_halfedge: 8,
            prop_vert: 3,
        },
        manifold_boolean::Halfedge {
            start_vert: 0,
            end_vert: 3,
            paired_halfedge: 0,
            prop_vert: 0,
        },
    ];

    // Create face edge array
    let face_edge = vec![0, 4];

    // Face normals
    mesh.face_normal = vec![
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(0.0, 0.0, -1.0),
    ];

    // Triangulate the faces
    let halfedge_ref = vec![
        manifold_boolean::TriRef {
            mesh_id: 0,
            original_id: 0,
            face_id: 0,
            coplanar_id: 0,
        };
        manifold_boolean::TriRef {
            mesh_id: 0,
            original_id: 0,
            face_id: 1,
            coplanar_id: 1,
        }
    ];

    triangulate_mesh_faces(&mut mesh, &face_edge, &halfedge_ref, true);

    // Verify the result
    assert!(
        !mesh.halfedge.is_empty(),
        "Triangulation should produce triangles"
    );
    assert!(
        mesh.halfedge.len() % 3 == 0,
        "Should have complete triangles"
    );
    assert!(
        mesh.face_normal.len() == mesh.halfedge.len() / 3,
        "Each triangle should have a normal"
    );
}
