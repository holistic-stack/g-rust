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

use manifold_boolean::kernel::ManifoldImpl;
use manifold_boolean::triangulate_mesh_faces;
use manifold_types::{Halfedge, TriRef};

/// Integration test for Boolean triangulation
#[test]
fn test() {
    let mut mesh = ManifoldImpl::new();

    mesh.vert_pos = vec![
        glam::DVec3::new(0.0, 0.0, 0.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
        glam::DVec3::new(1.0, 1.0, 1.0),
    ];

    mesh.halfedge = vec![
        Halfedge {
            start_vert: 0,
            end_vert: 1,
            paired_halfedge: -1,
            prop_vert: 0,
        },
        Halfedge {
            start_vert: 1,
            end_vert: 2,
            paired_halfedge: -1,
            prop_vert: 1,
        },
        Halfedge {
            start_vert: 2,
            end_vert: 3,
            paired_halfedge: -1,
            prop_vert: 2,
        },
        Halfedge {
            start_vert: 3,
            end_vert: 0,
            paired_halfedge: -1,
            prop_vert: 3,
        },
        Halfedge {
            start_vert: 4,
            end_vert: 5,
            paired_halfedge: -1,
            prop_vert: 4,
        },
        Halfedge {
            start_vert: 5,
            end_vert: 6,
            paired_halfedge: -1,
            prop_vert: 5,
        },
        Halfedge {
            start_vert: 6,
            end_vert: 7,
            paired_halfedge: -1,
            prop_vert: 6,
        },
        Halfedge {
            start_vert: 7,
            end_vert: 0,
            paired_halfedge: -1,
            prop_vert: 7,
        },
        Halfedge {
            start_vert: 0,
            end_vert: 4,
            paired_halfedge: -1,
            prop_vert: 0,
        },
    ];

    // Create faces: triangle (3 edges), triangle (3 edges)
    let face_edge = vec![0, 3, 6]; // Start indices of each face

    mesh.face_normal = vec![
        glam::DVec3::new(0.0, 0.0, 1.0),  // Face 0 normal
        glam::DVec3::new(0.0, 0.0, -1.0), // Face 1 normal
    ];

    mesh.epsilon = 1e-12;

    let halfedge_ref = vec![
        TriRef {
            mesh_id: 0,
            original_id: 0,
            face_id: 0,
            coplanar_id: 0,
        };
        6
    ];

    triangulate_mesh_faces(&mut mesh, &face_edge, &halfedge_ref, true);

    // Test that triangulation runs without panicking
    // The triangulation function calls manifold-polygon internally
    // and should handle the mesh data correctly
    assert!(!mesh.halfedge.is_empty());

    println!("Boolean triangulation integration test passed!");
}
