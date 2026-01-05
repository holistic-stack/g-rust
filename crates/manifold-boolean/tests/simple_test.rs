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
        manifold_boolean::Halfedge {
            start_vert: 0,
            end_vert: 1,
            paired_halfedge: -1,
            prop_vert: 0,
        },
        manifold_boolean::Halfedge {
            start_vert: 1,
            end_vert: 2,
            paired_halfedge: -1,
            prop_vert: 1,
        },
        manifold_boolean::Halfedge {
            start_vert: 2,
            end_vert: 3,
            paired_halfedge: -1,
            prop_vert: 2,
        },
        manifold_boolean::Halfedge {
            start_vert: 3,
            end_vert: 0,
            paired_halfedge: -1,
            prop_vert: 3,
        },
        manifold_boolean::Halfedge {
            start_vert: 4,
            end_vert: 5,
            paired_halfedge: -1,
            prop_vert: 4,
        },
        manifold_boolean::Halfedge {
            start_vert: 5,
            end_vert: 6,
            paired_halfedge: -1,
            prop_vert: 5,
        },
        manifold_boolean::Halfedge {
            start_vert: 6,
            end_vert: 7,
            paired_halfedge: -1,
            prop_vert: 6,
        },
        manifold_boolean::Halfedge {
            start_vert: 7,
            end_vert: 0,
            paired_halfedge: -1,
            prop_vert: 7,
        },
        manifold_boolean::Halfedge {
            start_vert: 0,
            end_vert: 4,
            paired_halfedge: -1,
            prop_vert: 0,
        },
    ];

    let face_edge = vec![0, 1, 2, 3, 4, 5, 6];

    mesh.face_normal = vec![
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(0.0, 0.0, -1.0),
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, -1.0, 1.0),
    ];

    mesh.epsilon = 1e-12;

    let halfedge_ref = vec![
        manifold_boolean::TriRef {
            mesh_id: 0,
            original_id: 0,
            face_id: 0,
            coplanar_id: 0,
        };
        6
    ];

    triangulate_mesh_faces(&mut mesh, &face_edge, &halfedge_ref, true);

    assert!(!mesh.halfedge.is_empty());
    assert!(mesh.halfedge.len() % 3 == 0);
    assert!(mesh.face_normal.len() == mesh.halfedge.len() / 3);

    println!("Boolean triangulation integration test passed!");
}
