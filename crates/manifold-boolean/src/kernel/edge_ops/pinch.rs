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

use crate::kernel::ManifoldImpl;
use manifold_types::next_halfedge;

/// Split pinched vertices (verts with multiple disconnected edge cycles).
/// Port of C++ Manifold::Impl::SplitPinchedVerts in edge_op.cpp:735-837 (serial path)
pub fn split_pinched_verts(mesh: &mut ManifoldImpl) {
    let nb_edges = mesh.halfedge.len();
    if nb_edges == 0 {
        return;
    }

    let num_vert = mesh.vert_pos.len();
    let mut vert_processed = vec![false; num_vert];
    let mut halfedge_processed = vec![false; nb_edges];

    for i in 0..nb_edges {
        if halfedge_processed[i] {
            continue;
        }

        let vert = mesh.halfedge[i].start_vert;
        if vert < 0 {
            continue;
        }

        if vert_processed[vert as usize] {
            mesh.vert_pos.push(mesh.vert_pos[vert as usize]);
            let new_vert = (mesh.vert_pos.len() - 1) as i32;

            let mut current = i;
            loop {
                let paired = mesh.halfedge[current].paired_halfedge;
                debug_assert!(paired >= 0);
                current = next_halfedge(paired) as usize;

                halfedge_processed[current] = true;
                mesh.halfedge[current].start_vert = new_vert;
                let paired = mesh.halfedge[current].paired_halfedge;
                debug_assert!(paired >= 0);
                mesh.halfedge[paired as usize].end_vert = new_vert;

                if current == i {
                    break;
                }
            }
        } else {
            vert_processed[vert as usize] = true;

            let mut current = i;
            loop {
                let paired = mesh.halfedge[current].paired_halfedge;
                debug_assert!(paired >= 0);
                current = next_halfedge(paired) as usize;

                halfedge_processed[current] = true;

                if current == i {
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_pinched_verts_empty_mesh() {
        use crate::kernel::ManifoldImpl;
        use glam::DVec3;
        use manifold_types::Halfedge;

        let mut mesh = ManifoldImpl {
            vert_pos: vec![],
            vert_normal: vec![],
            halfedge: vec![],
            face_normal: vec![],
            num_prop: 0,
            properties: vec![],
            mesh_relation: Default::default(),
            collider: None,
            bbox: Default::default(),
            epsilon: 1e-6,
            tolerance: 1e-6,
            halfedge_tangent: vec![],
        };

        split_pinched_verts(&mut mesh);

        assert_eq!(mesh.vert_pos.len(), 0);
        assert_eq!(mesh.halfedge.len(), 0);
    }

    #[test]
    fn test_split_pinched_verts_simple_triangle() {
        use crate::kernel::ManifoldImpl;
        use glam::DVec3;
        use manifold_types::Halfedge;

        let mut mesh = ManifoldImpl {
            vert_pos: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(2.0, 0.0, 0.0),
            ],
            vert_normal: vec![],
            halfedge: vec![
                Halfedge {
                    start_vert: 0,
                    end_vert: 1,
                    paired_halfedge: 3,
                    prop_vert: 0,
                },
                Halfedge {
                    start_vert: 1,
                    end_vert: 2,
                    paired_halfedge: 4,
                    prop_vert: 0,
                },
                Halfedge {
                    start_vert: 2,
                    end_vert: 0,
                    paired_halfedge: 5,
                    prop_vert: 0,
                },
                Halfedge {
                    start_vert: 1,
                    end_vert: 0,
                    paired_halfedge: 0,
                    prop_vert: 0,
                },
                Halfedge {
                    start_vert: 2,
                    end_vert: 1,
                    paired_halfedge: 1,
                    prop_vert: 0,
                },
                Halfedge {
                    start_vert: 0,
                    end_vert: 2,
                    paired_halfedge: 2,
                    prop_vert: 0,
                },
            ],
            face_normal: vec![],
            num_prop: 0,
            properties: vec![],
            mesh_relation: Default::default(),
            collider: None,
            bbox: Default::default(),
            epsilon: 1e-6,
            tolerance: 1e-6,
            halfedge_tangent: vec![],
        };

        let original_vert_count = mesh.vert_pos.len();
        split_pinched_verts(&mut mesh);

        assert_eq!(mesh.vert_pos.len(), original_vert_count);
    }
}
