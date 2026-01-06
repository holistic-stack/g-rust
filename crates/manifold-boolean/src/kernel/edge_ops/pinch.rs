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
use manifold_parallel::{auto_policy, for_each, stable_sort};
use manifold_types::next_halfedge;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Split pinched vertices (verts with multiple disconnected edge cycles).
/// Port of C++ Manifold::Impl::SplitPinchedVerts in edge_op.cpp
pub fn split_pinched_verts(mesh: &mut ManifoldImpl) {
    let nb_edges = mesh.halfedge.len();
    if nb_edges == 0 {
        return;
    }

    if nb_edges > 10_000 {
        let num_vert = mesh.vert_pos.len();
        let largest_edge = (0..num_vert)
            .map(|_| AtomicUsize::new(usize::MAX))
            .collect::<Vec<_>>();

        let pinched = Mutex::new(Vec::new());

        for_each(auto_policy(nb_edges, 10000), 0..nb_edges, |i| {
            if mesh.halfedge[i].start_vert < 0 {
                return;
            }
            // Find min edge in cycle
            let mut current = i;
            let mut min_edge = i;
            loop {
                let paired = mesh.halfedge[current].paired_halfedge;
                if paired < 0 { break; } // Should be manifold
                current = next_halfedge(paired) as usize;
                if current < min_edge {
                    min_edge = current;
                }
                if current == i { break; }
            }

            // Now try to record this min_edge for the vertex
            let vert = mesh.halfedge[i].start_vert as usize;

            // Loop for CAS
            let _ = largest_edge[vert].fetch_update(Ordering::SeqCst, Ordering::Relaxed, |val| {
                if val == usize::MAX {
                    Some(min_edge)
                } else if val == min_edge {
                    Some(min_edge) // Same cycle, ok
                } else {
                    None
                }
            });

            let val = largest_edge[vert].load(Ordering::SeqCst);
            if val != usize::MAX && val != min_edge {
                // Pinched. Add both to list.
                let mut p = pinched.lock().unwrap();
                p.push(min_edge);
                p.push(val);
            }
        });

        let mut pinched_vec = pinched.into_inner().unwrap();
        if !pinched_vec.is_empty() {
            stable_sort(auto_policy(pinched_vec.len(), 10000), &mut pinched_vec, |a, b| a.cmp(b));
            pinched_vec.dedup();

            let mut processed_verts = std::collections::HashSet::new();
            for idx in pinched_vec {
                let start_vert = mesh.halfedge[idx].start_vert;
                if processed_verts.contains(&start_vert) {
                    // Duplicate vertex
                    mesh.vert_pos.push(mesh.vert_pos[start_vert as usize]);
                    let new_vert = (mesh.vert_pos.len() - 1) as i32;

                    let mut current = idx;
                    loop {
                        mesh.halfedge[current].start_vert = new_vert;
                        let paired = mesh.halfedge[current].paired_halfedge as usize;
                        mesh.halfedge[paired].end_vert = new_vert;

                        let next = next_halfedge(paired as i32) as usize;
                        current = next;
                        if current == idx { break; }
                    }
                } else {
                    processed_verts.insert(start_vert);
                }
            }
        }

    } else {
        // Serial path
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
                    halfedge_processed[current] = true;
                    mesh.halfedge[current].start_vert = new_vert;
                    let paired = mesh.halfedge[current].paired_halfedge as usize;
                    debug_assert!(paired >= 0 as usize);
                    mesh.halfedge[paired].end_vert = new_vert;

                    let paired = mesh.halfedge[current].paired_halfedge;
                    current = next_halfedge(paired) as usize;

                    if current == i {
                        break;
                    }
                }
            } else {
                vert_processed[vert as usize] = true;

                let mut current = i;
                loop {
                    halfedge_processed[current] = true;
                    let paired = mesh.halfedge[current].paired_halfedge;
                    debug_assert!(paired >= 0);
                    current = next_halfedge(paired) as usize;

                    if current == i {
                        break;
                    }
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
