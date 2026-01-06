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
use manifold_types::Halfedge;

/// Deduplicate edges by duplicating vertices.
/// Port of C++ Manifold::Impl::DedupeEdges in edge_op.cpp:839-965 (serial path)
pub fn dedupe_edges(mesh: &mut ManifoldImpl) {
    loop {
        let nb_edges = mesh.halfedge.len();
        if nb_edges == 0 {
            break;
        }

        let mut processed = vec![false; nb_edges];
        let mut duplicates = Vec::new();

        for i in 0..nb_edges {
            if processed[i] {
                continue;
            }
            if mesh.halfedge[i].start_vert < 0 || mesh.halfedge[i].end_vert < 0 {
                continue;
            }

            let mut end_verts: Vec<(i32, usize)> = Vec::new();
            // HashMap used only for O(1) lookup, never iterated for output.
            // Determinism preserved: iteration order doesn't affect algorithm results.
            let mut end_vert_set: std::collections::HashMap<i32, usize> =
                std::collections::HashMap::new();

            let mut current = i;
            loop {
                let paired = mesh.halfedge[current].paired_halfedge;
                if paired < 0 {
                    break;
                }
                current = next_halfedge(paired) as usize;
                if current >= mesh.halfedge.len() {
                    break;
                }
                processed[current] = true;

                if mesh.halfedge[current].start_vert < 0 || mesh.halfedge[current].end_vert < 0 {
                    continue;
                }

                let end_v = mesh.halfedge[current].end_vert;

                if end_vert_set.is_empty() {
                    let mut found = false;
                    for pair in end_verts.iter_mut() {
                        if pair.0 == end_v {
                            pair.1 = pair.1.min(current);
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        end_verts.push((end_v, current));
                        if end_verts.len() > 32 {
                            for &(v, idx) in &end_verts {
                                end_vert_set.insert(v, idx);
                            }
                            end_verts.clear();
                        }
                    }
                } else {
                    end_vert_set
                        .entry(end_v)
                        .and_modify(|idx| *idx = (*idx).min(current))
                        .or_insert(current);
                }

                if current == i {
                    break;
                }
            }

            let mut current = i;
            loop {
                let paired = mesh.halfedge[current].paired_halfedge;
                if paired < 0 {
                    break;
                }
                current = next_halfedge(paired) as usize;
                if current >= mesh.halfedge.len() {
                    break;
                }

                if mesh.halfedge[current].start_vert < 0 || mesh.halfedge[current].end_vert < 0 {
                    continue;
                }

                let end_v = mesh.halfedge[current].end_vert;

                let keep = if end_vert_set.is_empty() {
                    end_verts
                        .iter()
                        .find(|(v, _)| *v == end_v)
                        .map(|&(_, idx)| idx)
                        .unwrap_or(current)
                } else {
                    end_vert_set.get(&end_v).copied().unwrap_or(current)
                };

                if current != keep {
                    duplicates.push(current);
                }

                if current == i {
                    break;
                }
            }
        }

        if duplicates.is_empty() {
            break;
        }

        duplicates.sort_unstable();
        duplicates.dedup();

        for idx in &duplicates {
            if *idx < mesh.halfedge.len()
                && mesh.halfedge[*idx].start_vert >= 0
                && mesh.halfedge[*idx].end_vert >= 0
            {
                dedupe_edge(mesh, *idx);
            }
        }
    }
}

/// Deduplicate a single 4-manifold edge by duplicating endVert.
/// Port of C++ Manifold::Impl::DedupeEdge in edge_op.cpp:342-427
fn dedupe_edge(mesh: &mut ManifoldImpl, edge: usize) {
    if edge >= mesh.halfedge.len() {
        return;
    }
    if mesh.halfedge[edge].paired_halfedge < 0 {
        return;
    }

    let start_vert = mesh.halfedge[edge].start_vert;
    let end_vert = mesh.halfedge[edge].end_vert;
    if start_vert < 0 || end_vert < 0 {
        return;
    }

    let next_edge = next_halfedge(edge as i32) as usize;
    let end_prop = mesh.halfedge[next_edge].prop_vert;

    let mut current = {
        let paired = mesh.halfedge[next_edge].paired_halfedge;
        if paired < 0 {
            return;
        }
        paired as usize
    };

    let mut split_with_faces = false;
    let opposite = {
        let paired = mesh.halfedge[next_edge].paired_halfedge;
        if paired < 0 {
            return;
        }
        paired as usize
    };

    loop {
        if current == edge {
            break;
        }

        let vert = mesh.halfedge[current].start_vert;
        if vert == start_vert {
            split_with_faces = true;

            let new_vert = mesh.vert_pos.len() as i32;
            mesh.vert_pos.push(mesh.vert_pos[end_vert as usize]);
            if !mesh.vert_normal.is_empty() {
                mesh.vert_normal.push(mesh.vert_normal[end_vert as usize]);
            }

            current = next_halfedge(current as i32) as usize;
            let paired = mesh.halfedge[current].paired_halfedge;
            if paired < 0 {
                break;
            }
            current = paired as usize;

            update_vert(mesh, new_vert, current, opposite);

            let mut new_halfedge = mesh.halfedge.len();
            let old_face = current / 3;
            let outside_vert = mesh.halfedge[current].start_vert;
            mesh.halfedge.push(Halfedge {
                start_vert: end_vert,
                end_vert: new_vert,
                paired_halfedge: -1,
                prop_vert: end_prop,
            });
            mesh.halfedge.push(Halfedge {
                start_vert: new_vert,
                end_vert: outside_vert,
                paired_halfedge: -1,
                prop_vert: end_prop,
            });
            mesh.halfedge.push(Halfedge {
                start_vert: outside_vert,
                end_vert: end_vert,
                paired_halfedge: -1,
                prop_vert: mesh.halfedge[current].prop_vert,
            });
            pair_up(
                mesh,
                new_halfedge + 2,
                mesh.halfedge[current].paired_halfedge as usize,
            );
            pair_up(mesh, new_halfedge + 1, current);
            if !mesh.mesh_relation.tri_ref.is_empty() {
                let tri_ref = mesh.mesh_relation.tri_ref[old_face];
                mesh.mesh_relation.tri_ref.push(tri_ref);
            }
            if !mesh.face_normal.is_empty() {
                let face = mesh.face_normal[old_face];
                mesh.face_normal.push(face);
            }

            new_halfedge += 3;
            let old_face = opposite / 3;
            let outside_vert = mesh.halfedge[opposite].start_vert;
            mesh.halfedge.push(Halfedge {
                start_vert: new_vert,
                end_vert: end_vert,
                paired_halfedge: -1,
                prop_vert: end_prop,
            });
            mesh.halfedge.push(Halfedge {
                start_vert: end_vert,
                end_vert: outside_vert,
                paired_halfedge: -1,
                prop_vert: end_prop,
            });
            mesh.halfedge.push(Halfedge {
                start_vert: outside_vert,
                end_vert: new_vert,
                paired_halfedge: -1,
                prop_vert: mesh.halfedge[opposite].prop_vert,
            });
            pair_up(
                mesh,
                new_halfedge + 2,
                mesh.halfedge[opposite].paired_halfedge as usize,
            );
            pair_up(mesh, new_halfedge + 1, opposite);
            pair_up(mesh, new_halfedge, new_halfedge - 3);
            if !mesh.mesh_relation.tri_ref.is_empty() {
                let tri_ref = mesh.mesh_relation.tri_ref[old_face];
                mesh.mesh_relation.tri_ref.push(tri_ref);
            }
            if !mesh.face_normal.is_empty() {
                let face = mesh.face_normal[old_face];
                mesh.face_normal.push(face);
            }

            break;
        }

        let next = next_halfedge(current as i32) as usize;
        let paired = mesh.halfedge[next].paired_halfedge;
        if paired < 0 {
            break;
        }
        current = paired as usize;
    }

    if !split_with_faces {
        let new_vert = mesh.vert_pos.len() as i32;
        mesh.vert_pos.push(mesh.vert_pos[end_vert as usize]);
        if !mesh.vert_normal.is_empty() {
            mesh.vert_normal.push(mesh.vert_normal[end_vert as usize]);
        }

        let start_edge = next_halfedge(current as i32) as usize;
        let mut current = start_edge;
        loop {
            mesh.halfedge[current].start_vert = new_vert;
            let paired = mesh.halfedge[current].paired_halfedge;
            debug_assert!(paired >= 0);
            mesh.halfedge[paired as usize].end_vert = new_vert;
            current = next_halfedge(paired) as usize;
            if current == start_edge {
                break;
            }
        }
    }

    let pair = mesh.halfedge[edge].paired_halfedge;
    if pair < 0 {
        return;
    }

    let mut current = {
        let next = next_halfedge(pair) as usize;
        let paired = mesh.halfedge[next].paired_halfedge;
        if paired < 0 {
            return;
        }
        paired as usize
    };

    while current != pair as usize {
        let vert = mesh.halfedge[current].start_vert;
        if vert == end_vert {
            break;
        }
        let next = next_halfedge(current as i32) as usize;
        let paired = mesh.halfedge[next].paired_halfedge;
        if paired < 0 {
            break;
        }
        current = paired as usize;
    }

    if current == pair as usize {
        let new_vert = mesh.vert_pos.len() as i32;
        mesh.vert_pos.push(mesh.vert_pos[end_vert as usize]);
        if !mesh.vert_normal.is_empty() {
            mesh.vert_normal.push(mesh.vert_normal[end_vert as usize]);
        }

        let start_edge = next_halfedge(current as i32) as usize;
        let mut current = start_edge;
        loop {
            mesh.halfedge[current].start_vert = new_vert;
            let paired = mesh.halfedge[current].paired_halfedge;
            debug_assert!(paired >= 0);
            mesh.halfedge[paired as usize].end_vert = new_vert;
            current = next_halfedge(paired) as usize;
            if current == start_edge {
                break;
            }
        }
    }
}

fn pair_up(mesh: &mut ManifoldImpl, edge0: usize, edge1: usize) {
    mesh.halfedge[edge0].paired_halfedge = edge1 as i32;
    mesh.halfedge[edge1].paired_halfedge = edge0 as i32;
}

fn update_vert(mesh: &mut ManifoldImpl, vert: i32, start_edge: usize, end_edge: usize) {
    let mut current = start_edge;
    loop {
        mesh.halfedge[current].end_vert = vert;
        let next = next_halfedge(current as i32) as usize;
        mesh.halfedge[next].start_vert = vert;
        current = mesh.halfedge[next].paired_halfedge as usize;
        debug_assert!(current != start_edge, "infinite loop in decimator!");
        if current == end_edge {
            break;
        }
    }
}
