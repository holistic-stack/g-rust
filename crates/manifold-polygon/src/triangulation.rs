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

use crate::bbox::Rect;
use crate::polygon::{triangulate_convex, PolyVert, PolygonsIdx};
use crate::tree2d::{build_two_d_tree, query_two_d_tree};
use glam::{DVec2, IVec3};
use std::collections::BTreeSet;

/// Triangulates a set of ε-valid polygons. If the input is not
/// ε-valid, triangulation may overlap, but will always return a
/// manifold result that matches the input edge directions.
pub fn triangulate(polys: &PolygonsIdx, epsilon: f64, allow_convex: bool) -> Vec<IVec3> {
    if allow_convex && is_convex(polys, epsilon) {
        return triangulate_convex(polys);
    }

    let mut ear_clip = EarClip::new(polys, epsilon);
    ear_clip.triangulate()
}

/// Index-based triangulation entry point
pub fn triangulate_idx(polys: &PolygonsIdx, epsilon: f64) -> Vec<IVec3> {
    triangulate(polys, epsilon, true)
}

/// Counter-clockwise test for three points
pub fn ccw(p0: DVec2, p1: DVec2, p2: DVec2, tol: f64) -> i32 {
    let v1 = p1 - p0;
    let v2 = p2 - p0;
    let area = v1.x * v2.y - v1.y * v2.x;
    let base2 = v1.length_squared().max(v2.length_squared());

    if area * area * 4.0 <= base2 * tol * tol {
        0
    } else {
        if area > 0.0 {
            1
        } else {
            -1
        }
    }
}

fn determinant_2x2(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

fn is_convex(polys: &PolygonsIdx, epsilon: f64) -> bool {
    for poly in polys {
        if poly.is_empty() { continue; }
        let first_edge = poly[0].pos - poly[poly.len() - 1].pos;
        let mut last_edge = first_edge.normalize();

        for v in 0..poly.len() {
            let edge = if v + 1 < poly.len() {
                poly[v + 1].pos - poly[v].pos
            } else {
                first_edge
            };

            let det = determinant_2x2(last_edge, edge);
            if det <= 0.0 || (det.abs() < epsilon && last_edge.dot(edge) < 0.0) {
                return false;
            }
            last_edge = edge.normalize();
        }
    }
    true
}

#[derive(Clone, Copy, Debug)]
struct Vert {
    mesh_idx: i32,
    cost: f64,
    pos: DVec2,
    right_dir: DVec2,
    left: usize,
    right: usize,
    ear_token: Option<EarToken>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct EarToken {
    cost_bits: u64, // To sort by cost (f64 bits)
    vert_idx: usize,
}

impl EarToken {
    fn new(cost: f64, vert_idx: usize) -> Self {
        Self {
            cost_bits: cost.to_bits(),
            vert_idx,
        }
    }
}

struct IdxCollider {
    points: Vec<PolyVert>,
    itr: Vec<usize>, // Maps PolyVert idx (index in points) to Vert index in polygon_
}

struct EarClip {
    polygon: Vec<Vert>,
    holes: Vec<usize>, // Vert indices
    outers: Vec<usize>,
    simples: Vec<usize>,
    hole2_bbox: std::collections::HashMap<usize, Rect>,
    ears_queue: BTreeSet<EarToken>,
    triangles: Vec<IVec3>,
    b_box: Rect,
    epsilon: f64,
}

impl EarClip {
    fn new(polys: &PolygonsIdx, mut epsilon: f64) -> Self {
        let mut num_vert = 0;
        for poly in polys {
            num_vert += poly.len();
        }

        let mut ec = EarClip {
            polygon: Vec::with_capacity(num_vert + 2 * polys.len()),
            holes: Vec::new(),
            outers: Vec::new(),
            simples: Vec::new(),
            hole2_bbox: std::collections::HashMap::new(),
            ears_queue: BTreeSet::new(),
            triangles: Vec::with_capacity(num_vert),
            b_box: Rect::default(),
            epsilon,
        };

        let starts = ec.initialize(polys);

        for i in 0..ec.polygon.len() {
            ec.clip_if_degenerate(i);
        }

        for start in starts {
            ec.find_start(start);
        }

        if epsilon < 0.0 {
            // Epsilon updated in initialize if needed, relying on self.epsilon
        }

        ec
    }

    fn triangulate(&mut self) -> Vec<IVec3> {
        let holes = self.holes.clone();
        for start in holes {
            self.cut_keyhole(start);
        }

        let simples = self.simples.clone();
        for start in simples {
            self.triangulate_poly(start);
        }

        self.triangles.clone()
    }

    fn initialize(&mut self, polys: &PolygonsIdx) -> Vec<usize> {
        let mut starts = Vec::new();

        for poly in polys {
            if poly.is_empty() { continue; }

            let start_idx = self.polygon.len();
            for (i, v) in poly.iter().enumerate() {
                self.b_box.union_point(v.pos);
                let current = self.polygon.len();
                let left = if i == 0 { start_idx + poly.len() - 1 } else { current - 1 };
                let right = if i == poly.len() - 1 { start_idx } else { current + 1 };

                self.polygon.push(Vert {
                    mesh_idx: v.idx,
                    cost: 0.0,
                    pos: v.pos,
                    right_dir: DVec2::ZERO,
                    left,
                    right,
                    ear_token: None,
                });
            }

            let mut curr = start_idx;
            for _ in 0..poly.len() {
                let right = self.polygon[curr].right;
                self.link(curr, right);
                curr = right;
            }

            starts.push(start_idx);
        }

        if self.epsilon < 0.0 {
            self.epsilon = self.b_box.scale() * 1e-5;
        }

        starts
    }

    fn link(&mut self, left: usize, right: usize) {
        self.polygon[left].right = right;
        self.polygon[right].left = left;
        let delta = self.polygon[right].pos - self.polygon[left].pos;
        let dist = delta.length();
        self.polygon[left].right_dir = if dist > 0.0 { delta / dist } else { DVec2::ZERO };
    }

    fn clipped(&self, v: usize) -> bool {
        self.polygon[self.polygon[v].right].left != v
    }

    fn loop_verts<F>(&self, first: usize, mut func: F) -> usize
    where F: FnMut(usize)
    {
        let mut v = first;
        loop {
            if self.clipped(v) {
                let new_first = self.polygon[self.polygon[v].right].left;
                if !self.clipped(new_first) {
                    v = new_first;
                    if self.polygon[v].right == self.polygon[v].left {
                        return usize::MAX;
                    }
                    func(v);
                }
            } else {
                if self.polygon[v].right == self.polygon[v].left {
                    return usize::MAX;
                }
                func(v);
            }

            v = self.polygon[v].right;
            if v == first { break; }
        }
        v
    }

    fn collect_verts(&self, first: usize) -> Vec<usize> {
        let mut verts = Vec::new();
        self.loop_verts(first, |v| verts.push(v));
        verts
    }

    fn clip_ear(&mut self, ear: usize) {
        let left = self.polygon[ear].left;
        let right = self.polygon[ear].right;
        self.link(left, right);

        let mesh_idx = self.polygon[ear].mesh_idx;
        let left_idx = self.polygon[left].mesh_idx;
        let right_idx = self.polygon[right].mesh_idx;

        if left_idx != mesh_idx && mesh_idx != right_idx && right_idx != left_idx {
            self.triangles.push(IVec3::new(left_idx, mesh_idx, right_idx));
        }
    }

    fn clip_if_degenerate(&mut self, ear: usize) {
        if self.clipped(ear) { return; }
        let left = self.polygon[ear].left;
        let right = self.polygon[ear].right;
        if left == right { return; }

        let edge = self.polygon[right].pos - self.polygon[ear].pos;
        let is_short = edge.length_squared() * 4.0 < self.epsilon * self.epsilon;

        let ccw_val = ccw(self.polygon[left].pos, self.polygon[ear].pos, self.polygon[right].pos, self.epsilon);
        let dot = (self.polygon[left].pos - self.polygon[ear].pos).dot(self.polygon[right].pos - self.polygon[ear].pos);

        if is_short || (ccw_val == 0 && dot > 0.0) {
            self.clip_ear(ear);
            self.clip_if_degenerate(left);
            self.clip_if_degenerate(right);
        }
    }

    fn find_start(&mut self, first: usize) {
        let origin = self.polygon[first].pos;
        let mut start = first;
        let mut max_x = f64::NEG_INFINITY;
        let mut bbox = Rect::default();
        let mut area = 0.0;

        let res = self.loop_verts(first, |v_idx| {
            let v = &self.polygon[v_idx];
            bbox.union_point(v.pos);
            let right = &self.polygon[v.right];
            area += determinant_2x2(v.pos - origin, right.pos - origin);

            if v.pos.x > max_x {
                max_x = v.pos.x;
                start = v_idx;
            }
        });

        if res == usize::MAX { return; }

        let size = bbox.size();
        let min_area = self.epsilon * size.x.max(size.y);

        if max_x.is_finite() && area < -min_area {
            self.holes.push(start);
            self.hole2_bbox.insert(start, bbox);
        } else {
            self.simples.push(start);
            if area > min_area {
                self.outers.push(start);
            }
        }
    }

    fn cut_keyhole(&mut self, start: usize) {
        let bbox = self.hole2_bbox[&start];
        let p = self.polygon[start].pos;

        let on_top = if p.y >= bbox.max.y - self.epsilon { 1 }
                     else if p.y <= bbox.min.y + self.epsilon { -1 }
                     else { 0 };

        let mut connector = usize::MAX;

        for &first in &self.outers {
            self.loop_verts(first, |edge_idx| {
                let edge = &self.polygon[edge_idx];

                let x = self.interp_y2x(edge_idx, p, on_top);
                if x.is_finite() && self.inside_edge(start, edge_idx, true) {
                    if connector == usize::MAX
                        || ccw(DVec2::new(x, p.y), self.polygon[connector].pos, self.polygon[self.polygon[connector].right].pos, self.epsilon) == 1
                        || (if self.polygon[connector].pos.y < edge.pos.y {
                                self.inside_edge(edge_idx, connector, false)
                            } else {
                                !self.inside_edge(connector, edge_idx, false)
                            })
                    {
                        connector = edge_idx;
                    }
                }
            });
        }

        if connector == usize::MAX {
            self.simples.push(start);
            return;
        }

        connector = self.find_closer_bridge(start, connector);
        self.join_polygons(start, connector);
    }

    fn interp_y2x(&self, edge_idx: usize, start_pos: DVec2, on_top: i32) -> f64 {
        let edge = &self.polygon[edge_idx];
        let right = &self.polygon[edge.right];

        if (edge.pos.y - start_pos.y).abs() <= self.epsilon {
            if right.pos.y <= start_pos.y + self.epsilon || on_top == 1 {
                f64::NAN
            } else {
                edge.pos.x
            }
        } else if edge.pos.y < start_pos.y - self.epsilon {
            if right.pos.y > start_pos.y + self.epsilon {
                edge.pos.x + (start_pos.y - edge.pos.y) * (right.pos.x - edge.pos.x) / (right.pos.y - edge.pos.y)
            } else if right.pos.y < start_pos.y - self.epsilon || on_top == -1 {
                f64::NAN
            } else {
                right.pos.x
            }
        } else {
            f64::NAN
        }
    }

    fn inside_edge(&self, vert_idx: usize, tail_idx: usize, to_left: bool) -> bool {
        let p2 = self.epsilon * self.epsilon;
        let mut next_l = self.polygon[self.polygon[vert_idx].left].right;
        let mut next_r = self.polygon[tail_idx].right;
        let mut center = tail_idx;
        let mut last = center;

        while next_l != next_r && tail_idx != next_r
            && next_l != (if to_left { self.polygon[vert_idx].right } else { self.polygon[vert_idx].left })
        {
            let edge_l = self.polygon[next_l].pos - self.polygon[center].pos;
            let l2 = edge_l.length_squared();
            if l2 <= p2 {
                next_l = if to_left { self.polygon[next_l].left } else { self.polygon[next_l].right };
                continue;
            }

            let edge_r = self.polygon[next_r].pos - self.polygon[center].pos;
            let r2 = edge_r.length_squared();
            if r2 <= p2 {
                next_r = self.polygon[next_r].right;
                continue;
            }

            let vec_lr = self.polygon[next_r].pos - self.polygon[next_l].pos;
            let lr2 = vec_lr.length_squared();
            if lr2 <= p2 {
                last = center;
                center = next_l;
                next_l = if to_left { self.polygon[next_l].left } else { self.polygon[next_l].right };
                if next_l == next_r { break; }
                next_r = self.polygon[next_r].right;
                continue;
            }

            let mut convexity = ccw(self.polygon[next_l].pos, self.polygon[center].pos, self.polygon[next_r].pos, self.epsilon);
            if center != last {
                convexity += ccw(self.polygon[last].pos, self.polygon[center].pos, self.polygon[next_l].pos, self.epsilon)
                           + ccw(self.polygon[next_r].pos, self.polygon[center].pos, self.polygon[last].pos, self.epsilon);
            }
            if convexity != 0 {
                return convexity > 0;
            }

            if l2 < r2 {
                center = next_l;
                next_l = if to_left { self.polygon[next_l].left } else { self.polygon[next_l].right };
            } else {
                center = next_r;
                next_r = self.polygon[next_r].right;
            }
            last = center;
        }
        true
    }

    fn find_closer_bridge(&self, start: usize, edge: usize) -> usize {
        let edge_pos = self.polygon[edge].pos;
        let start_pos = self.polygon[start].pos;
        let right = self.polygon[edge].right;
        let right_pos = self.polygon[right].pos;

        let mut connector = if edge_pos.x < start_pos.x { right }
            else if right_pos.x < start_pos.x { edge }
            else if right_pos.y - start_pos.y > start_pos.y - edge_pos.y { edge }
            else { right };

        if (self.polygon[connector].pos.y - start_pos.y).abs() <= self.epsilon {
            return connector;
        }

        let above = if self.polygon[connector].pos.y > start_pos.y { 1.0 } else { -1.0 };

        for &first in &self.outers {
            self.loop_verts(first, |vert| {
                let v_pos = self.polygon[vert].pos;
                let c_pos = self.polygon[connector].pos;

                let inside = above * ccw(start_pos, v_pos, c_pos, self.epsilon) as f64;
                if v_pos.x > start_pos.x - self.epsilon
                    && v_pos.y * above > start_pos.y * above - self.epsilon
                    && (inside > 0.0 || (inside == 0.0 && v_pos.x < c_pos.x && v_pos.y * above < c_pos.y * above))
                    && self.inside_edge(vert, edge, true) && self.is_reflex(vert)
                {
                    connector = vert;
                }
            });
        }

        connector
    }

    fn is_reflex(&self, vert: usize) -> bool {
        let left = self.polygon[vert].left;
        let right = self.polygon[left].right;
        !self.inside_edge(left, right, true)
    }

    fn join_polygons(&mut self, start: usize, connector: usize) {
        let new_start = self.polygon.len();
        self.polygon.push(self.polygon[start]);
        let new_connector = self.polygon.len();
        self.polygon.push(self.polygon[connector]);

        let start_right = self.polygon[start].right;
        let connector_left = self.polygon[connector].left;

        self.polygon[start_right].left = new_start;
        self.polygon[connector_left].right = new_connector;

        self.link(start, connector);
        self.link(new_connector, new_start);

        self.clip_if_degenerate(start);
        self.clip_if_degenerate(new_start);
        self.clip_if_degenerate(connector);
        self.clip_if_degenerate(new_connector);
    }

    fn triangulate_poly(&mut self, start: usize) {
        let vert_collider = self.vert_collider(start);
        if vert_collider.points.is_empty() { return; }

        let mut num_tri = -2i32;
        self.ears_queue.clear();

        let verts = self.collect_verts(start);
        for v in verts {
            self.process_ear(v, &vert_collider);
            num_tri += 1;
        }

        while num_tri > 0 {
            // Pop best ear
            let ear_token = if let Some(&t) = self.ears_queue.iter().next() {
                t
            } else {
                break;
            };
            self.ears_queue.remove(&ear_token);
            let v = ear_token.vert_idx;

            self.clip_ear(v);
            num_tri -= 1;

            let left = self.polygon[v].left;
            let right = self.polygon[v].right;

            self.process_ear(left, &vert_collider);
            self.process_ear(right, &vert_collider);
        }
    }

    fn vert_collider(&self, start: usize) -> IdxCollider {
        let mut points = Vec::new();
        let mut itr = Vec::new();

        self.loop_verts(start, |v| {
            points.push(PolyVert {
                pos: self.polygon[v].pos,
                idx: itr.len() as i32,
            });
            itr.push(v);
        });

        build_two_d_tree(&mut points);
        IdxCollider { points, itr }
    }

    fn process_ear(&mut self, v: usize, collider: &IdxCollider) {
        if let Some(token) = self.polygon[v].ear_token {
            self.ears_queue.remove(&token);
            self.polygon[v].ear_token = None;
        }

        let is_short = {
            let right = self.polygon[v].right;
            let edge = self.polygon[right].pos - self.polygon[v].pos;
            edge.length_squared() * 4.0 < self.epsilon * self.epsilon
        };

        let left = self.polygon[v].left;
        let right = self.polygon[v].right;

        let is_convex = ccw(self.polygon[left].pos, self.polygon[v].pos, self.polygon[right].pos, 2.0 * self.epsilon) >= 0;

        if is_short {
            let token = EarToken::new(f64::NEG_INFINITY, v);
            self.ears_queue.insert(token);
            self.polygon[v].cost = f64::NEG_INFINITY;
            self.polygon[v].ear_token = Some(token);
        } else if is_convex {
            let cost = self.ear_cost(v, collider);
            let token = EarToken::new(cost, v);
            self.ears_queue.insert(token);
            self.polygon[v].cost = cost;
            self.polygon[v].ear_token = Some(token);
        } else {
            self.polygon[v].cost = 1.0;
        }
    }

    fn ear_cost(&self, v: usize, collider: &IdxCollider) -> f64 {
        let left = self.polygon[v].left;
        let right = self.polygon[v].right;

        let mut open_side = self.polygon[left].pos - self.polygon[right].pos;
        let center = (self.polygon[left].pos + self.polygon[right].pos) * 0.5;
        let radius = open_side.length() * 0.5;
        let _scale = 4.0 / open_side.length_squared();
        let _ = open_side.normalize_or_zero();

        let total_cost = self.polygon[left].right_dir.dot(self.polygon[v].right_dir) - 1.0 - self.epsilon;

        if ccw(self.polygon[v].pos, self.polygon[left].pos, self.polygon[right].pos, self.epsilon) == 0 {
            return total_cost;
        }

        let mut ear_box = Rect {
            min: DVec2::new(center.x - radius, center.y - radius),
            max: DVec2::new(center.x + radius, center.y + radius),
        };
        ear_box.union_point(self.polygon[v].pos);
        ear_box.min -= self.epsilon;
        ear_box.max += self.epsilon;

        let mesh_idx = self.polygon[v].mesh_idx;
        let lid = self.polygon[left].mesh_idx;
        let rid = self.polygon[right].mesh_idx;

        let mut max_cost = total_cost;

        query_two_d_tree(&collider.points, ear_box, |point| {
            let test = collider.itr[point.idx as usize];
            if !self.clipped(test) && self.polygon[test].mesh_idx != mesh_idx
                && self.polygon[test].mesh_idx != lid
                && self.polygon[test].mesh_idx != rid
            {
                // Simple cost check approximation for now
                // Full cost check would involve self.polygon[test]
                let test_vert = &self.polygon[test];
                let _ = test_vert; // Suppress unused

                // If we implement full cost:
                // let cost = self.cost(test, open_side, epsilon);
                // if cost > max_cost { max_cost = cost; }
            }
        });

        max_cost
    }
}
