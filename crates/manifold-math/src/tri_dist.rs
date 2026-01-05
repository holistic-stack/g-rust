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

// Triangle distance module
//!
//! This module implements triangle distance functions from tri_dist.h
//! for computing distances between points and triangles.
//!
//! Ported exactly from C++ tri_dist.h:39-225 (PhysX algorithm)
//! with no simplification, no deviation from C++ reference implementation.

use glam::DVec3;

/// Edge-edge distance result
///
/// Port of C++ EdgeEdgeDist output from tri_dist.h:39
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeEdgeDist {
    pub x: DVec3,
    pub y: DVec3,
}

impl EdgeEdgeDist {
    pub fn new(x: DVec3, y: DVec3) -> Self {
        Self { x, y }
    }
}

/// Returns distance between two line segments.
///
/// Port of C++ EdgeEdgeDist from tri_dist.h:39-83
///
/// C++ Reference: tri_dist.h:39-83
pub fn edge_edge_dist(
    p: DVec3,
    a: DVec3, // seg 1 origin, vector
    q: DVec3,
    b: DVec3, // seg 2 origin, vector
) -> EdgeEdgeDist {
    let t_vec = q - p;
    let a_dot_a = a.dot(a);
    let b_dot_b = b.dot(b);
    let a_dot_b = a.dot(b);
    let a_dot_t = a.dot(t_vec);
    let b_dot_t = b.dot(t_vec);

    // t parameterizes ray (p, a)
    // u parameterizes ray (q, b)

    // Compute t for the closest point on ray (p, a) to ray (q, b)
    let denom = a_dot_a * b_dot_b - a_dot_b * a_dot_b;

    let mut t = if denom != 0.0 {
        ((a_dot_t * b_dot_b - b_dot_t * a_dot_b) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // find u for point on ray (q, b) closest to point at t
    let mut u;
    if b_dot_b != 0.0 {
        u = (t * a_dot_b - b_dot_t) / b_dot_b;

        // if u is on segment (q, b), t and u correspond to closest points,
        // otherwise, clamp u, recompute and clamp t
        if u < 0.0 {
            u = 0.0;
            t = if a_dot_a != 0.0 {
                (a_dot_t / a_dot_a).clamp(0.0, 1.0)
            } else {
                0.0
            };
        } else if u > 1.0 {
            u = 1.0;
            t = if a_dot_a != 0.0 {
                ((a_dot_b + a_dot_t) / a_dot_a).clamp(0.0, 1.0)
            } else {
                0.0
            };
        }
    } else {
        u = 0.0;
        t = if a_dot_a != 0.0 {
            (a_dot_t / a_dot_a).clamp(0.0, 1.0)
        } else {
            0.0
        };
    }
    let x = p + a * t;
    let y = q + b * u;
    EdgeEdgeDist::new(x, y)
}

/// Returns the minimum squared distance between two triangles.
///
/// Port of C++ DistanceTriangleTriangleSquared from tri_dist.h:90-225
///
/// C++ Reference: tri_dist.h:90-225
pub fn distance_triangle_triangle_squared(p: &[DVec3; 3], q: &[DVec3; 3]) -> f64 {
    let mut sv = [DVec3::ZERO; 3];
    sv[0] = p[1] - p[0];
    sv[1] = p[2] - p[1];
    sv[2] = p[0] - p[2];

    let mut tv = [DVec3::ZERO; 3];
    tv[0] = q[1] - q[0];
    tv[1] = q[2] - q[1];
    tv[2] = q[0] - q[2];

    let mut shown_disjoint = false;
    let mut mindd = f64::MAX;

    for i in 0..3 {
        for j in 0..3 {
            let dists = edge_edge_dist(p[i], sv[i], q[j], tv[j]);
            let v = dists.y - dists.x;
            let dd = v.dot(v);

            if dd <= mindd {
                mindd = dd;

                let mut id = i + 2;
                if id >= 3 {
                    id -= 3;
                }
                let z_p = p[id] - dists.x;
                let mut a = z_p.dot(v);

                id = j + 2;
                if id >= 3 {
                    id -= 3;
                }
                let z_q = q[id] - dists.y;
                let mut b = z_q.dot(v);

                if (a <= 0.0) && (b >= 0.0) {
                    return v.dot(v);
                }

                if a <= 0.0 {
                    a = 0.0;
                } else if b > 0.0 {
                    b = 0.0;
                }

                if (mindd - a + b) > 0.0 {
                    shown_disjoint = true;
                }
            }
        }
    }

    let sn = sv[0].cross(sv[1]);
    let snl = sn.dot(sn);

    if snl > 1e-15 {
        let tp = DVec3::new(
            (p[0] - q[0]).dot(sn),
            (p[0] - q[1]).dot(sn),
            (p[0] - q[2]).dot(sn),
        );

        let mut index = -1;
        if (tp.x > 0.0) && (tp.y > 0.0) && (tp.z > 0.0) {
            index = if tp.x < tp.y { 0 } else { 1 };
            if tp.z < tp[index as usize] {
                index = 2;
            }
        } else if (tp.x < 0.0) && (tp.y < 0.0) && (tp.z < 0.0) {
            index = if tp.x > tp.y { 0 } else { 1 };
            if tp.z > tp[index as usize] {
                index = 2;
            }
        }

        if index >= 0 {
            shown_disjoint = true;
            let q_index = q[index as usize];

            let mut v = q_index - p[0];
            let mut z = sn.cross(sv[0]);
            if v.dot(z) > 0.0 {
                v = q_index - p[1];
                z = sn.cross(sv[1]);
                if v.dot(z) > 0.0 {
                    v = q_index - p[2];
                    z = sn.cross(sv[2]);
                    if v.dot(z) > 0.0 {
                        let cp = q_index + sn * (tp[index as usize] / snl);
                        let cq = q_index;
                        return (cp - cq).dot(cp - cq);
                    }
                }
            }
        }
    }

    let tn = tv[0].cross(tv[1]);
    let tnl = tn.dot(tn);

    if tnl > 1e-15 {
        let sp = DVec3::new(
            (q[0] - p[0]).dot(tn),
            (q[0] - p[1]).dot(tn),
            (q[0] - p[2]).dot(tn),
        );

        let mut index = -1;
        if (sp.x > 0.0) && (sp.y > 0.0) && (sp.z > 0.0) {
            index = if sp.x < sp.y { 0 } else { 1 };
            if sp.z < sp[index as usize] {
                index = 2;
            }
        } else if (sp.x < 0.0) && (sp.y < 0.0) && (sp.z < 0.0) {
            index = if sp.x > sp.y { 0 } else { 1 };
            if sp.z > sp[index as usize] {
                index = 2;
            }
        }

        if index >= 0 {
            shown_disjoint = true;
            let p_index = p[index as usize];

            let mut v = p_index - q[0];
            let mut z = tn.cross(tv[0]);
            if v.dot(z) > 0.0 {
                v = p_index - q[1];
                z = tn.cross(tv[1]);
                if v.dot(z) > 0.0 {
                    v = p_index - q[2];
                    z = tn.cross(tv[2]);
                    if v.dot(z) > 0.0 {
                        let cp = p_index;
                        let cq = p_index + tn * (sp[index as usize] / tnl);
                        return (cp - cq).dot(cp - cq);
                    }
                }
            }
        }
    }

    if shown_disjoint {
        mindd
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_edge_dist_colinear() {
        let p = DVec3::ZERO;
        let a = DVec3::X;
        let q = DVec3::new(1.0, 1.0, 0.0);
        let b = DVec3::new(-1.0, 0.0, 0.0);

        let result = edge_edge_dist(p, a, q, b);

        assert_eq!(result.x, DVec3::ZERO);
        assert_eq!(result.y, DVec3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn test_edge_edge_dist_perpendicular() {
        let p = DVec3::ZERO;
        let a = DVec3::X;
        let q = DVec3::new(0.0, 1.0, 0.0);
        let b = DVec3::Y;

        let result = edge_edge_dist(p, a, q, b);

        // Should find closest point at origin
        assert_eq!(result.x, DVec3::ZERO);
        assert_eq!(result.y, DVec3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn test_distance_triangle_triangle_squared_identical() {
        let p = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let q = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];

        let result = distance_triangle_triangle_squared(&p, &q);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_distance_triangle_triangle_squared_parallel() {
        let p = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let q = [
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 1.0),
            DVec3::new(0.0, 1.0, 1.0),
        ];

        let result = distance_triangle_triangle_squared(&p, &q);
        assert_eq!(result, 1.0);
    }

    #[test]
    fn test_distance_triangle_triangle_squared_disjoint() {
        let p = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let q = [
            DVec3::new(0.0, 0.0, 10.0),
            DVec3::new(1.0, 0.0, 10.0),
            DVec3::new(0.0, 1.0, 10.0),
        ];

        let result = distance_triangle_triangle_squared(&p, &q);
        // Distance should be 10.0 (z-distance), squared is 100.0
        assert_eq!(result, 100.0);
    }
}
