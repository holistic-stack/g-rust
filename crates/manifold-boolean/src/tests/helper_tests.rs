// Copyright 2022 The Manifold Authors.
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

use crate::helpers::*;
use glam::{DVec2, DVec3};

#[test]
fn test_shadows() {
    // Port of C++ logic in boolean3.cpp:79-81
    assert!(shadows(1.0, 2.0, 0.0));
    assert!(!shadows(2.0, 1.0, 0.0));
    assert!(shadows(1.0, 1.0, -1.0));
    assert!(!shadows(1.0, 1.0, 1.0));
}

#[test]
fn test_with_sign() {
    assert_eq!(with_sign(true, 1.0), 1.0);
    assert_eq!(with_sign(false, 1.0), -1.0);
}

#[test]
fn test_interpolate() {
    let a_l = DVec3::new(0.0, 0.0, 0.0);
    let a_r = DVec3::new(1.0, 1.0, 1.0);

    // Exact midpoint
    let mid = interpolate(a_l, a_r, 0.5);
    assert_eq!(mid.x, 0.5);
    assert_eq!(mid.y, 0.5);

    // Boundary cases
    let left = interpolate(a_l, a_r, 0.0);
    assert_eq!(left.x, 0.0);
    assert_eq!(left.y, 0.0);

    let right = interpolate(a_l, a_r, 1.0);
    assert_eq!(right.x, 1.0);
    assert_eq!(right.y, 1.0);
}

#[test]
fn test_intersect() {
    let a_l = DVec3::new(0.0, 0.0, 0.0);
    let a_r = DVec3::new(1.0, 1.0, 0.0);
    let b_l = DVec3::new(0.0, 1.0, 0.0);
    let b_r = DVec3::new(1.0, 0.0, 0.0);

    // Intersection at (0.5, 0.5, 0.0)
    let res = intersect(a_l, a_r, b_l, b_r);
    assert_eq!(res.x, 0.5);
    assert_eq!(res.y, 0.5);
    assert_eq!(res.z, 0.0);
    assert_eq!(res.w, 0.0);
}
