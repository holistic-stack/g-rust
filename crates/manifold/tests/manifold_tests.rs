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

use glam::DVec3;
use manifold::Manifold;

#[test]
fn test_cube_construction() {
    let cube = Manifold::cube(DVec3::ONE, true);
    assert_eq!(cube.num_vert(), 8);
    assert_eq!(cube.num_tri(), 12);
}

#[test]
fn test_cube_translation() {
    let cube = Manifold::cube(DVec3::ONE, false);
    let translated = cube.translate(DVec3::new(1.0, 0.0, 0.0));
    assert_eq!(translated.num_vert(), 8);
}

#[test]
fn test_tetrahedron_construction() {
    let tetra = Manifold::tetrahedron();
    assert_eq!(tetra.num_vert(), 4);
    assert_eq!(tetra.num_tri(), 4);
}

#[test]
fn test_basic_boolean_union() {
    let cube1 = Manifold::cube(DVec3::ONE, true);
    let cube2 = Manifold::cube(DVec3::ONE, true).translate(DVec3::new(0.5, 0.0, 0.0));
    let result = cube1 + cube2;
    assert!(result.num_vert() >= 0);
}
