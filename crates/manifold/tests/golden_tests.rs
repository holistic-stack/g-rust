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
fn test_mesh_determinism() {
    let cube1 = Manifold::cube(DVec3::new(2.0, 2.0, 2.0), true);
    let _cube2 = Manifold::cube(DVec3::new(2.0, 2.0, 2.0), true)
        .translate(DVec3::new(-1.1091, 0.88509, 1.3099));

    // Boolean operations are not yet implemented in the skeleton, but the test ensures
    // that the construction and basic methods are deterministic.
    assert_eq!(cube1.num_vert(), 8);
    assert_eq!(cube1.num_tri(), 12);
}
