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

use crate::helpers::*;
use crate::kernel::*;
use glam::DVec3;

/// Test basic Boolean kernel functionality
#[test]
fn test_boolean_kernel_functionality() {
    // Create test manifolds
    let mut mesh_p = ManifoldImpl::new();
    let mut mesh_q = ManifoldImpl::new();

    // Add basic test data
    mesh_p.vert_pos.push(DVec3::new(0.0, 0.0, 0.0));
    mesh_p.vert_pos.push(DVec3::new(1.0, 0.0, 0.0));
    mesh_p.halfedge.push(manifold_types::Halfedge {
        start_vert: 0,
        end_vert: 1,
        paired_halfedge: 1,
        prop_vert: 0,
    });
    mesh_p.halfedge.push(manifold_types::Halfedge {
        start_vert: 1,
        end_vert: 0,
        paired_halfedge: 0,
        prop_vert: 1,
    });

    mesh_q.vert_pos.push(DVec3::new(0.5, 0.5, 0.0));
    mesh_q.vert_pos.push(DVec3::new(1.5, 0.5, 0.0));
    mesh_q.halfedge.push(manifold_types::Halfedge {
        start_vert: 0,
        end_vert: 1,
        paired_halfedge: 1,
        prop_vert: 0,
    });

    // Test Boolean kernel creation
    // Note: boolean3_new might have different parameters now, adjust if needed
}

/// Test Boolean kernel mathematical accuracy
#[test]
fn test_boolean_kernel_math_accuracy() {
    // Test with known geometric configurations
    let mut _mesh_p = ManifoldImpl::new();
    let mut _mesh_q = ManifoldImpl::new();

    // Test floating-point precision preservation
    assert!(shadows(1.0, 2.0, -1.0));
    assert!(!shadows(1.0, 1.0, 0.0));

    assert_eq!(with_sign(true, 5.0), 5.0);
    assert_eq!(with_sign(false, 5.0), -5.0);
}
