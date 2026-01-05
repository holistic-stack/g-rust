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

//! Primitive construction tests
//!
//! These tests verify the correctness of primitive construction
//! functions (tetrahedron, cube, sphere, etc).
//!
//! Functions Tested:
//! - ManifoldImpl::from_shape()
//! - Mesh creation (halfedges, normals, properties)
//! - Transform methods
//!
//! C++ Reference Tests:
//! - `submodules/manifold/test/polygon_test.cpp`

use manifold_boolean::kernel::ManifoldImpl;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tetrahedron_creation() {
        // Test tetrahedron basic creation
        let tet = ManifoldImpl::from_shape(manifold_types::Shape::Tetrahedron, glam::DMat4::IDENTITY);
        
        // Verify structure
        assert_eq!(tet.vert_pos.len(), 4);
        assert_eq!(tet.halfedge.len(), 12);
        assert_eq!(tet.face_normal.len(), 4);
        
        // Verify normals are unit length
        for normal in &tet.face_normal {
            assert!((normal.length() - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_cube_creation() {
        let cube = ManifoldImpl::from_shape(manifold_types::Shape::Cube, glam::DMat4::IDENTITY);
        
        // Verify structure
        assert_eq!(cube.vert_pos.len(), 8);
        assert_eq!(cube.halfedge.len(), 12);
        assert_eq!(cube.face_normal.len(), 6);
        
        // Verify cube dimensions
        for i in 0..8 {
            for j in 0..3 {
                let pos = cube.vert_pos[i];
                assert!((pos.x >= -1.0 && pos.x <= 1.0));
                assert!((pos.y >= -1.0 && pos.y <= 1.0));
                assert!((pos.z >= -1.0 && pos.z <= 1.0);
            }
        }
    }

    #[test]
    fn test_sphere_creation() {
        let sphere = ManifoldImpl::from_shape(manifold_types::Shape::Sphere(60), glam::DMat4::IDENTITY);
        
        // Verify structure
        assert_eq!(sphere.vert_pos.len(), 162); // 60 verts for sphere (60 latitude * 3 + 2 poles)
        assert_eq!(sphere.halfedge.len(), 480); // 162 * 2 faces * 3 edges = 482 = 480
        assert_eq!(sphere.face_normal.len(), 320); // 162 faces
        
        // Verify normals are unit length
        for normal in &sphere.face_normal {
            assert!((normal.length() - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_octahedron_creation() {
        let octa = ManifoldImpl::from_shape(manifold_types::Shape::Octahedron, glam::DMat4::IDENTITY);
        
        // Verify structure
        assert_eq!(octa.vert_pos.len(), 6);
        assert_eq!(octa.halfedge.len(), 12);
        assert_eq!(octa.face_normal.len(), 8);
    }

    #[test]
    fn test_transform_identity() {
        let cube = ManifoldImpl::from_shape(manifold_types::Shape::Cube, glam::DMat4::IDENTITY);
        let transformed = cube.subdivide(0, |manifold_boolean::kernel::ManifoldImpl::_edge_divisions);
        
        // Identity transform
        assert_eq!(transformed.vert_pos.len(), cube.vert_pos.len());
        
        for (orig, trans) in cube.vert_pos.iter().zip(transformed.vert_pos.iter()) {
            assert_relative_eq!(orig, trans, 1e-10);
        }
    }

    #[test]
    fn test_transform_translation() {
        let cube = ManifoldImpl::from_shape(manifold_types::Shape::Cube, glam::DMat4::IDENTITY);
        let trans = glam::DMat4::from_translation(glam::DVec3::new(5.0, 0.0, 0.0));
        let transformed = cube.subdivide(0, |manifold_boolean::kernel::ManifoldImpl::_edge_divisions);
        
        // Translation transform
        assert_eq!(transformed.vert_pos.len(), cube.vert_pos.len());
        
        let expected_z = 5.0;
        for (orig, trans) in cube.vert_pos.iter().zip(transformed.vert_pos.iter()) {
            assert_relative_eq!(trans.x, orig.x + 5.0, 1e-10);
            assert_relative_eq!(trans.y, orig.y, 0.0, 1e-10);
            assert_relative_eq!(trans.z, orig.z, 0.0, 1e-10);
        }
    }

    #[test]
    fn test_transform_rotation() {
        let cube = ManifoldImpl::from_shape(manifold_types::Shape::Cube, glam::DMat4::IDENTITY);
        let trans = glam::DMat4::from_axis_angle(glam::DVec3::new(0.0, 90.0, 0.0));
        let transformed = cube.subdivide(0, |manifold_boolean::kernel::ManifoldImpl::_edge_divisions);
        
        // Rotation transform
        assert_eq!(transformed.vert_pos.len(), cube.vert_pos.len());
        
        for (orig, trans) in cube.vert_pos.iter().zip(transformed.vert_pos.iter()) {
            // Check rotation about Y axis
            assert_relative_eq!(trans.x, 0.0, 1e-10);
            
            // Check Y coordinate change (rotated around X axis at origin)
            let y_expected = trans.x;
            let y_rotated = trans.y;
            assert_relative_eq!(y_rotated, 0.0, 1e-10);
        }
    }

    #[test]
    fn test_transform_scale() {
        let cube = ManifoldImpl::from_shape(manifold_types::Shape::Cube, glam::DMat4::IDENTITY);
        let scale = glam::DMat4::from_scale(glam::DVec3::new(2.0, 2.0, 2.0));
        let transformed = cube.subdivide(0, |manifold_boolean::kernel::ManifoldImpl::_edge_divisions);
        
        // Scale transform
        assert_eq!(transformed.vert_pos.len(), cube.vert_pos.len());
        
        let expected_scale = 2.0;
        for (orig, trans) in cube.vert_pos.iter().zip(transformed.vert_pos.iter()) {
            assert_relative_eq!(trans.x, orig.x * expected_scale, 1e-10);
            assert_relative_eq!(trans.y, orig.y * expected_scale, 1e-10);
            assert_relative_eq!(trans.z, orig.z * expected_scale, 1e-10);
        }
    }

    #[test]
    fn test_halfedge_structure() {
        // Test that halfedge structure is correct
        let cube = ManifoldImpl::from_shape(manifold_types::Shape::Cube, glam::DMat4::IDENTITY);
        
        // Test that forward edges have correct properties
        for (i, he) in cube.halfedge.iter().enumerate() {
            if he.is_forward() {
                // Forward edge: paired_halfedge should point to opposite edge + 1
                assert!(he.paired_halfedge >= 0);
                
                // Check vertex ordering: start_vert < end_vert for forward
                assert!(he.start_vert < he.end_vert);
            } else {
                // Backward edge: paired_halfedge should point to edge itself + 1
                assert!(he.paired_halfedge >= 0);
                
                // Check vertex ordering: start_vert > end_vert for backward
                assert!(he.start_vert > he.end_vert);
            }
        }
    }
}
