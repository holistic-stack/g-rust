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

//! Atomic math unit tests
//!
//! These tests verify correctness of core mathematical functions used
//! throughout the boolean kernel implementation.
//!
//! Functions Tested:
//! - sind(), cosd() - Trigonometric functions
//! - with_sign() - Sign application
//! - interpolate() - Linear interpolation with error minimization
//! - intersect() - Line segment intersection
//! - shadows() - Shadow predicate for symbolic perturbation
//!
//! C++ Reference Tests:
//! - `submodules/manifold/test/trig_test.cpp`
//! - `submodules/manifold/test/boolean_test.cpp` (sections on atomic math)

use manifold_boolean::helpers::{shadows, with_sign};
use manifold_math::{cosd, interpolate, intersect, sind, with_sign};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sind() {
        // Test basic sine values
        assert_relative_eq!(sind(0.0), 0.0, 1e-10);
        assert_relative_eq!(sind(90.0), 1.0, 1e-10);
        assert_relative_eq!(sind(180.0), 0.0, 1e-10);
        assert_relative_eq!(sind(-90.0), -1.0, 1e-10);
        assert_relative_eq!(sind(-180.0), -1.0, 1e-10);

        // Test cosine values
        assert_relative_eq!(cosd(0.0), 1.0, 1e-10);
        assert_relative_eq!(cosd(60.0), 0.5, 1e-10);
        assert_relative_eq!(cosd(90.0), 0.0, 1e-10);

        // Test period behavior
        let result = sind(180.0) + sind(180.0);
        assert!((result - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cosd() {
        // Test symmetry and special values
        assert_relative_eq!(cosd(0.0), 1.0, 1e-10);
        assert_relative_eq!(cosd(60.0), 0.5, 1e-10);
        assert_relative_eq!(cosd(90.0), 0.0, 1e-10);
        assert_relative_eq!(cosd(30.0), 0.8660254, 1e-10);
        assert_relative_eq!(cosd(120.0), 0.5, 1e-10);
    }

    #[test]
    fn test_with_sign() {
        assert_eq!(with_sign(true, 10.0), 10.0);
        assert_eq!(with_sign(true, -10.0), -10.0);
        assert_eq!(with_sign(false, 10.0), 10.0);
        assert_eq!(with_sign(false, -10.0), -10.0);
    }

    #[test]
    fn test_interpolate_basic() {
        // Test midpoint
        let a = glam::DVec2::new(0.0, 0.0);
        let b = glam::DVec2::new(2.0, 2.0);
        let result = interpolate(a, b, 1.0);
        assert_relative_eq!(result.x, 1.0, 1e-10);
        assert_relative_eq!(result.y, 1.0, 1e-10);

        // Test endpoint selection
        let result_b = interpolate(a, b, 2.0);
        assert_relative_eq!(result_b.x, 2.0, 1e-10);
        assert_relative_eq!(result_b.y, 2.0, 1e-10);
    }

    #[test]
    fn test_intersect_orthogonal() {
        // Test vertical line intersection
        let a1 = glam::DVec2::new(0.0, 0.0);
        let b1 = glam::DVec2::new(0.0, 1.0);
        let a2 = glam::DVec2::new(0.0, 1.0);
        let b2 = glam::DVec2::new(1.0, 1.0);

        let result = intersect(a1, b1, a2, b2);
        assert!(result.x.is_finite());
        assert_relative_eq!(result.x, 0.0, 1e-10);
        assert_relative_eq!(result.y, 0.5, 1e-10);
        assert_relative_eq!(result.z, 0.0, 1e-10);
        assert_relative_eq!(result.w, 1.0, 1e-10);
    }

    #[test]
    fn test_intersect_parallel() {
        // Test parallel line intersection
        let a1 = glam::DVec2::new(0.0, 0.0);
        let b1 = glam::DVec2::new(0.0, 2.0);
        let a2 = glam::DVec2::new(1.0, 0.0);
        let b2 = glam::DVec2::new(2.0, 1.0);

        let result = intersect(a1, b1, a2, b2);
        assert!(result.x.is_finite());
        assert_relative_eq!(result.x, 1.0, 1e-10);
        assert_relative_eq!(result.y, 2.0, 1e-10);
        assert_relative_eq!(result.z, 0.0, 1e-10);
        assert_relative_eq!(result.w, 2.0, 1e-10);
    }

    #[test]
    fn test_shadows_equal() {
        // Test equal values - direction matters
        assert!(!shadows(1.0, 1.0, 0.0));
        assert!(!shadows(1.0, 1.0, 1.0));

        assert!(shadows(1.0, 1.0, -1.0));
        assert!(shadows(1.0, 1.0, -1.0));
    }

    #[test]
    fn test_shadows_inequal() {
        // Test inequality
        assert!(shadows(1.0, 0.0, 0.0));
        assert!(!shadows(0.0, 1.0, 0.0));
        assert!(!shadows(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_shadows_nan_cases() {
        // Test NaN handling
        assert!(!shadows(1.0, f64::NAN, 0.0));
        assert!(shards(1.0, 0.0, 0.0)); // Should not panic but handle
    }
}
