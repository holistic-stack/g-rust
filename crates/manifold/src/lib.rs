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

use glam::{DMat4, DVec3, DVec4};
use manifold_boolean::kernel::{create_boolean_result, Boolean3, ManifoldImpl, Shape};
use manifold_types::OpType;
use std::sync::Arc;

pub struct Manifold {
    impl_: Arc<ManifoldImpl>,
}

impl Manifold {
    pub fn new(impl_: ManifoldImpl) -> Self {
        Self {
            impl_: Arc::new(impl_),
        }
    }

    pub fn cube(size: DVec3, center: bool) -> Self {
        if size.x < 0.0 || size.y < 0.0 || size.z < 0.0 || size.length() == 0.0 {
            return Self::invalid();
        }

        let m = DMat4::from_cols_array_2d(&[
            [size.x, 0.0, 0.0, 0.0],
            [0.0, size.y, 0.0, 0.0],
            [0.0, 0.0, size.z, 0.0],
            [
                if center { -size.x / 2.0 } else { 0.0 },
                if center { -size.y / 2.0 } else { 0.0 },
                if center { -size.z / 2.0 } else { 0.0 },
                1.0,
            ],
        ]);

        Self::new(ManifoldImpl::from_shape(Shape::Cube, m))
    }

    pub fn sphere(radius: f64, circular_segments: i32) -> Self {
        if radius <= 0.0 {
            return Self::invalid();
        }

        let n = if circular_segments > 0 {
            (circular_segments + 3) / 4
        } else {
            8
        };

        let mut impl_ = ManifoldImpl::from_shape(Shape::Octahedron, DMat4::IDENTITY);
        impl_.subdivide(|_, _, _| n - 1);

        for v in &mut impl_.vert_pos {
            let cos_v = DVec3::new(
                (std::f64::consts::FRAC_PI_2 * (1.0 - v.x)).cos(),
                (std::f64::consts::FRAC_PI_2 * (1.0 - v.y)).cos(),
                (std::f64::consts::FRAC_PI_2 * (1.0 - v.z)).cos(),
            );
            *v = radius * cos_v.normalize();
            if v.x.is_nan() {
                *v = DVec3::ZERO;
            }
        }

        impl_.finish();
        impl_.initialize_original();
        Self::new(impl_)
    }

    pub fn tetrahedron() -> Self {
        Self::new(ManifoldImpl::from_shape(
            Shape::Tetrahedron,
            DMat4::IDENTITY,
        ))
    }

    pub fn invalid() -> Self {
        Self::new(ManifoldImpl::default())
    }

    pub fn num_vert(&self) -> usize {
        self.impl_.vert_pos.len()
    }

    pub fn num_tri(&self) -> usize {
        self.impl_.halfedge.len() / 3
    }

    pub fn translate(&self, v: DVec3) -> Self {
        let m = DMat4::from_translation(v);
        self.transform(m)
    }

    pub fn scale(&self, v: DVec3) -> Self {
        let m = DMat4::from_scale(v);
        self.transform(m)
    }

    pub fn transform(&self, m: DMat4) -> Self {
        let mut new_impl = (*self.impl_).clone();
        for v in &mut new_impl.vert_pos {
            let v_h: DVec4 = m * v.extend(1.0);
            *v = DVec3::new(v_h.x, v_h.y, v_h.z);
        }
        Self::new(new_impl)
    }

    pub fn boolean(&self, other: &Self, op: OpType) -> Self {
        let boolean = Boolean3::new(&self.impl_, &other.impl_, op);
        let result_impl = create_boolean_result(
            &self.impl_,
            &other.impl_,
            &boolean.xv12_,
            &boolean.xv21_,
            &boolean.w03_,
            &boolean.w30_,
            op,
        );
        Self::new(result_impl.mesh_r)
    }
}

impl std::ops::Add for Manifold {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self.boolean(&other, OpType::Add)
    }
}

impl std::ops::Sub for Manifold {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self.boolean(&other, OpType::Subtract)
    }
}

impl std::ops::BitXor for Manifold {
    type Output = Self;
    fn bitxor(self, other: Self) -> Self {
        self.boolean(&other, OpType::Intersect)
    }
}
