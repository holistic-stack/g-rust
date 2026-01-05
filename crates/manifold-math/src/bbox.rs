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

use glam::{DMat3, DMat4, DVec2, DVec3, Vec3Swizzles, Vec4Swizzles};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Box {
    pub min: DVec3,
    pub max: DVec3,
}

impl Box {
    pub fn new(p1: DVec3, p2: DVec3) -> Self {
        Self {
            min: p1.min(p2),
            max: p1.max(p2),
        }
    }

    pub fn size(&self) -> DVec3 {
        self.max - self.min
    }

    pub fn center(&self) -> DVec3 {
        0.5 * (self.max + self.min)
    }

    pub fn scale(&self) -> f64 {
        let abs_max = self.min.abs().max(self.max.abs());
        abs_max.max_element()
    }

    pub fn contains(&self, p: DVec3) -> bool {
        p.cmpge(self.min).all() && self.max.cmpge(p).all()
    }

    pub fn contains_box(&self, other: &Box) -> bool {
        other.min.cmpge(self.min).all() && self.max.cmpge(other.max).all()
    }

    pub fn union_box(&self, other: &Box) -> Box {
        Box {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }

    pub fn union_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    pub fn transform(&self, transform: DMat4) -> Box {
        if !self.is_finite() {
            return *self;
        }
        let mut out = Box::default();
        let corners = [
            DVec3::new(self.min.x, self.min.y, self.min.z),
            DVec3::new(self.min.x, self.min.y, self.max.z),
            DVec3::new(self.min.x, self.max.y, self.min.z),
            DVec3::new(self.min.x, self.max.y, self.max.z),
            DVec3::new(self.max.x, self.min.y, self.min.z),
            DVec3::new(self.max.x, self.min.y, self.max.z),
            DVec3::new(self.max.x, self.max.y, self.min.z),
            DVec3::new(self.max.x, self.max.y, self.max.z),
        ];
        for corner in corners {
            let transformed = (transform * corner.extend(1.0)).xyz();
            out.union_point(transformed);
        }
        out
    }

    pub fn does_overlap(&self, other: &Box) -> bool {
        self.min.x <= other.max.x
            && self.min.y <= other.max.y
            && self.min.z <= other.max.z
            && self.max.x >= other.min.x
            && self.max.y >= other.min.y
            && self.max.z >= other.min.z
    }

    pub fn does_overlap_point(&self, p: DVec3) -> bool {
        p.x <= self.max.x && p.x >= self.min.x && p.y <= self.max.y && p.y >= self.min.y
    }

    pub fn is_finite(&self) -> bool {
        self.min.is_finite() && self.max.is_finite()
    }
}

impl Default for Box {
    fn default() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect {
    pub min: DVec2,
    pub max: DVec2,
}

impl Rect {
    pub fn new(a: DVec2, b: DVec2) -> Self {
        Self {
            min: a.min(b),
            max: a.max(b),
        }
    }

    pub fn size(&self) -> DVec2 {
        self.max - self.min
    }

    pub fn area(&self) -> f64 {
        let sz = self.size();
        sz.x * sz.y
    }

    pub fn scale(&self) -> f64 {
        let abs_max = self.min.abs().max(self.max.abs());
        abs_max.max_element()
    }

    pub fn center(&self) -> DVec2 {
        0.5 * (self.max + self.min)
    }

    pub fn contains(&self, p: DVec2) -> bool {
        p.cmpge(self.min).all() && self.max.cmpge(p).all()
    }

    pub fn contains_rect(&self, rect: &Rect) -> bool {
        rect.min.cmpge(self.min).all() && self.max.cmpge(rect.max).all()
    }

    pub fn does_overlap(&self, rect: &Rect) -> bool {
        self.min.x <= rect.max.x
            && self.min.y <= rect.max.y
            && self.max.x >= rect.min.x
            && self.max.y >= rect.min.y
    }

    pub fn is_empty(&self) -> bool {
        self.max.y <= self.min.y || self.max.x <= self.min.x
    }

    pub fn is_finite(&self) -> bool {
        self.min.is_finite() && self.max.is_finite()
    }

    pub fn union_(&mut self, p: DVec2) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    pub fn union_rect(&self, rect: &Rect) -> Rect {
        Rect {
            min: self.min.min(rect.min),
            max: self.max.max(rect.max),
        }
    }

    pub fn transform(&self, m: DMat3) -> Rect {
        let min_t = (m * self.min.extend(1.0)).xy();
        let max_t = (m * self.max.extend(1.0)).xy();
        Rect {
            min: min_t.min(max_t),
            max: min_t.max(max_t),
        }
    }
}

impl Default for Rect {
    fn default() -> Self {
        Self {
            min: DVec2::splat(f64::INFINITY),
            max: DVec2::splat(f64::NEG_INFINITY),
        }
    }
}
