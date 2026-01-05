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

use glam::DVec2;

/// 2D axis-aligned bounding box
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Rect {
    pub min: DVec2,
    pub max: DVec2,
}

impl Rect {
    pub fn new(min: DVec2, max: DVec2) -> Self {
        Self {
            min: min.min(max),
            max: max.max(min),
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

    pub fn union_point(&mut self, p: DVec2) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    pub fn union_rect(&self, rect: &Rect) -> Rect {
        Rect {
            min: self.min.min(rect.min),
            max: self.max.max(rect.max),
        }
    }
}
