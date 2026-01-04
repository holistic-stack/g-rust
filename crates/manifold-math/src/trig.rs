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

use crate::constants::K_PI;

pub fn radians(a: f64) -> f64 {
    a * K_PI / 180.0
}

pub fn degrees(a: f64) -> f64 {
    a * 180.0 / K_PI
}

pub fn sind(x: f64) -> f64 {
    if !x.is_finite() {
        return x.sin();
    }
    if x < 0.0 {
        return -sind(-x);
    }
    let quo = (x / 90.0).round() as i32;
    let x_rem = x - (quo as f64) * 90.0;
    match quo % 4 {
        0 => radians(x_rem).sin(),
        1 => radians(x_rem).cos(),
        2 => -radians(x_rem).sin(),
        3 => -radians(x_rem).cos(),
        _ => 0.0,
    }
}

pub fn cosd(x: f64) -> f64 {
    sind(x + 90.0)
}
