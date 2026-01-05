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

use crate::Box;
use glam::DVec3;

pub const K_NO_CODE: u32 = 0xFFFFFFFF;

pub fn spread_bits_3(mut v: u32) -> u32 {
    v = 0xFF0000FF & (v.wrapping_mul(0x00010001));
    v = 0x0F00F00F & (v.wrapping_mul(0x00000101));
    v = 0xC30C30C3 & (v.wrapping_mul(0x00000011));
    v = 0x49249249 & (v.wrapping_mul(0x00000005));
    v
}

pub fn morton_code(position: DVec3, bbox: Box) -> u32 {
    let xyz = (position - bbox.min) / (bbox.max - bbox.min);
    let x = f64::min(1023.0, f64::max(0.0, 1024.0 * xyz.x)) as u32;
    let y = f64::min(1023.0, f64::max(0.0, 1024.0 * xyz.y)) as u32;
    let z = f64::min(1023.0, f64::max(0.0, 1024.0 * xyz.z)) as u32;
    spread_bits_3(x) * 4 + spread_bits_3(y) * 2 + spread_bits_3(z)
}
