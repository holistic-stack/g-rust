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

pub mod bbox;
pub mod constants;
pub mod interp;
pub mod morton;
pub mod tri_dist;
pub mod trig;
pub mod utils;

pub use bbox::{Box, Rect};
pub use constants::*;
pub use interp::smoothstep;
pub use morton::*;
pub use tri_dist::{distance_triangle_triangle_squared, edge_edge_dist, EdgeEdgeDist};
pub use trig::*;
pub use utils::*;

pub use utils::*;
