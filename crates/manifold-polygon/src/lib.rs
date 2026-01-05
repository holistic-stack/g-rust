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
pub mod polygon;
pub mod triangulation;

pub use bbox::Rect;
pub use polygon::triangulate_convex;
pub use polygon::{Polygons, PolygonsIdx, SimplePolygon, SimplePolygonIdx};
pub use triangulation::{ccw, triangulate_idx};
