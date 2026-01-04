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

use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum Error {
    #[error("No Error")]
    NoError,
    #[error("Non Finite Vertex")]
    NonFiniteVertex,
    #[error("Not Manifold")]
    NotManifold,
    #[error("Vertex Out Of Bounds")]
    VertexOutOfBounds,
    #[error("Properties Wrong Length")]
    PropertiesWrongLength,
    #[error("Missing Position Properties")]
    MissingPositionProperties,
    #[error("Merge Vectors Different Lengths")]
    MergeVectorsDifferentLengths,
    #[error("Merge Index Out Of Bounds")]
    MergeIndexOutOfBounds,
    #[error("Transform Wrong Length")]
    TransformWrongLength,
    #[error("Run Index Wrong Length")]
    RunIndexWrongLength,
    #[error("Face ID Wrong Length")]
    FaceIDWrongLength,
    #[error("Invalid Construction")]
    InvalidConstruction,
    #[error("Result Too Large")]
    ResultTooLarge,
}
