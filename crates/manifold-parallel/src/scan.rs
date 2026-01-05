// Copyright 2022 The Manifold Authors.
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

use crate::policy::{ExecutionPolicy, K_SEQ_THRESHOLD};
use rayon::prelude::*;

/// Compute the exclusive prefix sum for the range.
/// Port of C++ exclusive_scan from parallel.h:550
pub fn exclusive_scan<T, F>(
    policy: ExecutionPolicy,
    input: &[T],
    output: &mut [T],
    init: T,
    op: F,
    identity: T,
) where
    T: Send + Sync + Copy,
    F: Fn(T, T) -> T + Send + Sync + Copy,
{
    assert_eq!(input.len(), output.len());
    if input.is_empty() {
        return;
    }

    match policy {
        ExecutionPolicy::Seq => {
            let mut curr = init;
            for i in 0..input.len() {
                output[i] = curr;
                curr = op(curr, input[i]);
            }
        }
        ExecutionPolicy::Par => {
            let n = input.len();
            if n <= K_SEQ_THRESHOLD {
                exclusive_scan(ExecutionPolicy::Seq, input, output, init, op, identity);
                return;
            }

            let num_threads = rayon::current_num_threads();
            let block_size = (n + num_threads - 1) / num_threads;

            let mut block_sums: Vec<T> = (0..num_threads)
                .into_par_iter()
                .map(|i| {
                    let start = i * block_size;
                    if start >= n {
                        return identity;
                    }
                    let end = std::cmp::min(start + block_size, n);
                    let mut sum = identity;
                    for j in start..end {
                        sum = op(sum, input[j]);
                    }
                    sum
                })
                .collect();

            let mut curr = init;
            for i in 0..num_threads {
                let tmp = block_sums[i];
                block_sums[i] = curr;
                curr = op(curr, tmp);
            }

            output
                .par_chunks_mut(block_size)
                .enumerate()
                .for_each(|(i, chunk)| {
                    let start = i * block_size;
                    let mut curr = block_sums[i];
                    for j in 0..chunk.len() {
                        chunk[j] = curr;
                        curr = op(curr, input[start + j]);
                    }
                });
        }
    }
}

/// Compute the inclusive prefix sum for the range.
/// Port of C++ inclusive_scan from parallel.h:486
pub fn inclusive_scan<T, F>(policy: ExecutionPolicy, input: &[T], output: &mut [T], op: F)
where
    T: Send + Sync + Copy,
    F: Fn(T, T) -> T + Send + Sync + Copy,
{
    assert_eq!(input.len(), output.len());
    if input.is_empty() {
        return;
    }

    match policy {
        ExecutionPolicy::Seq => {
            let mut curr = input[0];
            output[0] = curr;
            for i in 1..input.len() {
                curr = op(curr, input[i]);
                output[i] = curr;
            }
        }
        ExecutionPolicy::Par => {
            let n = input.len();
            if n <= K_SEQ_THRESHOLD {
                inclusive_scan(ExecutionPolicy::Seq, input, output, op);
                return;
            }

            let num_threads = rayon::current_num_threads();
            let block_size = (n + num_threads - 1) / num_threads;

            let block_sums: Vec<T> = (0..num_threads)
                .into_par_iter()
                .map(|i| {
                    let start = i * block_size;
                    if start >= n {
                        return input[0];
                    }
                    let end = std::cmp::min(start + block_size, n);
                    let mut sum = input[start];
                    for j in (start + 1)..end {
                        sum = op(sum, input[j]);
                    }
                    sum
                })
                .collect();

            let mut block_offsets = Vec::with_capacity(num_threads);
            let mut curr_offset = None;
            for i in 0..num_threads {
                block_offsets.push(curr_offset);
                let s = block_sums[i];
                match curr_offset {
                    None => curr_offset = Some(s),
                    Some(prev) => curr_offset = Some(op(prev, s)),
                }
            }

            output
                .par_chunks_mut(block_size)
                .enumerate()
                .for_each(|(i, chunk)| {
                    let start = i * block_size;
                    let mut curr = match block_offsets[i] {
                        None => input[start],
                        Some(offset) => op(offset, input[start]),
                    };
                    chunk[0] = curr;
                    for j in 1..chunk.len() {
                        curr = op(curr, input[start + j]);
                        chunk[j] = curr;
                    }
                });
        }
    }
}
