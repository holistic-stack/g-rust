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

/// Sort the range in ascending order using the comparison function.
/// Port of C++ stable_sort from parallel.h:1028
pub fn stable_sort<T, F>(policy: ExecutionPolicy, data: &mut [T], compare: F)
where
    T: Send + Sync + Copy,
    F: Fn(&T, &T) -> std::cmp::Ordering + Send + Sync + Copy,
{
    match policy {
        ExecutionPolicy::Seq => {
            data.sort_by(compare);
        }
        ExecutionPolicy::Par => {
            let n = data.len();
            if n <= K_SEQ_THRESHOLD {
                data.sort_by(compare);
            } else {
                let mut tmp = data.to_vec();
                merge_sort_rec(&mut tmp, data, compare);
            }
        }
    }
}

fn merge_sort_rec<T, F>(src: &mut [T], dest: &mut [T], compare: F)
where
    T: Send + Sync + Copy,
    F: Fn(&T, &T) -> std::cmp::Ordering + Send + Sync + Copy,
{
    let n = src.len();
    if n <= K_SEQ_THRESHOLD {
        dest.copy_from_slice(src);
        dest.sort_by(compare);
    } else {
        let mid = n / 2;
        let (src_left, src_right) = src.split_at_mut(mid);
        let (dest_left, dest_right) = dest.split_at_mut(mid);
        rayon::join(
            || merge_sort_rec(dest_left, src_left, compare),
            || merge_sort_rec(dest_right, src_right, compare),
        );
        merge_rec(src_left, src_right, dest, compare);
    }
}

fn merge_rec<T, F>(src1: &[T], src2: &[T], dest: &mut [T], compare: F)
where
    T: Send + Sync + Copy,
    F: Fn(&T, &T) -> std::cmp::Ordering + Send + Sync + Copy,
{
    let n1 = src1.len();
    let n2 = src2.len();
    if n1 < n2 {
        merge_rec(src2, src1, dest, compare);
        return;
    }
    if n1 == 0 {
        return;
    }
    if n1 + n2 <= K_SEQ_THRESHOLD {
        let mut i = 0;
        let mut j = 0;
        let mut k = 0;
        while i < n1 && j < n2 {
            if compare(&src1[i], &src2[j]) != std::cmp::Ordering::Greater {
                dest[k] = src1[i];
                i += 1;
            } else {
                dest[k] = src2[j];
                j += 1;
            }
            k += 1;
        }
        if i < n1 {
            dest[k..k + n1 - i].copy_from_slice(&src1[i..n1]);
        }
        if j < n2 {
            dest[k..k + n2 - j].copy_from_slice(&src2[j..n2]);
        }
    } else {
        let q1 = n1 / 2;
        let val = src1[q1];
        let q2 = binary_search_lower_bound(src2, &val, compare);
        let q3 = q1 + q2;
        dest[q3] = val;

        let (dest_left, dest_right_all) = dest.split_at_mut(q3);
        let (_, dest_right) = dest_right_all.split_at_mut(1);

        rayon::join(
            || merge_rec(&src1[..q1], &src2[..q2], dest_left, compare),
            || merge_rec(&src1[q1 + 1..], &src2[q2..], dest_right, compare),
        );
    }
}

fn binary_search_lower_bound<T, F>(data: &[T], val: &T, compare: F) -> usize
where
    F: Fn(&T, &T) -> std::cmp::Ordering,
{
    let mut low = 0;
    let mut high = data.len();
    while low < high {
        let mid = low + (high - low) / 2;
        if compare(&data[mid], val) == std::cmp::Ordering::Less {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low
}
