// Copyright 2025 The Manifold Authors.
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

use crate::bbox::Rect;
use crate::polygon::PolyVert;
use manifold_parallel::{stable_sort, ExecutionPolicy};

fn build_two_d_tree_impl(points: &mut [PolyVert], sort_x: bool) {
    if sort_x {
        stable_sort(ExecutionPolicy::Seq, points, |a, b| a.pos.x.partial_cmp(&b.pos.x).unwrap());
    } else {
        stable_sort(ExecutionPolicy::Seq, points, |a, b| a.pos.y.partial_cmp(&b.pos.y).unwrap());
    }

    if points.len() < 2 {
        return;
    }

    let mid = points.len() / 2;
    // partition at mid.
    let (left, rest) = points.split_at_mut(mid);
    let (_pivot, right) = rest.split_at_mut(1);

    build_two_d_tree_impl(left, !sort_x);
    build_two_d_tree_impl(right, !sort_x);
}

pub fn build_two_d_tree(points: &mut [PolyVert]) {
    if points.len() <= 8 {
        return;
    }
    build_two_d_tree_impl(points, true);
}

pub fn query_two_d_tree<F>(points: &[PolyVert], r: Rect, mut f: F)
where
    F: FnMut(&PolyVert),
{
    if points.len() <= 8 {
        for p in points {
            if r.contains(p.pos) {
                f(p);
            }
        }
        return;
    }

    let mut current = Rect {
        min: glam::DVec2::splat(f64::NEG_INFINITY),
        max: glam::DVec2::splat(f64::INFINITY),
    };

    struct StackItem {
        level: i32,
        view_start: usize,
        view_len: usize,
        rect: Rect,
    }

    let mut stack: Vec<StackItem> = Vec::with_capacity(64);

    let mut level = 0;
    let mut view_start = 0;
    let mut view_len = points.len();

    loop {
        let view = &points[view_start..view_start + view_len];
        if view.len() <= 8 {
            for p in view {
                if r.contains(p.pos) {
                    f(p);
                }
            }
            if let Some(item) = stack.pop() {
                level = item.level;
                view_start = item.view_start;
                view_len = item.view_len;
                current = item.rect;
                continue;
            } else {
                break;
            }
        }

        let mut left = current;
        let mut right = current;
        let mid_idx = view.len() / 2;
        let middle = &view[mid_idx];

        if level % 2 == 0 {
            left.max.x = middle.pos.x;
            right.min.x = middle.pos.x;
        } else {
            left.max.y = middle.pos.y;
            right.min.y = middle.pos.y;
        }

        if r.contains(middle.pos) {
            f(middle);
        }

        let process_left = left.does_overlap(&r);
        let process_right = right.does_overlap(&r);

        if process_left {
            if process_right {
                stack.push(StackItem {
                    level: level + 1,
                    view_start: view_start + mid_idx + 1,
                    view_len: view.len() - (mid_idx + 1),
                    rect: right,
                });
            }
            current = left;
            view_len = mid_idx;
            level += 1;
        } else if process_right {
            current = right;
            view_start = view_start + mid_idx + 1;
            view_len = view.len() - (mid_idx + 1);
            level += 1;
        } else {
            if let Some(item) = stack.pop() {
                level = item.level;
                view_start = item.view_start;
                view_len = item.view_len;
                current = item.rect;
            } else {
                break;
            }
        }
    }
}
