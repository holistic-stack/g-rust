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

use crate::policy::ExecutionPolicy;
use rayon::prelude::*;

pub fn for_each<I, F, T>(policy: ExecutionPolicy, iter: I, f: F)
where
    I: IntoParallelIterator<Item = T> + IntoIterator<Item = T>,
    F: Fn(T) + Sync + Send + Copy,
    T: Send,
{
    match policy {
        ExecutionPolicy::Par => iter.into_par_iter().for_each(f),
        ExecutionPolicy::Seq => iter.into_iter().for_each(f),
    }
}

pub fn for_each_n<F>(policy: ExecutionPolicy, start: usize, n: usize, f: F)
where
    F: Fn(usize) + Sync + Send + Copy,
{
    match policy {
        ExecutionPolicy::Par => (start..start + n).into_par_iter().for_each(f),
        ExecutionPolicy::Seq => (start..start + n).into_iter().for_each(f),
    }
}

pub fn transform<I, O, F, T>(policy: ExecutionPolicy, input: I, output: &mut [O], f: F)
where
    I: IntoParallelIterator<Item = T> + IntoIterator<Item = T>,
    I::Iter: IndexedParallelIterator,
    F: Fn(T) -> O + Sync + Send + Copy,
    O: Send,
    T: Send,
{
    match policy {
        ExecutionPolicy::Par => {
            input
                .into_par_iter()
                .zip(output.par_iter_mut())
                .for_each(|(in_val, out_val)| {
                    *out_val = f(in_val);
                });
        }
        ExecutionPolicy::Seq => {
            input
                .into_iter()
                .zip(output.iter_mut())
                .for_each(|(in_val, out_val)| {
                    *out_val = f(in_val);
                });
        }
    }
}

pub fn reduce<T, I, F, R>(policy: ExecutionPolicy, iter: I, identity: T, f: F, reduce_op: R) -> T
where
    I: IntoParallelIterator<Item = T> + IntoIterator<Item = T>,
    F: Fn(T, T) -> T + Sync + Send + Copy,
    R: Fn(T, T) -> T + Sync + Send + Copy,
    T: Send + Sync + Copy,
{
    match policy {
        ExecutionPolicy::Par => iter.into_par_iter().reduce(|| identity, reduce_op),
        ExecutionPolicy::Seq => iter.into_iter().fold(identity, f),
    }
}

pub fn copy_if<T, P>(policy: ExecutionPolicy, input: &[T], output: &mut Vec<T>, pred: P)
where
    T: Send + Sync + Copy,
    P: Fn(&T) -> bool + Sync + Send + Copy,
{
    match policy {
        ExecutionPolicy::Seq => {
            for x in input {
                if pred(x) {
                    output.push(*x);
                }
            }
        }
        ExecutionPolicy::Par => {
            let mut result: Vec<T> = input.into_par_iter().filter(|x| pred(x)).copied().collect();
            output.append(&mut result);
        }
    }
}
