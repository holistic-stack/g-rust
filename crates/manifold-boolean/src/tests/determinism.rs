use glam::DVec3;
use glam::DVec3;
use manifold::Manifold;
use rayon::ThreadPoolBuilder;

use crate::tests::util::{assert_deterministic_runs, assert_same_signature};

fn scenario() -> Manifold {
    let base = Manifold::cube(DVec3::splat(1.0), true);
    let shifted = base.translate(DVec3::new(0.35, -0.2, 0.45));
    let sphere = Manifold::sphere(0.75, 0);
    (base + shifted) - sphere
}

#[test]
fn determinism_multi_thread_runs() {
    assert_deterministic_runs(3, || scenario());
}

#[test]
fn determinism_thread_pool_parity() {
    let multi = scenario();
    let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let single = pool.install(|| scenario());
    assert_same_signature(|| multi, || single);
}
