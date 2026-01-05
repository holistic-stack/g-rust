use crate::kernel::kernel12;
use crate::kernel::winding03;
use crate::kernel::{Boolean3, Intersections, ManifoldImpl};
use manifold_types::OpType;

/// Boolean3 constructor - main kernel entry point
/// From boolean3.cpp lines 530-553
pub fn boolean3_new(in_p: &ManifoldImpl, in_q: &ManifoldImpl, op: OpType) -> Boolean3 {
    let expand_p = matches!(op, OpType::Add);

    let xv12 = Intersections::new();
    let xv21 = Intersections::new();

    let valid = true;

    let w03 = winding03::winding03_true_true(in_p, in_q, &xv12, expand_p);
    let w30 = winding03::winding03_false_true(in_p, in_q, &xv21, expand_p);

    Boolean3 {
        in_p,
        in_q,
        expand_p,
        xv12_: xv12,
        xv21_: xv21,
        w03_: w03,
        w30_: w30,
        valid,
    }
}
