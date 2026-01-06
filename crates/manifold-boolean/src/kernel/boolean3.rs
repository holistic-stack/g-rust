use crate::kernel::kernel12;
use crate::kernel::winding03;
use crate::kernel::{Boolean3, Intersections, ManifoldImpl};
use manifold_types::OpType;

/// Boolean3 constructor - main kernel entry point
/// From boolean3.cpp lines 530-553
pub fn boolean3_new(in_p: &ManifoldImpl, in_q: &ManifoldImpl, op: OpType) -> Boolean3 {
    let expand_p = matches!(op, OpType::Add);

    if in_p.is_empty() || in_q.is_empty() || !in_p.b_box.does_overlap(in_q.b_box) {
        return Boolean3 {
            in_p,
            in_q,
            expand_p,
            xv12_: Intersections::new(),
            xv21_: Intersections::new(),
            w03_: vec![0; in_p.num_vert()],
            w30_: vec![0; in_q.num_vert()],
            valid: true,
        };
    }

    let xv12 = kernel12::intersect12_true_true(in_p, in_q, expand_p);
    let xv21 = kernel12::intersect12_false_true(in_p, in_q, expand_p);

    let int_max_sz = i32::MAX as usize;
    if xv12.x12.len() > int_max_sz || xv21.x12.len() > int_max_sz {
        return Boolean3 {
            in_p,
            in_q,
            expand_p,
            xv12_: xv12,
            xv21_: xv21,
            w03_: Vec::new(),
            w30_: Vec::new(),
            valid: false,
        };
    }

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
        valid: true,
    }
}
