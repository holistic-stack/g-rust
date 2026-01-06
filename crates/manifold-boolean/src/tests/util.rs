use glam::DVec3;
use manifold::Manifold;

fn to_bits(v: &DVec3) -> [u64; 3] {
    [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
}

pub fn canonical_signature(manifold: &Manifold) -> (Vec<[u64; 3]>, Vec<[i32; 3]>) {
    let (verts, tris) = manifold.mesh_data();
    let mut indexed: Vec<(usize, DVec3)> = verts.into_iter().enumerate().collect();
    indexed.sort_by(|a, b| {
        a.1.x
            .total_cmp(&b.1.x)
            .then(a.1.y.total_cmp(&b.1.y))
            .then(a.1.z.total_cmp(&b.1.z))
            .then(a.0.cmp(&b.0))
    });

    let mut old_to_new = vec![0i32; indexed.len()];
    let mut ordered_bits = Vec::with_capacity(indexed.len());
    for (new_idx, (old_idx, v)) in indexed.into_iter().enumerate() {
        old_to_new[old_idx] = new_idx as i32;
        ordered_bits.push(to_bits(&v));
    }

    let mut mapped_tris: Vec<[i32; 3]> = tris
        .into_iter()
        .map(|t| {
            [
                old_to_new[t[0] as usize],
                old_to_new[t[1] as usize],
                old_to_new[t[2] as usize],
            ]
        })
        .collect();
    mapped_tris.sort_by(|a, b| a.cmp(b));

    (ordered_bits, mapped_tris)
}

pub fn counts(manifold: &Manifold) -> (usize, usize) {
    let (verts, tris) = manifold.mesh_data();
    (verts.len(), tris.len())
}

fn hash_u64(h: u64, v: u64) -> u64 {
    h.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(v)
}

pub fn signature_hash(manifold: &Manifold) -> u64 {
    let (verts, tris) = canonical_signature(manifold);
    let mut h = 0u64;
    for v in verts {
        h = hash_u64(h, v[0]);
        h = hash_u64(h, v[1]);
        h = hash_u64(h, v[2]);
    }
    for t in tris {
        h = hash_u64(h, t[0] as u64);
        h = hash_u64(h, t[1] as u64);
        h = hash_u64(h, t[2] as u64);
    }
    h
}

pub fn assert_deterministic_runs<F>(runs: usize, f: F)
where
    F: Fn() -> Manifold,
{
    let mut first: Option<u64> = None;
    for _ in 0..runs {
        let m = f();
        let h = signature_hash(&m);
        if let Some(prev) = first {
            assert_eq!(h, prev);
        } else {
            first = Some(h);
        }
    }
}

pub fn assert_same_signature<F, G>(a: F, b: G)
where
    F: FnOnce() -> Manifold,
    G: FnOnce() -> Manifold,
{
    let ha = signature_hash(&a());
    let hb = signature_hash(&b());
    assert_eq!(ha, hb);
}
