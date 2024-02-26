use std::os::unix::ffi::OsStrExt;

use argweavers::{ffi, sites::Sites};
use autocxx::prelude::*;
use criterion::{criterion_group, criterion_main, Criterion};

fn ffi_read_sites() {
    let mut s = ffi::Sites::new("", c_int(0), c_int(0)).within_unique_ptr();
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let _ok = unsafe {
        ffi::read_sites1(
            path.join("examples/sim1/sim1.sites")
                .as_mut_os_string()
                .as_bytes()
                .as_ptr() as *const i8,
            std::pin::Pin::<&mut ffi::Sites>::into_inner_unchecked(s.pin_mut()),
            c_int(-1),
            c_int(-1),
            false,
        )
    };
}

fn read_sites() {
    let path = std::path::PathBuf::from("examples/sim1/sim1.sites");
    let _sites = Sites::from_path(&path).unwrap();
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("read_sites");
    group.bench_function("ffi_read_sites", |b| b.iter(ffi_read_sites));
    group.bench_function("read_sites", |b| b.iter(read_sites));
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
