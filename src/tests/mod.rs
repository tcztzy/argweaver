extern crate test;

use std::fs::read_to_string;
use std::os::unix::ffi::OsStrExt;

use crate::ffi;
use crate::sites::Sites;
use autocxx::prelude::*;

#[test]
fn test_ffi_read_sites() {
    let mut s = ffi::Sites::new("", c_int(0), c_int(0)).within_unique_ptr();
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let ok = unsafe {
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
    assert!(ok);
}

#[bench]
fn bench_ffi_read_sites(b: &mut test::Bencher) {
    b.iter(test_ffi_read_sites);
}
#[test]
fn test_from_path() {
    let path = std::path::PathBuf::from("examples/sim1/sim1.sites");
    let sites = Sites::from_path(path).unwrap();
    assert_eq!(sites.chrom, "chr");
    assert_eq!(sites.start, 1);
    assert_eq!(sites.end, 100000);
    assert_eq!(sites.data.shape(), (170, 9)); // 170 sites, 8 samples + pos
}

#[test]
fn test_from_path_and_to_string() {
    let path = std::path::PathBuf::from("examples/sim1/sim1.sites");
    let sites = Sites::from_path(path).unwrap();
    let s = sites.to_string();
    assert_eq!(
        s.to_string(),
        read_to_string("examples/sim1/sim1.sites").unwrap()
    );
}

#[test]
fn test_try_into_ffi_sites() {
    let path = std::path::PathBuf::from("examples/sim1/sim1.sites");
    let sites = Sites::from_path(path).unwrap();
    let ffi_sites: UniquePtr<ffi::Sites> = sites.try_into().unwrap();
    let num_seq: i32 = ffi_sites.get_num_seqs().into();
    assert_eq!(num_seq, 8);
    let num_sites: i32 = ffi_sites.get_num_sites().into();
    assert_eq!(num_sites, 170);
}

#[bench]
fn bench_from_path(b: &mut test::Bencher) {
    b.iter(test_from_path);
}
