fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(feature = "extension-module")]
    pyo3_build_config::add_extension_module_link_args();
    let path = std::path::PathBuf::from("src");
    let mut b = autocxx_build::Builder::new("src/lib.rs", [&path]).build()?;
    let source_files = glob::glob("src/argweaver/*.cpp")?.map(|x| x.unwrap());
    b.flag_if_supported("-std=c++14")
        .files(source_files)
        .flag("-w")
        .compile("argweavers");
    println!("cargo:rerun-if-changed=src/lib.rs");
    Ok(())
}
