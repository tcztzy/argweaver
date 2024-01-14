use miette::IntoDiagnostic;

fn main() -> miette::Result<()> {
    pyo3_build_config::add_extension_module_link_args();
    cmake::build(".");
    let path = std::path::PathBuf::from("src");
    let mut b = autocxx_build::Builder::new("src/lib.rs", [&path]).build()?;
    let source_files = glob::glob("src/argweaver/*.cpp")
        .into_diagnostic()?
        .map(|p| p.unwrap());
    b.flag_if_supported("-std=c++14")
        .files(source_files)
        .compile("argweavers");
    println!("cargo:rerun-if-changed=src/lib.rs");
    Ok(())
}
