#[cfg(feature = "extension-module")]
use pyo3::prelude::*;

pub mod binaries;
pub mod de;
pub mod ffi {
    use autocxx::prelude::*;
    include_cpp! {
        #include "argweaver/ffi.h"
        #include "argweaver/local_tree.h"
        #include "argweaver/model.h"
        #include "argweaver/sequences.h"
        safety!(unsafe_ffi)
        // ffi.h
        generate!("empty_seqnames")
        generate!("vector_push_string")
        generate!("empty_times")
        generate!("vector_push_double")
        generate!("read_local_trees")
        generate!("write_local_trees_as_bed")
        // local_tree.h
        block!("argweaver::LocalTrees_iterator")
        block!("argweaver::LocalTrees_const_iterator")
        block!("argweaver::LocalTrees_reverse_iterator")
        block!("argweaver::LocalTrees_const_reverse_iterator")
        generate!("argweaver::LocalTrees")
        // model.h
        generate!("argweaver::ArgModel")
        // sequences.h
        generate!("argweaver::Sites")
        generate!("argweaver::read_sites")
        generate!("argweaver::read_sites1")
    }
    pub use ffi::{argweaver::*, spidir::*, *};
}
pub mod fs;
pub mod io;
pub mod ser;
pub mod sites;
#[cfg(test)]
mod tests;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[cfg(feature = "extension-module")]
#[pyfunction]
fn smc2bed(args: Option<Vec<String>>) -> PyResult<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    binaries::smc2bed(Some(args.unwrap_or(argv))).unwrap();
    Ok(())
}

#[cfg(feature = "extension-module")]
/// Ancestral recombination graph sampling method
#[pymodule]
fn argweavers(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<sites::Sites>()?;
    m.add_function(wrap_pyfunction!(sites::read_sites, m)?)?;
    m.add_function(wrap_pyfunction!(smc2bed, m)?)?;
    Ok(())
}
