#![feature(test)]
#[cfg(feature = "extension-module")]
use pyo3::prelude::*;

pub mod ffi {
    use autocxx::prelude::*;
    include_cpp! {
        #include "argweaver/local_tree.h"
        #include "argweaver/sequences.h"
        safety!(unsafe_ffi)
        // local_tree.h
        block!("argweaver::LocalTrees_iterator")
        block!("argweaver::LocalTrees_const_iterator")
        block!("argweaver::LocalTrees_reverse_iterator")
        block!("argweaver::LocalTrees_const_reverse_iterator")
        generate!("argweaver::LocalTrees")
        // sequences.h
        generate!("argweaver::Sites")
        generate!("argweaver::read_sites")
        generate!("argweaver::read_sites1")
    }
    pub use ffi::argweaver::*;
}
pub mod ser;
pub mod sites;
#[cfg(test)]
mod tests;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[cfg(feature = "extension-module")]
/// Ancestral recombination graph sampling method
#[pymodule]
fn s(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<sites::Sites>()?;
    #[cfg(not(test))]
    m.add_function(wrap_pyfunction!(sites::read_sites, m)?)?;
    Ok(())
}
