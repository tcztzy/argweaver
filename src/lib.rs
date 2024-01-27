use pyo3::prelude::*;

pub mod ffi {
    use autocxx::prelude::*;
    include_cpp! {
        #include "argweaver/local_tree.h"
        safety!(unsafe_ffi)
        block!("argweaver::LocalTrees_iterator")
        block!("argweaver::LocalTrees_const_iterator")
        block!("argweaver::LocalTrees_reverse_iterator")
        block!("argweaver::LocalTrees_const_reverse_iterator")
        generate!("argweaver::LocalTrees")
    }
    pub use ffi::argweaver::*;
}
pub mod ser;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// Ancestral recombination graph sampling method
#[pymodule]
fn s(_py: Python, _m: &PyModule) -> PyResult<()> {
    Ok(())
}
