use std::convert::TryInto;
use std::fs::read_to_string;

use autocxx::prelude::{c_int, UniquePtr, WithinUniquePtr};
use nom::{
    bytes::complete::tag,
    character::complete::{alphanumeric1, digit1, newline, one_of, tab},
    combinator::map_res,
    multi::{many1, separated_list1},
    sequence::{delimited, separated_pair, terminated},
    IResult,
};
use polars::{
    chunked_array::ChunkedArray,
    frame::DataFrame,
    series::{IntoSeries, Series},
};
#[cfg(feature = "extension-module")]
use pyo3::{exceptions::PyIndexError, prelude::*, types::PySlice};
#[cfg(feature = "extension-module")]
use pyo3_polars::PyDataFrame;

use crate::{ffi, Result};

#[cfg(not(feature = "extension-module"))]
pub struct Sites {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub data: DataFrame,
}

#[cfg(feature = "extension-module")]
#[pyclass]
pub struct Sites {
    #[pyo3(get, set)]
    pub chrom: String,
    #[pyo3(get, set)]
    pub start: usize,
    #[pyo3(get, set)]
    pub end: usize,
    pub data: DataFrame,
}

#[cfg(feature = "extension-module")]
#[derive(FromPyObject)]
enum SliceOrIsize<'a> {
    Slice(&'a PySlice),
    Isize(isize),
}

#[cfg(feature = "extension-module")]
#[pymethods]
impl Sites {
    #[new]
    fn new(chrom: String, start: usize, end: usize, data: PyDataFrame) -> Self {
        Self {
            chrom,
            start,
            end,
            data: data.0,
        }
    }
    #[getter]
    fn data(&self) -> PyResult<PyDataFrame> {
        Ok(PyDataFrame(self.data.clone()))
    }

    fn __getitem__(&self, idx: SliceOrIsize) -> PyResult<Self> {
        use polars::prelude::{col, lit, IntoLazy};
        let (df, start, end) = match idx {
            SliceOrIsize::Slice(slice) => {
                let indices = slice.indices(self.end as i64)?;
                let start = indices.start;
                let stop = indices.stop;
                if indices.step != 1 {
                    return Err(PyIndexError::new_err("step must be 1"));
                }
                let df = self.data.clone().lazy().filter(
                    col("pos")
                        .gt_eq(lit(start as u32))
                        .and(col("pos").lt_eq(lit(stop as u32))),
                );
                (df.collect().unwrap(), start as usize, stop as usize)
            }
            SliceOrIsize::Isize(_) => {
                return Err(PyIndexError::new_err("not supported yet"));
            }
        };
        Ok(Self {
            chrom: self.chrom.clone(),
            start,
            end,
            data: df,
        })
    }
}

impl TryInto<UniquePtr<ffi::Sites>> for Sites {
    type Error = Box<dyn std::error::Error>;

    fn try_into(self) -> Result<UniquePtr<ffi::Sites>> {
        use std::io::Write;
        use std::os::unix::ffi::OsStrExt;
        let mut s = ffi::Sites::new(
            self.chrom.clone(),
            c_int(self.start as i32),
            c_int(self.end as i32),
        )
        .within_unique_ptr();
        let mut temp_file = tempfile::NamedTempFile::new()?;
        write!(temp_file, "{}", self.to_string())?;
        unsafe {
            ffi::read_sites1(
                temp_file.path().as_os_str().as_bytes().as_ptr() as *const i8,
                std::pin::Pin::<&mut ffi::Sites>::into_inner_unchecked(s.pin_mut()),
                c_int(-1),
                c_int(-1),
                false,
            );
        }
        Ok(s)
    }
}

impl std::string::ToString for Sites {
    fn to_string(&self) -> String {
        let mut s = String::new();
        s.push_str("NAMES\t");
        let names: Vec<&str> = self
            .data
            .get_column_names()
            .iter()
            .filter(|&&c| c != "pos")
            .map(|&c| c)
            .collect();
        s.push_str(names.join("\t").as_str());
        s.push('\n');
        s.push_str(&format!(
            "REGION\t{}\t{}\t{}\n",
            self.chrom, self.start, self.end
        ));
        let pos = self.data.column("pos").unwrap().u32().unwrap();
        let mut alleles = self
            .data
            .columns(names)
            .unwrap()
            .iter()
            .map(|s| Ok(s.u32().unwrap().into_iter()))
            .collect::<Result<Vec<_>>>()
            .unwrap();
        for row in 0..self.data.height() {
            let p = pos.get(row).unwrap();
            s.push_str(&format!("{}\t", p));
            for individual in &mut alleles {
                s.push(std::char::from_u32(individual.next().unwrap().unwrap()).unwrap());
            }
            s.push('\n');
        }
        s
    }
}

fn parse_names(input: &str) -> IResult<&str, Vec<&str>> {
    delimited(
        terminated(tag("NAMES"), tab),
        separated_list1(tab, alphanumeric1),
        newline,
    )(input)
}

#[test]
fn test_parse_names() {
    let input = "NAMES\tA\tB\tC\n";
    let res = parse_names(input);
    assert_eq!(res, Ok(("", vec!["A", "B", "C"])));
}

fn parse_usize(input: &str) -> IResult<&str, usize> {
    map_res(digit1, |s: &str| s.parse::<usize>())(input)
}

fn parse_u32(input: &str) -> IResult<&str, u32> {
    map_res(digit1, |s: &str| s.parse::<u32>())(input)
}

#[test]
fn test_parse_region() {
    let input = "REGION\tchr1\t1\t100\n";
    let res = parse_region(input);
    assert_eq!(res, Ok(("", ("chr1", (1, 100)))));
}

fn parse_region(input: &str) -> IResult<&str, (&str, (usize, usize))> {
    delimited(
        terminated(tag("REGION"), tab),
        separated_pair(
            alphanumeric1,
            tab,
            separated_pair(parse_usize, tab, parse_usize),
        ),
        newline,
    )(input)
}

fn parse_locs(input: &str) -> IResult<&str, (Vec<u32>, Vec<Vec<u32>>)> {
    let (input, res) = separated_list1(
        newline,
        separated_pair(
            parse_u32,
            tab,
            map_res(many1(one_of("ACGT")), |r| {
                Ok::<Vec<_>, nom::error::Error<&str>>(
                    r.into_iter().map(|c| c as u32).collect::<Vec<u32>>(),
                )
            }),
        ),
    )(input)?;
    let (pos, seqs): (Vec<_>, Vec<_>) = res.into_iter().unzip();
    Ok((input, (pos, seqs)))
}

impl Sites {
    pub fn from_path(path: std::path::PathBuf) -> Result<Self> {
        let content = read_to_string(path)?;
        let (input, names) = parse_names(&content).map_err(|e| e.map_input(|s| s.to_owned()))?;
        let num_of_cols = names.len();
        let (input, (chrom, (start, end))) =
            parse_region(&input).map_err(|e| e.map_input(|s| s.to_owned()))?;
        let (input, (pos, rows)) = parse_locs(input).map_err(|e| e.map_input(|s| s.to_owned()))?;
        assert_eq!(input, "\n");
        let pos = ChunkedArray::<polars::datatypes::UInt32Type>::from_vec("pos", pos).into_series();
        for (i, row) in rows.iter().enumerate() {
            assert_eq!(
                row.len(),
                num_of_cols,
                "row {} has {} columns",
                i,
                row.len()
            );
        }
        let samples = names
            .into_iter()
            .enumerate()
            .map(|(i, name)| {
                let col = rows.iter().map(|r| r[i]).collect::<Vec<_>>();
                let series = ChunkedArray::<polars::datatypes::UInt32Type>::from_vec(name, col)
                    .into_series();
                series
            })
            .collect::<Vec<Series>>();
        Ok(Self {
            data: DataFrame::new(std::iter::once(pos).chain(samples).collect())?,
            chrom: chrom.to_string(),
            start,
            end,
        })
    }
}

#[cfg(feature = "extension-module")]
#[pyfunction]
pub fn read_sites(path: std::path::PathBuf) -> PyResult<Sites> {
    let sites = Sites::from_path(std::path::PathBuf::from(path)).unwrap();
    Ok(sites)
}
