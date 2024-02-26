use std::convert::TryInto;
use std::fs::read_to_string;

use autocxx::prelude::{c_int, UniquePtr, WithinUniquePtr};
use bio::{bio_types::genome::AbstractLocus, io::fasta};
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
use rust_htslib::bcf::{self, Read};

use crate::{de::parse_names, ffi, Result};

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
            map_res(many1(one_of("ACGTN-")), |r| {
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
    pub fn from_path(path: &std::path::PathBuf) -> Result<Self> {
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
    /// Read a sites file from a multi-sequence alignment FASTA file.
    /// All sequences must be the same length.
    pub fn from_msa(path: &std::path::PathBuf) -> Result<Self> {
        let reader = fasta::Reader::from_file(path)?;
        let mut records = reader.records();
        let first = records.next().unwrap()?;
        let seq = first.seq();
        let len = seq.len();
        let mut seqs = vec![Series::from_iter(seq.iter().map(|b| *b as u32)).with_name(first.id())];
        for record in records {
            let record = record?;
            let seq = record.seq();
            assert_eq!(seq.len(), len);
            seqs.push(Series::from_iter(seq.iter().map(|b| *b as u32)).with_name(record.id()));
        }
        let pos = Series::from_iter(1..=len as u32).with_name("pos");
        Ok(Self {
            data: DataFrame::new(std::iter::once(pos).chain(seqs).collect())?,
            chrom: "chr".to_string(),
            start: 1,
            end: len,
        })
    }

    pub fn from_vcf(path: &std::path::PathBuf) -> Result<Self> {
        let mut reader = bcf::Reader::from_path(path)?;
        let header = reader.header();
        let names: Vec<String> = header
            .samples()
            .into_iter()
            .map(|s| std::str::from_utf8(s).unwrap().to_owned())
            .collect();
        let samples = names
            .into_iter()
            .enumerate()
            .map(|(i, name)| {
                let col = reader
                    .records()
                    .into_iter()
                    .map(|r| {
                        let r = r.unwrap();
                        let allele = std::str::from_utf8(r.alleles()[i + 1]).unwrap();
                        assert!(
                            allele.len() == 1,
                            "allele {} has length {}",
                            allele,
                            allele.len()
                        );
                        allele.chars().next().unwrap() as u32
                    })
                    .collect::<Vec<_>>();
                let series = ChunkedArray::<polars::datatypes::UInt32Type>::from_vec(&name, col)
                    .into_series();
                series
            })
            .collect::<Vec<Series>>();
        let pos = ChunkedArray::<polars::datatypes::UInt32Type>::from_vec(
            "pos",
            reader
                .records()
                .into_iter()
                .map(|r| r.unwrap().pos() as u32)
                .collect::<Vec<u32>>(),
        )
        .into_series();
        let chrom = reader
            .records()
            .into_iter()
            .next()
            .unwrap()
            .unwrap()
            .contig()
            .to_string();
        Ok(Self {
            data: DataFrame::new(std::iter::once(pos).chain(samples).collect())?,
            chrom,
            start: 1,
            end: 100000,
        })
    }

    pub fn from_vcfs(paths: &[std::path::PathBuf]) -> Result<Self> {
        let mut sites = Self::from_vcf(&paths[0])?;
        for path in &paths[1..] {
            let other = Self::from_vcf(path)?;
            sites.hstack_mut(&other);
        }
        Ok(sites)
    }

    fn hstack_mut(&mut self, other: &Self) {
        let self_names = self.data.get_column_names();
        let names: Vec<&str> = other
            .data
            .get_column_names()
            .iter()
            .filter(|&&c| !self_names.contains(&c))
            .map(|&c| c)
            .collect();
        let columns: Vec<Series> = other
            .data
            .columns(names)
            .unwrap()
            .into_iter()
            .map(|s| s.to_owned())
            .collect();
        self.data.hstack_mut(&columns).unwrap();
    }
}

#[cfg(feature = "extension-module")]
#[pyfunction]
pub fn read_sites(path: std::path::PathBuf) -> PyResult<Sites> {
    let sites = Sites::from_path(&path).unwrap();
    Ok(sites)
}
