use std::path::PathBuf;

use crate::{de::parse_names, ffi, fs::read_to_string, Result};
use autocxx::prelude::*;
use clap::{arg, command, value_parser};
use nom::{
    bytes::complete::tag, multi::separated_list1, number::complete::double, sequence::delimited,
};

fn get_model_times(log_file: &PathBuf) -> Result<Vec<f64>> {
    let input = read_to_string(log_file)?;
    for line in input.lines() {
        if line.starts_with("  times") {
            let (_, times) = delimited(
                tag::<&str, &str, nom::error::Error<&str>>("  times = ["),
                separated_list1(tag(","), double),
                tag("]"),
            )(line)
            .unwrap();
            return Ok(times);
        }
    }
    Ok(Vec::new())
}

pub fn smc2bed(args: Option<Vec<String>>) -> Result<()> {
    let cmd = command!()
        .arg(arg!(
            ["smc-file"] "smc-file can be gzipped"
        ).required(true)
        .value_parser(value_parser!(PathBuf)
        ))
        .arg(arg!(
            --region <region> "Process only these coordinates (1-based)"
        ))
        .arg(arg!(
            --sample <sample> "Give the sample number for this file; this is important when combining multiple smc files."
        ).value_parser(value_parser!(i32)).default_value("0"))
        .arg(arg!(
            --"log-file" <log_file> "Log file from arg-sample run; this is used as input to read model parameters. If not provided, smc2bed will look for log file in directory with smc file."
        ).value_parser(value_parser!(PathBuf)));
    let matches = if let Some(args) = args {
        cmd.get_matches_from(args)
    } else {
        cmd.get_matches()
    };
    let smc_file = matches.get_one::<PathBuf>("smc-file").unwrap();
    if matches.get_one::<String>("region").is_some() {
        return Err(
            "Would not implement region option. You can partition BED file after this command"
                .into(),
        );
    };
    let log_file = match matches.get_one::<PathBuf>("log-file") {
        Some(log_file) => log_file.to_owned(),
        None => {
            let log_file = smc_file
                .with_extension("")
                .with_extension("")
                .with_extension("log");
            if !log_file.exists() {
                return Err(format!(
                    "Could not guess log file name, provide with -l {}",
                    log_file.display()
                )
                .into());
            }
            log_file
        }
    };
    let model = unsafe {
        ffi::ArgModel::new6(log_file.to_str().unwrap().as_ptr() as *const i8).within_unique_ptr()
    };
    let model_times = get_model_times(&log_file)?;
    let mut seqnames = ffi::empty_seqnames();
    let cntnt: String = read_to_string(&smc_file)?;
    let (_, names) = parse_names(&cntnt).map_err(|e| e.map_input(|s| s.to_owned()))?;
    for name in names {
        ffi::vector_push_string(seqnames.pin_mut(), ffi::make_string(name).pin_mut());
    }
    let mut times = ffi::empty_times();
    for time in model_times {
        ffi::vector_push_double(times.pin_mut(), time);
    }
    let trees = ffi::read_local_trees(smc_file.to_str().unwrap(), times.pin_mut());
    ffi::write_local_trees_as_bed(
        trees,
        seqnames,
        model,
        c_int(matches.get_one::<i32>("sample").unwrap().to_owned()),
    );
    Ok(())
}
