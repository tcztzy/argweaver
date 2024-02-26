use std::fs::read_to_string;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{builder::ArgAction, Args, Parser};
use log::{error, info, warn, LevelFilter};
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config as LogConfig, Root};
use log4rs::filter::threshold::ThresholdFilter;
use polars::prelude::*;

use argweavers::{
    ser::{StatsRecord, StatsStage, StatsWriter},
    sites::Sites,
    Result,
};

#[derive(Args)]
#[group(required = true, multiple = false)]
struct InputOptions {
    /// sequence alignment in sites format
    #[arg(short, long = "sites", value_name = "sites alignment")]
    sites_file: Option<PathBuf>,
    /// sequence alignment in FASTA format
    #[arg(short, long = "fasta", value_name = "fasta alignment")]
    fasta_file: Option<PathBuf>,
    /// sequence alignment in gzipped vcf format.
    /// Must also supply --region in format chr:start-end, and tabix index file
    /// (.vcf.gz.tbi) must also be present. tabix program needs to be installed
    /// and in PATH or path must be provided with --tabix-dir. Assumes samples
    /// are diploid and unphased; each individual will have _1 and _2 appended
    /// to its name for its two haploid lineages. Indel-type variants are
    /// skipped. Any positions not specified in VCF are assumed to be invariant
    #[arg(long = "vcf", value_name = ".vcf.gz file")]
    vcf_file: Option<PathBuf>,
    /// Same as --vcf, but loads all individuals from multiple .vcf.gz files.
    /// The argument should be a file containing list of files to read (one per
    /// line). All files should be aligned to same reference genome
    #[arg(long = "vcf-files", value_name = "vcf_file_list.tx")]
    vcf_list_file: Option<PathBuf>,
}

impl InputOptions {
    fn sites(&self) -> Result<Sites> {
        if let Some(sites_file) = &self.sites_file {
            return Sites::from_path(sites_file);
        }
        if let Some(fasta_file) = &self.fasta_file {
            return Sites::from_msa(fasta_file);
        }
        if let Some(vcf_file) = &self.vcf_file {
            return Sites::from_vcf(vcf_file);
        }
        if let Some(vcf_list_file) = &self.vcf_list_file {
            let vcfs: Vec<PathBuf> = read_to_string(vcf_list_file)?
                .lines()
                .map(PathBuf::from)
                .collect();
            return Sites::from_vcfs(&vcfs);
        }
        Err("No input file provided".into())
    }
}

#[derive(Args)]
struct IOArgs {
    /// prefix for all output filenames
    #[arg(
        short,
        long = "output",
        default_value = "arg-sample",
        value_name = "output prefix"
    )]
    out_prefix: PathBuf,
    #[command(flatten)]
    input_options: InputOptions,

    /// Used to rename sequences (usually from cryptic names in VCF files to
    /// more meaningful names). The file should contain two columns, with the
    /// old sequence name, followed by the new sequence name. Names can be
    /// diploid or haploid, and any old names not matching the current sequence
    /// set will be ignored.
    #[arg(long = "rename-seqs", value_name = "name_map_file.txt")]
    rename_file: Option<PathBuf>,

    /// sample ARG for only a region of the sites (optional). Note the [chr:]
    /// prefix should be added ONLY if alignment is in vcf format. If this
    /// option is given, the vcf must also be indexed with tabix
    #[arg(long = "region", value_name = "[chr:]start-end")]
    region: Option<String>,

    /// file listing NAMES from sites file (or sequences from fasta) to keep;
    /// others will not be used. May be diploid or haploid names (i.e., ind will
    /// have same effect as ind_1 and ind_2).
    #[arg(long, visible_aliases = ["subsites", "keep"], value_delimiter = ' ', num_args = 1..)]
    keep_ids: Option<Vec<String>>,

    /// data is unphased (will integrate over phasings).
    #[arg(long)]
    unphased: bool,

    /// do not gzip output files
    #[clap(long = "no-compress-output", action = ArgAction::SetFalse)]
    compress_output: bool,

    #[clap(
        long = "compress-output",
        overrides_with = "compress_output",
        hide = true
    )]
    _no_compress_output: (),
}

#[derive(Args)]
struct SamplingArgs {
    /// resume a previous run
    #[arg(long = "resume", default_value_t = false)]
    resume: bool,
    /// force an overwrite of a previous run
    #[arg(long = "overwrite", default_value_t = false)]
    overwrite: bool,
}

#[derive(Args)]
struct MiscArgs {
    /// seed for random number generator (default=current time)
    #[arg(short = 'x', long = "randseed", value_name = "random seed")]
    seed: Option<u64>,
}

/// Sampler for large ancestral recombination graphs
#[derive(Parser)]
#[clap(author, version, about)]
struct Config {
    #[command(flatten, next_help_heading = "Input/output options")]
    io_args: IOArgs,
    #[command(flatten, next_help_heading = "Sampling")]
    sampling_args: SamplingArgs,
    #[command(flatten, next_help_heading = "Miscellaneous")]
    misc_args: MiscArgs,
}

impl Config {
    fn stats_filename(&self) -> PathBuf {
        self.io_args.out_prefix.with_extension("stats")
    }
    fn check_overwrite(&self) -> Result<()> {
        let stats_filename = self.stats_filename();
        if !self.sampling_args.resume && !self.sampling_args.overwrite && stats_filename.exists() {
            return Err(format!(
                "stats file already exists: {}\nTo force overwrite use --overwrite option",
                stats_filename.display()
            )
            .into());
        }
        Ok(())
    }
    fn out_postfix(&self) -> &str {
        if self.io_args.compress_output {
            ".gz"
        } else {
            ""
        }
    }
}

fn ensure_output_dir_exists<P: AsRef<std::path::Path>>(path: P) -> Result<()> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    Ok(())
}

fn setup_logging(log_file: &Path, resume: bool) -> Result<()> {
    let stdout = ConsoleAppender::builder().build();

    let logfile = FileAppender::builder().append(resume).build(log_file)?;

    let config = LogConfig::builder()
        .appender(Appender::builder().build("stdout", Box::new(stdout)))
        .appender(
            Appender::builder()
                .filter(Box::new(ThresholdFilter::new(LevelFilter::Debug)))
                .build("logfile", Box::new(logfile)),
        )
        .build(
            Root::builder()
                .appender("stdout")
                .appender("logfile")
                .build(LevelFilter::Trace),
        )?;

    let _handle = log4rs::init_config(config)?;
    Ok(())
}

fn find_previous_smc_file(out_prefix: &PathBuf) -> Result<(PathBuf, i64)> {
    let stats_filename = out_prefix.with_extension("stats");
    info!(
        "Checking previous run from stats file: {}",
        stats_filename.display()
    );
    let stats = CsvReader::from_path(stats_filename)?
        .has_header(true)
        .with_separator(b'\t')
        .finish()
        .map_err(|e| match e {
            polars::error::PolarsError::NoData(_) => {
                polars::error::PolarsError::NoData("stats file is empty".into())
            }
            _ => e,
        })?;
    let stage = stats.column("stage")?;
    let resample = stats
        .filter(&stage.equal("resample")?)?
        .sort(["iter"], true, true)?;
    for i in resample.column("iter")?.i64()? {
        let i = i.unwrap();
        let smc_file = out_prefix.with_extension(format!("{}.smc.gz", i));
        if smc_file.exists() {
            return Ok((smc_file, i));
        }
        let smc_file = smc_file.with_extension("");
        if smc_file.exists() {
            return Ok((smc_file, i));
        }
    }
    let msg = "Could not find any previously written SMC files. Try disabling resume";
    error!("{}", msg);
    Err(std::io::Error::new(std::io::ErrorKind::NotFound, msg).into())
}

fn setup_resume(config: &mut Config) -> Result<()> {
    info!("Resuming from previous run");
    let (arg_file, resume_iter) = find_previous_smc_file(&config.io_args.out_prefix)?;
    let sites_file = config.io_args.out_prefix.with_extension(format!(
        "{}.sites{}",
        resume_iter,
        config.out_postfix(),
    ));
    if sites_file.exists() {
        let test_sites = Sites::from_path(&sites_file)?;
        info!(
            "Detected phased output sites file. Using %s as input {} and assuming data is unphased",
            &sites_file.display()
        );
        config.io_args.input_options.sites_file = Some(sites_file);
        config.io_args.unphased = true;
        let mask_file = config
            .io_args
            .out_prefix
            .with_extension("masked_regions.bed");
    }
    let sites = config.io_args.input_options.sites()?;
    info!(
        "resuming at stage={}, iter={}, arg={}",
        "resample",
        resume_iter,
        arg_file.display()
    );
    Ok(())
}

fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
}

fn main() -> Result<()> {
    let mut config = Config::parse();
    ensure_output_dir_exists(&config.io_args.out_prefix)?;
    config.check_overwrite()?;
    if let Some(fasta_file) = &config.io_args.input_options.fasta_file {
        println!("fasta_file: {:?}", fasta_file);
    }
    let stats_filename = &config.stats_filename();
    setup_logging(
        &config.io_args.out_prefix.with_extension("log"),
        config.sampling_args.resume,
    )?;
    if config.sampling_args.resume {
        match setup_resume(&mut config) {
            Ok(_) => info!("RESUME"),
            Err(e) => {
                if config.sampling_args.overwrite {
                    warn!("Resume failed.  Sampling will start from scratch since overwrite is enabled.");
                    config.sampling_args.resume = false;
                } else {
                    error!("Could not resume: {}", e);
                    return Err(e);
                }
            }
        }
    }
    info!("arg-sample {}", env!("CARGO_PKG_VERSION"));
    info!(
        "next line for backward compatibility\ncommand: {}",
        std::env::args()
            .into_iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    );
    let seed = config.misc_args.seed.unwrap_or(now());
    info!("random seed: {}", seed);
    unsafe {
        libc::srand(seed as _);
    }
    if let Some(keep_ids) = config.io_args.keep_ids {
        if keep_ids.len() == 1 && PathBuf::from(&keep_ids[0]).exists() {
            let keep_ids = std::fs::read_to_string(&keep_ids[0])?;
            config.io_args.keep_ids =
                Some(keep_ids.split_whitespace().map(|s| s.to_string()).collect());
        }
    }
    let mut writer = StatsWriter::from_path(stats_filename)?;
    writer.serialize(&StatsRecord {
        stage: StatsStage::Resample,
        iter: 0,
        prior: 0.0,
    })?;
    writer.serialize(&StatsRecord {
        stage: StatsStage::Resample,
        iter: 1,
        prior: 0.0,
    })?;
    Ok(())
}
