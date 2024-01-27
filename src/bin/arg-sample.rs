use std::path::PathBuf;

use clap::{Args, Parser};

use argweavers::{
    ser::{StatsRecord, StatsStage, StatsWriter},
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

/// Sampler for large ancestral recombination graphs
#[derive(Parser)]
#[clap(author, version, about)]
struct Config {
    #[command(flatten, next_help_heading = "Input/output options")]
    io_args: IOArgs,
    #[command(flatten, next_help_heading = "Sampling")]
    sampling_args: SamplingArgs,
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
}

fn ensure_output_dir_exists<P: AsRef<std::path::Path>>(path: P) -> Result<()> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let config = Config::parse();
    ensure_output_dir_exists(&config.io_args.out_prefix)?;
    config.check_overwrite()?;
    if let Some(fasta_file) = &config.io_args.input_options.fasta_file {
        println!("fasta_file: {:?}", fasta_file);
    }
    let stats_filename = &config.stats_filename();
    let mut writer = StatsWriter::from_path(stats_filename)?;
    writer.serialize(&StatsRecord {
        stage: StatsStage::Resample,
        iter: 0,
        prior: 0.0,
    })?;
    println!("Hello, world!");
    Ok(())
}
