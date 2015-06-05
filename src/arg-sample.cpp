// C/C++ includes
#ifdef ARGWEAVER_MPI
#include "mpi.h"
#endif
#include <time.h>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>

// arghmm includes
#include "compress.h"
#include "ConfigParam.h"
#include "emit.h"
#include "fs.h"
#include "logging.h"
#include "mem.h"
#include "parsing.h"
#include "sample_arg.h"
#include "sequences.h"
#include "total_prob.h"
#include "track.h"
#include "est_popsize.h"
#include "mcmcmc.h"


using namespace argweaver;

// version info
#define VERSION_TEXT "0.8"
#define VERSION_INFO  "\
ARGweaver " VERSION_TEXT " \n\
Matt Rasmussen\n\
Sampler for large ancestral recombination graphs\n\
"

// file extensions
const char *SMC_SUFFIX = ".smc";
const char *SITES_SUFFIX = ".sites";
const char *STATS_SUFFIX = ".stats";
const char *LOG_SUFFIX = ".log";


// debug options level
const int DEBUG_OPT = 1;


const int EXIT_ERROR = 1;



// parsing command-line options
class Config
{
public:

    Config()
    {
        make_parser();

        resample_region[0] = -1;
        resample_region[1] = -1;
    }

    void make_parser()
    {
        config.clear();

        // input/output
        config.add(new ConfigParam<string>
                   ("-s", "--sites", "<sites alignment>", &sites_file,
                    "sequence alignment in sites format"));
        config.add(new ConfigParam<string>
                   ("-f", "--fasta", "<fasta alignment>", &fasta_file,
                    "sequence alignment in FASTA format"));
        config.add(new ConfigParam<string>
                   ("-o", "--output", "<output prefix>", &out_prefix,
                    "arg-sample",
                    "prefix for all output filenames (default='arg-sample')"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<SMC file>", &arg_file, "",
                    "initial ARG file (*.smc) for resampling (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--region", "<start>-<end>",
                    &subregion_str, "",
                    "sample ARG for only a region of the sites (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--maskmap", "<sites mask>",
                    &maskmap, "",
                    "mask map file (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--subsites", "<subsites file>", &subsites_file,
                    "file listing NAMES from sites file (or sequences from"
                    " fasta) to keep; others will not be used"));
#ifdef ARGWEAVER_MPI
        config.add(new ConfigSwitch
                   ("", "--mpi", &mpi, "this is an mpi run, add <rank>.sites"
                    " to sites file name and <rank>. to out root, and"
                    " <rank>.smc.gz to --arg option (if given)"));
#endif

        // model parameters
        config.add(new ConfigParamComment("Model parameters"));
        config.add(new ConfigParam<string>
                   ("-N", "--popsize", "<population size>", &popsize_str,
                    "10000",
                    "effective population size (default=1e4)"));
        config.add(new ConfigSwitch
                   ("", "--sample-popsize", &sample_popsize,
                    "sample population size for each time interval using"
                    " Metropolis-Hastings update"));
        config.add(new ConfigSwitch
                   ("", "--sample-popsize-recomb", &sample_popsize_recomb,
                    "do not integrate over recombination events when sampling"
                    "popsize", DEBUG_OPT));
        config.add(new ConfigParam<string>
                   ("", "--sample-popsize-config", "<popsize config file>",
                    &popsize_config_file, "",
                    "optional, for use with --sample-popsize: should have a"
                    " line for each time interval, starting with the most"
                    " recent. Each line can have up to three tab-separated"
                    " columns, but only the first is required. The first column"
                    " gives the name of the popsize parameter- any time"
                    " intervals with the same entry here will be constrained"
                    " to have the same popsize. The second line is the initial"
                    " value for the parameter, and the third should be a 1 or"
                    " 0 indicating whether to sample that parameter. The second"
                    " and third columns are optional; by default all parameters"
                    " will be sampled when --sample-popsize is used. For rows"
                    " with the same value in the first column, the second and"
                    " third columns should also be the same."));
        config.add(new ConfigParam<int>
                   ("", "--sample-popsize-num", "<num>", &sample_popsize_num, 1,
                    "number of times to sample popsize per threading operation"
                    " (default=1)",
                    DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--popsize-prior-neighbor", &popsize_prior_neighbor,
                    "(for use with --sample-popsize) use prior that encourages"
                    " neighboring popsizes to be similar",
                    DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--sample-popsize-const", &sample_popsize_const,
                    "sample popsize, keep constant across times", DEBUG_OPT));
        config.add(new ConfigParam<int>
                   ("", "--sample-popsize-buildup", "<n>",
                    &sample_popsize_buildup, 0,
                    "for use with --sample-popsize, alternative to"
                    " --sample-popsize-config. Start off estimating single"
                    " popsize for all times, and every <n> iterations split"
                    " interval in half, until each time interval estimated"
                    " separately",
                    DEBUG_OPT));
	config.add(new ConfigParam<double>
		   ("", "--epsilon", "<val>", &epsilon,
		    0.01, "(for use with --sample-popsize) epsilon value for"
		    "Hamiltonian population size updates", DEBUG_OPT));
        config.add(new ConfigParam<double>
		   ("", "--pseudocount", "<val>", &pseudocount,
		    1.0, "(for use with --sample-popsize) gives weight to prior",
		    DEBUG_OPT));
#ifdef ARGWEAVER_MPI
        config.add(new ConfigParam<int>
                   ("", "--mcmcmc", "<int>", &mcmcmc_numgroup,
                    1, "number of mcmcmc threads",
                    DEBUG_OPT));
        config.add(new ConfigParam<double>
                   ("", "--mcmcmc-heat", "<val>", &mcmcmc_heat,
                    0.05, "heat interval for each thread in (MC)^3 group",
                    DEBUG_OPT));
#endif
        config.add(new ConfigSwitch
                   ("", "--init-popsize-random", &init_popsize_random,
                    "(for use with --sample-popsize). Initialize each"
                    " population size parameter to a random number sampled"
                    " uniformly in [5000, 50000]"));
        config.add(new ConfigParam<double>
                   ("-m", "--mutrate", "<mutation rate>", &mu, 2.5e-8,
                    "mutations per site per generation (default=2.5e-8)"));
        config.add(new ConfigParam<double>
                   ("-r", "--recombrate", "<recombination rate>", &rho, 1.5e-8,
                    "recombination per site per generation (default=1.5e-8)"));
        config.add(new ConfigParam<int>
                   ("-t", "--ntimes", "<ntimes>", &ntimes, 20,
                    "number of time points (default=20)"));
        config.add(new ConfigParam<double>
                   ("", "--maxtime", "<maxtime>", &maxtime, 200e3,
                    "maximum time point in generations (default=200e3)"));
        config.add(new ConfigParam<double>
                   ("", "--time-step", "<time>", &time_step, 0,
                    "linear time step in generations (optional)"));
        config.add(new ConfigParam<double>
                   ("", "--delta", "<delta>", &delta, 0.01,
                    "delta value for choosing log times (bigger value-> more"
                    " dense time points at leaves", DEBUG_OPT));
        config.add(new ConfigParam<string>
                   ("", "--times-file", "<times filename>", &times_file, "",
                    "file containing time points (optional)"));
        config.add(new ConfigParam<string>
                   ("-M", "--mutmap", "<mutation rate map file>", &mutmap, "",
                    "mutation map file (optional)"));
        config.add(new ConfigParam<string>
                   ("-R", "--recombmap", "<recombination rate map file>",
                    &recombmap, "",
                    "recombination map file (optional)"));

        // sampling
        config.add(new ConfigParamComment("Sampling"));
        config.add(new ConfigParam<int>
                   ("-n", "--iters", "<# of iterations>", &niters, 1000,
                    "(default=1000)"));
        config.add(new ConfigParam<string>
                   ("", "--resample-region", "<start>-<end>",
                    &resample_region_str, "",
                    "region to resample of input ARG (optional)"));
        config.add(new ConfigSwitch
                   ("", "--resume", &resume, "resume a previous run"));
        config.add(new ConfigSwitch
                   ("", "--overwrite", &overwrite,
                    "force an overwrite of a previous run"));
	config.add(new ConfigSwitch
		   ("", "--no-sample-arg", &no_sample_arg,
		    "Do not sample the ARG; hold at initial value "
		    "(usually for use with --sample-popsize)",
		    DEBUG_OPT));

        // misc
        config.add(new ConfigParamComment("Miscellaneous"));
        config.add(new ConfigParam<int>
                   ("-c", "--compress-seq", "<compression factor>",
                    &compress_seq, 1,
                    "alignment compression factor (default=1)"));
        config.add(new ConfigParam<int>
                   ("", "--climb", "<# of climb iterations>", &nclimb, 0,
                    "(default=0)"));
        config.add(new ConfigParam<int>
                   ("", "--sample-step", "<sample step size>", &sample_step,
                    10, "number of iterations between steps (default=10)"));
        config.add(new ConfigSwitch
                   ("", "--no-compress-output", &no_compress_output,
                    "do not use compressed output"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default=current time)"));

        // advance options
        config.add(new ConfigParamComment("Advanced Options", DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--gibbs", &gibbs,
                    "use Gibbs sampling"));
        config.add(new ConfigParam<double>
                   ("", "--prob-path-switch", "<probability>",
                    &prob_path_switch, .1,
                    "removal path switch (default=.1)", DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--infsites", &infsites,
                    "assume infinite sites model (at most one mutation per site)",
                    DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--unphased", &unphased,
                    "data is unphased (will integrate over phasings). "
                    "Note: Experimental!", DEBUG_OPT));
        config.add(new ConfigParam<string>
                   ("", "--unphased-file", "<filename>", &unphased_file, "",
                    "use this file to identify haplotype pairs (file should"
                    " have two sequence names per line)", DEBUG_OPT));
        config.add(new ConfigParam<double>
                   ("", "--randomize-phase", "<frac_random>", &randomize_phase,
                    0.0, "randomize phasings at start", DEBUG_OPT));
        config.add(new ConfigParam<int>
                   ("", "--sample-phase", "<niters>", &sample_phase, 0,
                    "output phasings every <niters> samples", DEBUG_OPT));
        config.add(new ConfigParam<int>
                   ("", "--resample-window", "<window size>",
                    &resample_window, 100000,
                    "sliding window for resampling (default=100000)",
                    DEBUG_OPT));
        config.add(new ConfigParam<int>
                   ("", "--resample-window-iters", "<iterations>",
                    &resample_window_iters, 10,
                    "number of iterations per sliding window for resampling"
                    " (default=10)", DEBUG_OPT));

        // help information
        config.add(new ConfigParamComment("Information"));
        config.add(new ConfigParam<int>
                   ("-V", "--verbose", "<verbosity level>",
                    &verbose, LOG_LOW,
                    "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
        config.add(new ConfigSwitch
                   ("-q", "--quiet", &quiet, "suppress logging to stderr"));
        config.add(new ConfigSwitch
                   ("-v", "--version", &version, "display version information"));
        config.add(new ConfigSwitch
                   ("-h", "--help", &help,
                    "display help information"));
        config.add(new ConfigSwitch
                   ("", "--help-advanced", &help_debug,
                    "display help information about advanced options"));
    }

    int parse_args(int argc, char **argv)
    {
        // parse arguments
        if (!config.parse(argc, (const char**) argv)) {
            if (argc < 2)
                config.printHelp();
            return EXIT_ERROR;
        }

        // display help
        if (help) {
            config.printHelp();
            return EXIT_ERROR;
        }

        // display debug help
        if (help_debug) {
            config.printHelp(stderr, DEBUG_OPT);
            return EXIT_ERROR;
        }

        // display version info
        if (version) {
            printf(VERSION_INFO);
            return EXIT_ERROR;
        }
#ifdef ARGWEAVER_MPI
        mcmcmc_group = 0;
        int groupsize = MPI::COMM_WORLD.Get_size() / mcmcmc_numgroup;
        mcmcmc_group = MPI::COMM_WORLD.Get_rank() / groupsize;
        if (mcmcmc_group != 0) {
            char tmp[1000];
            sprintf(tmp, ".%i", mcmcmc_group);
            mcmcmc_prefix = string(tmp);
        }
        printLog(LOG_LOW, "mcmcmc_prefix = %s\n", mcmcmc_prefix.c_str());
        printLog(LOG_LOW, "mcmcmc_group=%i\n", mcmcmc_group);
#endif

        return 0;
    }

    ConfigParser config;

    // input/output
    string fasta_file;
    string sites_file;
    string subsites_file;
    string out_prefix;
    string arg_file;
    string subregion_str;

    // model parameters
    double popsize;
    string popsize_str;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    double time_step;
    double delta;
    string times_file;
    string mutmap;
    string recombmap;
    string maskmap;
    ArgModel model;
    bool sample_popsize;
    bool sample_popsize_const;
    bool popsize_prior_neighbor;
    int sample_popsize_buildup;
    bool sample_popsize_recomb;
    bool init_popsize_random;
    string popsize_config_file;
    int sample_popsize_num;
    double epsilon;
    double pseudocount;
#ifdef ARGWEAVER_MPI
    double mcmcmc_heat;
    int mcmcmc_group;
    int mcmcmc_numgroup;
    bool mpi;
#endif
    string mcmcmc_prefix;
    bool no_sample_arg;

    // search
    int nclimb;
    int niters;
    string resample_region_str;
    int resample_region[2];
    bool resume;
    bool overwrite;
    string resume_stage;
    int resume_iter;
    int resample_window;
    int resample_window_iters;
    bool gibbs;

    // misc
    int compress_seq;
    int sample_step;
    bool no_compress_output;
    int randseed;
    double prob_path_switch;
    bool infsites;
    bool unphased;
    string unphased_file;
    double randomize_phase;
    int sample_phase;

    // help/information
    bool quiet;
    int verbose;
    bool version;
    bool help;
    bool help_debug;

    // logging
    FILE *stats_file;
};



bool parse_region(const char *region, int *start, int *end)
{
    return sscanf(region, "%d-%d", start, end) == 2;
}

//=============================================================================
// logging

void set_up_logging(const Config &c, int level, const char *log_mode) {
    Logger *logger;
    printf("set_up_logging %s %s %s %i %i\n",
           c.out_prefix.c_str(), c.mcmcmc_prefix.c_str(), LOG_SUFFIX, c.quiet,
           level);
    setLogLevel(level);
    if (c.quiet) {
        // log only to file
        logger = &g_logger;
    } else {
        // log to both stdout and file
        logger = new Logger(NULL, level);
        g_logger.setChain(logger);
    }
    string log_filename = c.out_prefix + c.mcmcmc_prefix + LOG_SUFFIX;
    if (!logger->openLogFile(log_filename.c_str(), log_mode)) {
        printError("Could not open log file '%s'", log_filename.c_str());
        abort();
    }
}

// log the program version and start time
void log_intro(int level)
{
    time_t t = time(NULL);

    printLog(level, "arg-sample " VERSION_TEXT "\n");
    printLog(level, "start time: %s", ctime(&t));  // newline is in ctime
}

// log the command used
void log_prog_commands(int level, int argc, char **argv)
{
    printLog(level, "command:");
    for (int i=0; i<argc; i++) {
        printLog(level, " %s", argv[i]);
    }
    printLog(level, "\n");
}


// log the model used
void log_model(const ArgModel &model)
{
    printLog(LOG_LOW, "\n");
    printLog(LOG_LOW, "model: \n");
    printLog(LOG_LOW, "  mu = %e\n", model.mu);
    printLog(LOG_LOW, "  rho = %e\n", model.rho);
    printLog(LOG_LOW, "  ntimes = %d\n", model.ntimes);
    printLog(LOG_LOW, "  times = [");
    for (int i=0; i<model.ntimes-1; i++)
        printLog(LOG_LOW, "%f,", model.times[i]);
    printLog(LOG_LOW, "%f]\n", model.times[model.ntimes-1]);
    printLog(LOG_LOW, "  popsizes = [");
    for (int i=0; i<2*model.ntimes-2; i++)
        printLog(LOG_LOW, "%f,", model.popsizes[i]);
    printLog(LOG_LOW, "%f]\n", model.popsizes[2*model.ntimes-2]);

    if (isLogLevel(LOG_HIGH)) {
        printLog(LOG_HIGH, "mutmap = [\n");
        for (unsigned int i=0; i<model.mutmap.size(); i++) {
            printLog(LOG_HIGH, "%d\t%d\t%e\n",
                     model.mutmap[i].start, model.mutmap[i].end,
                     model.mutmap[i].value);
        }
        printLog(LOG_HIGH, "]\n");

        printLog(LOG_HIGH, "recombmap = [\n");
        for (unsigned int i=0; i<model.recombmap.size(); i++) {
            printLog(LOG_HIGH, "%d\t%d\t%e\n",
                     model.recombmap[i].start, model.recombmap[i].end,
                     model.recombmap[i].value);
        }
        printLog(LOG_HIGH, "]\n");
    }

    printLog(LOG_LOW, "\n");
}


//=============================================================================
// alignment compression

template<class T>
void compress_track(Track<T> &track, SitesMapping *sites_mapping,
                    double compress_seq, bool is_rate)
{
    Track<T> track2;

    if (sites_mapping) {
        // get block lengths
        vector<int> blocks;
        for (unsigned int i=0; i<track.size(); i++)
            blocks.push_back(track[i].length());

        // compress block lengths
        vector<int> blocks2;
        sites_mapping->compress_blocks(blocks, blocks2);

        // build compress track
        int start = sites_mapping->new_start;
        for (unsigned int i=0; i<track.size(); i++) {
            int end = start + blocks2[i];
            if (end > start)
                track2.append(track[i].chrom, start, end, track[i].value);
            start = end;
        }

        // replace track
        track.clear();
        track.insert(track.begin(), track2.begin(), track2.end());
    }


    // compress rate
    if (is_rate) {
        for (unsigned int i=0; i<track.size(); i++)
            track[i].value *= compress_seq;
    }
}


template<class T>
void compress_mask(Track<T> &track, SitesMapping *sites_mapping)
{
    if (sites_mapping) {
        int prev_start_orig = 0, prev_start_new = 0;
        for (unsigned int i=0; i<track.size(); i++) {
            if (track[i].start < prev_start_orig)
                prev_start_new = 0;
            prev_start_orig = track[i].start;
            track[i].start = sites_mapping->compress(track[i].start,
                                                     prev_start_new, 1);
            prev_start_new = track[i].start;
            track[i].end = sites_mapping->compress(track[i].end,
                                                   track[i].start, -1);
        }
    }
}


void compress_model(ArgModel *model, SitesMapping *sites_mapping,
                    double compress_seq)
{
    model->rho *= compress_seq;
    model->mu *= compress_seq;

    compress_track(model->mutmap, sites_mapping, compress_seq, true);
    compress_track(model->recombmap, sites_mapping, compress_seq, true);
}




//=============================================================================
// statistics output

void print_stats_header(Config *config) {
    fprintf(config->stats_file, "stage\titer\tprior\tlikelihood\tjoint\t"
            "recombs\tnoncompats\targlen");
    if (config->model.popsize_config.config_buildup) {
        for (int i=0; i < 2*config->model.ntimes-1; i++)
            fprintf(config->stats_file, "\tN%i", i);
    } else if (config->model.popsize_config.sample) {
        list<PopsizeConfigParam> l = config->model.popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin();
             it != l.end(); ++it) {
            fprintf(config->stats_file, "\t%s", it->name.c_str());
        }
    }
    fprintf(config->stats_file, "\n");
}


void print_stats(FILE *stats_file, const char *stage, int iter,
                 ArgModel *model,
                 const Sequences *sequences, LocalTrees *trees,
                 const SitesMapping* sites_mapping, const Config *config)
{
    // calculate number of recombinations
    int nrecombs = trees->get_num_trees() - 1;

    // calculate number of non-compatiable sites
    int nseqs = trees->get_num_leaves();
    char *seqs[nseqs];
    for (int i=0; i<nseqs; i++)
        seqs[i] = sequences->seqs[trees->seqids[i]];
    int noncompats = count_noncompat(trees, seqs, nseqs, sequences->length());

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;


    // calculate likelihood, prior, and joint probabilities
    // uncompressed local trees
    if (sites_mapping)
        uncompress_local_trees(trees, sites_mapping);

    double prior = calc_arg_prior(model, trees);
    double likelihood = calc_arg_likelihood(model, sequences, trees,
                                            sites_mapping);
    double joint = prior + likelihood;
    double arglen = get_arglen(trees, model->times);

    // recompress local trees
    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);

    // output stats
    fprintf(stats_file, "%s\t%d\t%f\t%f\t%f\t%d\t%d\t%f",
            stage, iter,
            prior, likelihood, joint, nrecombs, noncompats, arglen);
    if (model->popsize_config.config_buildup) {
        for (int i=0; i < 2*model->ntimes-1; i++)
            fprintf(stats_file, "\t%f", model->popsizes[i]);
    } else if (model->popsize_config.sample) {
        list<PopsizeConfigParam> l=model->popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin();
             it != l.end(); ++it) {
            set<int>::iterator it2 = it->pops.begin();
            int val = *it2;
            fprintf(stats_file, "\t%f", model->popsizes[val]);
            it2++;
            //just checking here; can delete later
            while (it2 != it->pops.end()) {
                assert(model->popsizes[*it2] == model->popsizes[*it2]);
                it2++;
            }
        }
    }
    fprintf(stats_file, "\n");
    fflush(stats_file);

    printLog(LOG_LOW, "\n"
             "prior:      %f\n"
             "likelihood: %f\n"
             "joint:      %f\n"
             "nrecombs:   %d\n"
             "noncompats: %d\n"
             "arglen:     %f\n"
             "max memory: %.1f MB\n\n",
             prior, likelihood, joint, nrecombs, noncompats, arglen, maxrss);

}

//=============================================================================
// sample output


// Returns the iteration-specific ARG filename
string get_out_arg_file(const Config &config, int iter)
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    return config.out_prefix + config.mcmcmc_prefix + iterstr + SMC_SUFFIX;
}

string get_out_sites_file(const Config &config, int iter)
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    return config.out_prefix + config.mcmcmc_prefix + iterstr + SITES_SUFFIX;
}



bool log_sequences(string chrom, const Sequences *sequences,
                   const Config *config,
                   const SitesMapping *sites_mapping, int iter) {
    Sites sites(chrom);
    string out_sites_file = get_out_sites_file(*config, iter);
    make_sites_from_sequences(sequences, &sites);
    if (!config->no_compress_output)
        out_sites_file += ".gz";
    if (sites_mapping)
        uncompress_sites(&sites, sites_mapping);
    CompressStream stream(out_sites_file.c_str(), "w");
    if (!stream.stream) {
        printError("cannot write '%s'", out_sites_file.c_str());
        return false;
    }
    write_sites(stream.stream, &sites);
    return true;
}

bool log_local_trees(const ArgModel *model, const Sequences *sequences,
                     LocalTrees *trees, const SitesMapping* sites_mapping,
                     const Config *config, int iter)
{
    string out_arg_file = get_out_arg_file(*config, iter);
    if (!config->no_compress_output)
        out_arg_file += ".gz";

    // write local trees uncompressed
    if (sites_mapping)
        uncompress_local_trees(trees, sites_mapping);

    // setup output stream
    CompressStream stream(out_arg_file.c_str(), "w");
    if (!stream.stream) {
        printError("cannot write '%s'", out_arg_file.c_str());
        return false;
    }


    write_local_trees(stream.stream, trees, sequences, model->times);

    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);

    return true;
}


//=============================================================================


// Read a list of doubles from a file.
bool read_doubles(const char *filename, vector<double> &values) {
    FILE *infile = fopen(filename, "r");
    if (infile == NULL) {
        return false;
    }

    // Read all doubles and sort them.
    double value;
    while (fscanf(infile, "%lf", &value) == 1)
        values.push_back(value);
    fclose(infile);

    return true;
}


bool read_init_arg(const char *arg_file, const ArgModel *model,
                   LocalTrees *trees, vector<string> &seqnames)
{
    CompressStream stream(arg_file, "r");
    if (!stream.stream) {
        printError("cannot read '%s'", arg_file);
        return false;
    }

    return read_local_trees(stream.stream, model->times, model->ntimes,
                            trees, seqnames);
}



//=============================================================================
// sampling methods

// build initial arg by sequential sampling
void seq_sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                    SitesMapping* sites_mapping, Config *config)
{
    if (trees->get_num_leaves() < sequences->get_num_seqs()) {
        printLog(LOG_LOW, "Sequentially Sample Initial ARG (%d sequences)\n",
                 sequences->get_num_seqs());
        printLog(LOG_LOW, "------------------------------------------------\n");
        sample_arg_seq(model, sequences, trees, true);
        print_stats(config->stats_file, "seq", trees->get_num_leaves(),
                    model, sequences, trees, sites_mapping, config);
    }
}


void climb_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
               SitesMapping* sites_mapping, Config *config)
{
    if (config->resume)
        return;

    printLog(LOG_LOW, "Climb Search (%d iterations)\n", config->nclimb);
    printLog(LOG_LOW, "-----------------------------\n");
    double recomb_preference = .9;
    for (int i=0; i<config->nclimb; i++) {
        printLog(LOG_LOW, "climb %d\n", i+1);
        resample_arg_climb(model, sequences, trees, recomb_preference);
        print_stats(config->stats_file, "climb", i, model, sequences, trees,
                    sites_mapping, config);
    }
    printLog(LOG_LOW, "\n");
}

void mcmcmc_swap(Config *config, ArgModel *model, const Sequences *sequences,
                 const LocalTrees *trees, const SitesMapping *sites_mapping) {
#ifdef ARGWEAVER_MPI
    printLog(LOG_LOW, "mcmcmc_swap model->mc3.max_group=%i\n", model->mc3.max_group);
    if (model->mc3.max_group == 0) return;
    Mc3Config *mc3 = &(model->mc3);
    int swap[2];
    int rank = MPI::COMM_WORLD.Get_rank();
    bool accept;
    double swapstats[4];
    //pick two groups to swap
    if (rank == 0) {
        swap[0] = irand(mc3->max_group + 1);
        swap[1] = irand(mc3->max_group);
        if (swap[1] >= swap[0]) swap[1]++;
        assert(swap[0] >= 0 && swap[0] <= mc3->max_group);
        assert(swap[1] >= 0 && swap[1] <= mc3->max_group);
        assert(swap[0] != swap[1]);
    }
    MPI::COMM_WORLD.Bcast(swap, 2, MPI::INT, 0);
    if (mc3->group == swap[0] || mc3->group == swap[1]) {
        double vals[2];
        vals[0] = calc_arg_prior(model, trees) +
            calc_arg_likelihood(model, sequences, trees, sites_mapping);
        vals[1] = mc3->heat;
        if (mc3->group_comm->Get_rank()==0)
            mc3->group_comm->Reduce(MPI_IN_PLACE, vals, 1, MPI::DOUBLE, MPI_SUM,
                                   0);
        else mc3->group_comm->Reduce(vals, vals, 1, MPI::DOUBLE, MPI_SUM, 0);
        if (mc3->group_comm->Get_rank()==0)
            MPI::COMM_WORLD.Isend(vals, 2, MPI::DOUBLE, 0, 301);
    }
    if (rank == 0) {
        double like[2];
        double heat[2];
        for (int i=0; i < 2; i++) {
            double vals[2];
            MPI::COMM_WORLD.Recv(vals, 2, MPI::DOUBLE, MPI_ANY_SOURCE, 301);
            like[i] = vals[0];
            heat[i] = vals[1];
        }
        double accept_ratio = (heat[0]-heat[1])*like[1] +
            (heat[1] - heat[0])*like[0];
        accept = (accept_ratio >= 0.0 ||
                  frand() < exp(accept_ratio));
        //Note heat[0] may correspond to swap[0] or swap[1]
        swapstats[0] = accept ? 1.0 : 0.0;
        swapstats[1] = heat[0];
        swapstats[2] = heat[1];
        swapstats[3] = accept_ratio;
    }
    MPI::COMM_WORLD.Bcast(swapstats, 4, MPI::DOUBLE, 0);
    accept = ( swapstats[0] > 0.5);
    printLog(LOG_LOW, "swap\t%i\t%i\t%f\t%f\t%f\t%s\n",
             swap[0], swap[1], swapstats[1], swapstats[2], swapstats[3],
             accept ? "accept" : "reject");

    // now all processes know accept. Do the switch
    if (accept && (mc3->group == swap[0] || mc3->group == swap[1])) {
        int other = (mc3->group == swap[1] ? 0 : 1);
        mc3->group = swap[other];
        mc3->heat = 1.0 - mc3->heat_interval * mc3->group;
        //need to switch output files as well, including stats_file, arg output,
        //phase output, log files.  First close all the files.
        fclose(config->stats_file);
        if (config->verbose) {
            Logger *chain = g_logger.getChain();
            chain->closeLogFile();
        }
        if (mc3->group == 0) config->mcmcmc_prefix = "";
        else {
            char tmp[1000];
            sprintf(tmp, ".%i", mc3->group);
            config->mcmcmc_prefix = string(tmp);
        }
    }
    MPI::COMM_WORLD.Barrier();
    //Now reopen in append mode
    if (accept && (mc3->group == swap[0] || mc3->group == swap[1])) {
        string stats_filename = config->out_prefix + config->mcmcmc_prefix
            + STATS_SUFFIX;
        if (!(config->stats_file = fopen(stats_filename.c_str(), "a"))) {
            printError("Error reopening stats file %s in mcmcmc_swap\n",
                    stats_filename.c_str());
            abort();
        }
        if (config->verbose) {
            string log_filename = config->out_prefix + config->mcmcmc_prefix
                + LOG_SUFFIX;
            Logger *chain = g_logger.getChain();
            if (!chain->openLogFile(log_filename.c_str(), "a")) {
                fprintf(stderr, "Error opening %s\n", log_filename.c_str());
                abort();
            }
        }
    }
#endif
}


void resample_arg_all(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                      SitesMapping* sites_mapping, Config *config)
{
    // setup search options
    bool do_leaf[config->niters+1];
    //int window = 100000;
    //int niters = 10;
    int window = config->resample_window;
    int niters = config->resample_window_iters;
    window /= config->compress_seq;
    int step = window / 2;

    // set iteration counter
    int iter = 1;
    if (config->resume)
        iter = config->resume_iter + 1;
    else {
        // save first ARG (iter=0)
        print_stats(config->stats_file, "resample", 0, model, sequences, trees,
                    sites_mapping, config);
        log_local_trees(model, sequences, trees, sites_mapping, config, 0);
    }


    printLog(LOG_LOW, "Resample All Branches (%d iterations)\n",
             config->niters);
    printLog(LOG_LOW, "--------------------------------------\n");

    double frac_leaf = 0.5;
    for (int i=0; i <= config->niters; i++) do_leaf[i] = (frand() < frac_leaf);

#ifdef ARGWEAVER_MPI
    MPI::COMM_WORLD.Bcast(do_leaf, config->niters+1, MPI::BOOL, 0);
#endif

    for (int i=iter; i<=config->niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i);
        Timer timer;
        double heat = model->mc3.heat;

	if ( ! config->no_sample_arg) {
	    if (config->gibbs)
		resample_arg(model, sequences, trees);
	    else
		resample_arg_mcmc_all(model, sequences, trees, do_leaf[i],
				      window, step, niters, heat);
	}

        if (model->popsize_config.sample) {
            if (model->popsize_config.config_buildup > 0 &&
                i > 0 && i % model->popsize_config.config_buildup == 0)
                model->popsize_config.split_config();
	    //            resample_popsizes(model, trees, config->sample_popsize_recomb, heat);
	    //	    mle_popsize(model, trees);
	    update_popsize_hmc(model, trees);
        }

        printTimerLog(timer, LOG_LOW, "sample time:");

        mcmcmc_swap(config, model, sequences, trees, sites_mapping);

        // logging
        print_stats(config->stats_file, "resample", i, model, sequences, trees,
                    sites_mapping, config);

        // sample saving
        if (i % config->sample_step == 0 && ! config->no_sample_arg)
            log_local_trees(model, sequences, trees, sites_mapping, config, i);

        if (config->sample_phase > 0 && i%config->sample_phase == 0)
            log_sequences(trees->chrom, sequences, config, sites_mapping, i);
    }
    printLog(LOG_LOW, "\n");
}


// overall sampling workflow
void sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                SitesMapping* sites_mapping, Config *config)
{
    if (!config->resume)
        print_stats_header(config);

    // build initial arg by sequential sampling
    seq_sample_arg(model, sequences, trees, sites_mapping, config);

    if (config->resample_region[0] != -1) {
        // region sampling
        printLog(LOG_LOW, "Resample Region (%d-%d, %d iterations)\n",
                 config->resample_region[0], config->resample_region[1],
                 config->niters);
        printLog(LOG_LOW, "--------------------------------------------\n");

        print_stats(config->stats_file, "resample_region", 0,
                    model, sequences, trees, sites_mapping, config);

        resample_arg_region(model, sequences, trees,
                            config->resample_region[0],
                            config->resample_region[1],
                            config->niters);

        // logging
        print_stats(config->stats_file, "resample_region", config->niters,
                    model, sequences, trees, sites_mapping, config);
        log_local_trees(model, sequences, trees, sites_mapping, config, 0);

    } else{
        // climb sampling
        climb_arg(model, sequences, trees, sites_mapping, config);
        // resample all branches
        resample_arg_all(model, sequences, trees, sites_mapping, config);
    }
}


//=============================================================================

bool parse_status_line(const char* line, const Config &config,
                       string &stage, int &iter, string &arg_file,
                       vector<string> header)
{
    // parse stage and last iter
    vector<string> tokens;
    split(line, "\t", tokens);
    if (tokens.size() < 2) {
        printError("incomplete line in status file");
        return false;
    }

    string stage2 = tokens[0];
    int iter2;
    if (sscanf(tokens[1].c_str(), "%d", &iter2) != 1) {
        printError("iter column is not an integer");
        return false;
    }

    // NOTE: only resume resample stage for now
    if (stage2 != "resample")
        return true;

    // see if ARG file exists
    string out_arg_file = get_out_arg_file(config, iter2);
    struct stat st;
    if (stat(out_arg_file.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_arg_file;
    }

    // try compress output
    out_arg_file += ".gz";
    if (stat(out_arg_file.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_arg_file;
    }

    //TODO: need to make this work with popsize_config.config_buildup
    if (config.model.popsize_config.sample) {
        list<PopsizeConfigParam> l=config.model.popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin(); it != l.end();
             ++it) {
            string popname=it->name;
            set<int> popset = it->pops;
            int found=false;
            for (unsigned int i=0; i < header.size(); i++) {
                if (header[i] == popname) {
                    found=true;
                    double tempN;
                    sscanf(tokens[i].c_str(), "%lf", &tempN);
                    for (set<int>::iterator it2=popset.begin();
                         it2 != popset.end(); ++it2) {
                        int pop = *it2;
                        if (pop < 0 || pop >= 2*config.model.ntimes-1) {
                            printError("Error in resume: popsize config does"
                                       " not match previous run\n");
                            abort();
                        }
                        config.model.popsizes[pop] = tempN;
                    }
                    break;
                }
            }
            if (!found) {
                printError("Error in resume: did not find pop %s in previous"
                           " run stats file\n",
                           popname.c_str());
                abort();
            }
        }
    }
    return true;
}


bool setup_resume(Config &config)
{
    if (!config.resume)
        return true;

    printLog(LOG_LOW, "Resuming previous run\n");

    // open stats file
    string stats_filename = config.out_prefix + config.mcmcmc_prefix
        + STATS_SUFFIX;
    printLog(LOG_LOW, "Checking previous run from stats file: %s\n",
             stats_filename.c_str());

    FILE *stats_file;
    if (!(stats_file = fopen(stats_filename.c_str(), "r"))) {
        printError("could not open stats file '%s'",
                   stats_filename.c_str());
        return false;
    }

    // find last line of stats file that has a written ARG
    char *line = NULL;

    line = fgetline(stats_file);
    if (!line) {
        printError("status file is empty");
        return false;
    }
    chomp(line);
    vector<string> header;
    split(line, "\t", header);
    delete [] line;

    // loop through status lines
    string arg_file = "";
    while ((line = fgetline(stats_file))) {
        if (!parse_status_line(
             line, config, config.resume_stage, config.resume_iter, arg_file,
             header)) {
            delete [] line;
            return false;
        }
        delete [] line;
    }

    if (arg_file == "") {
        printLog(LOG_LOW, "Could not find any previously written ARG files. "
                 "Try disabling resume\n");
        return false;
    }
    config.arg_file = arg_file;

    printLog(LOG_LOW, "resuming at stage=%s, iter=%d, arg=%s\n",
             config.resume_stage.c_str(), config.resume_iter,
             config.arg_file.c_str());

    // clean up
    fclose(stats_file);

    return true;
}


bool check_overwrite(Config &config)
{
    // check for stats file
    string stats_filename = config.out_prefix + config.mcmcmc_prefix
        + STATS_SUFFIX;
    bool exists = !access(stats_filename.c_str(), F_OK);
    if (config.resume || config.overwrite || !exists) {
        return true;
    } else {
        printError("Stats file already exists: %s", stats_filename.c_str());
        printError("To force overwrite use --overwrite option");
        return false;
    }
}


bool ensure_output_dir(const char *outdir)
{
    char *path = strdup(outdir);
    char *dir = dirname(path);
    bool result = true;
    if (!makedirs(dir)) {
        printError("could not make directory for output files '%s'", dir);
        result = false;
    }
    free(path);
    return result;
}


//=============================================================================


int main(int argc, char **argv)
{

#ifdef ARGWEAVER_MPI
    MPI::Init(argc, argv);
#endif

    // parse command line arguments
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
        return ret;

    // ensure output dir
    if (!ensure_output_dir(c.out_prefix.c_str()))
        return EXIT_ERROR;

    // check overwriting
    if (!check_overwrite(c))
        return EXIT_ERROR;

#ifdef ARGWEAVER_MPI
    if (c.mpi) {
        int numcore = MPI::COMM_WORLD.Get_size();
        if (numcore % c.mcmcmc_numgroup != 0) {
            fprintf(stderr, "Error: number of cores should be evenly divisible"
                    " by number of mcmcmc threads");
        }
        int groupsize = numcore / c.mcmcmc_numgroup;
        int sites_num = MPI::COMM_WORLD.Get_rank() % groupsize;
        char tmp[10000];
        sprintf(tmp, "%s%i.sites", c.sites_file.c_str(),
                sites_num);
        c.sites_file = (string)tmp;
        sprintf(tmp, "%s%i", c.out_prefix.c_str(),
                sites_num);
        c.out_prefix = (string)tmp;
	if (c.arg_file != "") {
	    sprintf(tmp, "%s%i.smc.gz", c.arg_file.c_str(), sites_num);
	    c.arg_file = (string)tmp;
	}
    }
#endif

    // setup logging
    set_up_logging(c, c.verbose, (c.resume ? "a" : "w"));

    // log intro
    if (c.resume)
        printLog(LOG_LOW, "RESUME\n");
    log_intro(LOG_LOW);
    log_prog_commands(LOG_LOW, argc, argv);
    Timer timer;


    // init random number generator
    if (c.randseed == 0)
        c.randseed = time(NULL);
#ifdef ARGWEAVER_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
        for (int i=1; i < MPI::COMM_WORLD.Get_size(); i++) {
            int seed = irand(12581020);
            MPI::COMM_WORLD.Send(&seed, 1, MPI::INT, i, 13);
        }
    } else {
        MPI::COMM_WORLD.Recv(&c.randseed, 1, MPI::INT, 0, 13);
    }
#endif
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n", c.randseed);

    // read sequences
    Sites sites;
    Sequences sequences;
    SitesMapping *sites_mapping = NULL;
    auto_ptr<SitesMapping> sites_mapping_ptr;
    Region seq_region;
    Region seq_region_compress;


    if (c.fasta_file != "") {
        // read FASTA file

        if (!read_fasta(c.fasta_file.c_str(), &sequences)) {
            printError("could not read fasta file");
            return EXIT_ERROR;
        }
        seq_region.set("chr", 0, sequences.length());

        printLog(LOG_LOW, "read input sequences (nseqs=%d, length=%d)\n",
                 sequences.get_num_seqs(), sequences.length());

        // if compress or subset requested, make sites object
        if (c.compress_seq > 1 || c.subsites_file != "")
            make_sites_from_sequences(&sequences, &sites);

    } else if (c.sites_file != "") {
        // read sites file

        // parse subregion if given
        int subregion[2] = {-1, -1};
        if (c.subregion_str != "") {
            if (!parse_region(c.subregion_str.c_str(),
                              &subregion[0], &subregion[1])) {
                printError("subregion is not specified as 'start-end'");
                return EXIT_ERROR;
            }
            subregion[0] -= 1; // convert to 0-index
        }

        // read sites
        CompressStream stream(c.sites_file.c_str());
        if (!stream.stream ||
            !read_sites(stream.stream, &sites, subregion[0], subregion[1])) {
            printError("could not read sites file");
            return EXIT_ERROR;
        }
        stream.close();

        printLog(LOG_LOW, "read input sites (chrom=%s, start=%d, end=%d, "
                 "length=%d, nseqs=%d, nsites=%d)\n",
                 sites.chrom.c_str(), sites.start_coord, sites.end_coord,
                 sites.length(), sites.get_num_seqs(),
                 sites.get_num_sites());

        // sanity check for sites
        if (sites.get_num_sites() == 0) {
            printLog(LOG_LOW, "no sites given.  terminating.\n");
            return EXIT_ERROR;
        }
        seq_region.set(sites.chrom, sites.start_coord, sites.end_coord);

    } else {
        // no input sequence specified
        printError("must specify sequences (use --fasta or --sites)");
        return EXIT_ERROR;
    }

    if (c.subsites_file != "") {
        set<string> keep;
        FILE *infile;
        if (!(infile = fopen(c.subsites_file.c_str(), "r"))) {
            printError("Could not open subsites file %s",
                       c.subsites_file.c_str());
            return false;
        }
        char name[10000];
        while (EOF != fscanf(infile, "%s", name))
            keep.insert(string(name));
        fclose(infile);
        if (sites.subset(keep)) {
            printError("Error subsetting sites\n");
            return false;
        }
    }

    //read in mask
    TrackNullValue maskmap;
    if (c.maskmap != "") {
        //read mask
        CompressStream stream(c.maskmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &maskmap, seq_region)) {
            printError("cannot read mask map '%s'",
                       c.maskmap.c_str());
            return EXIT_ERROR;
        }
    }

    // compress sequences
    if (sites.get_num_sites() > 0) {
        // first remove any sites that fall under mask

        sites_mapping = new SitesMapping();
        sites_mapping_ptr = auto_ptr<SitesMapping>(sites_mapping);

        if (!find_compress_cols(&sites, c.compress_seq, sites_mapping)) {
            printError("unable to compress sequences at given compression level"
                       " (--compress-seq)");
            return EXIT_ERROR;
        }
        compress_sites(&sites, sites_mapping);
        make_sequences_from_sites(&sites, &sequences);
    }
    seq_region_compress.set(seq_region.chrom, 0, sequences.length());

    // compress mask
    if (c.maskmap != "") {
        // apply mask
        // TODO: when mask is compressed, only allow it to apply to
        // whole compressed regions
        if (sites_mapping)
            compress_mask(maskmap, sites_mapping);
        apply_mask_sequences(&sequences, maskmap);

        // report number of masked sites
        bool *masked = new bool [sequences.length()];
        find_masked_sites(sequences.get_seqs(), sequences.get_num_seqs(),
                          sequences.length(), masked);
        int nmasked = 0;
        for (int i=0; i<sequences.length(); i++)
            nmasked += int(masked[i]);
        delete [] masked;
        printLog(LOG_LOW, "masked %d (%.1f%%) sites\n", nmasked,
                 100.0 * nmasked / double(sequences.length()));
    }


    // setup model parameters
    if (c.times_file != "")
        c.model.set_times_from_file(c.times_file);
    else if (c.time_step)
        c.model.set_linear_times(c.time_step, c.ntimes);
    else
        c.model.set_log_times(c.maxtime, c.ntimes, c.delta);
    c.model.rho = c.rho;
    c.model.mu = c.mu;
    const double infsites_penalty = 1e-100; // TODO: make configurable
    if (c.infsites)
        c.model.infsites_penalty = infsites_penalty;
    c.model.set_popsizes(c.popsize_str, c.model.ntimes);
    if (c.unphased_file != "")
        c.model.unphased_file = c.unphased_file;
    sequences.set_pairs(&c.model);
    if (c.randomize_phase > 0.0) {
        sequences.randomize_phase(c.randomize_phase);
    }
    if (c.unphased)
        c.model.unphased = true;
    c.model.sample_phase = c.sample_phase;
    if (c.sample_popsize_const)
        c.sample_popsize=true;
    if (c.sample_popsize) {
        if (c.sample_popsize_buildup) {
            if (c.popsize_config_file != "" || c.sample_popsize_const) {
                printError("Error: cannot use --sample-popsize-buildup with"
                           " --sample-popsize-const or"
                           " --sample-popsize-config\n");
                return 1;
            }
            c.model.popsize_config = PopsizeConfig(c.ntimes, true, true);
            c.model.popsize_config.config_buildup = c.sample_popsize_buildup;
        } else if (c.sample_popsize_const) {
            c.model.popsize_config =
                PopsizeConfig(c.ntimes, true, true);
        } else {
            c.model.popsize_config =
                PopsizeConfig(c.popsize_config_file, c.model.ntimes,
                              c.model.popsizes);
        }
        c.model.popsize_config.numsample = c.sample_popsize_num;
        c.model.popsize_config.neighbor_prior = c.popsize_prior_neighbor;
	c.model.popsize_config.epsilon = c.epsilon;
	c.model.popsize_config.pseudocount = c.pseudocount;
    }
#ifdef ARGWEAVER_MPI
    c.model.mc3 = Mc3Config(c.mcmcmc_group, c.mcmcmc_heat);

    printf("MPI rank=%i size=%i\n",
           MPI::COMM_WORLD.Get_rank(),
           MPI::COMM_WORLD.Get_size());
    MPI::COMM_WORLD.Barrier();
#endif
    if (c.init_popsize_random)
        c.model.set_popsizes_random();

    // read model parameter maps if given
    if (c.mutmap != "") {
        CompressStream stream(c.mutmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model.mutmap, seq_region)) {
            printError("cannot read mutation rate map '%s'", c.mutmap.c_str());
            return EXIT_ERROR;
        }
    }
    if (c.recombmap != "") {
        CompressStream stream(c.recombmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model.recombmap, seq_region)) {
            printError("cannot read recombination rate map '%s'",
                       c.recombmap.c_str());
            return EXIT_ERROR;
        }
    }

    // try to resume a previous run
    if (!setup_resume(c)) {
        printError("resume failed.");
        if (c.overwrite) {
            c.resume = false;
            printLog(LOG_LOW, "Resume failed.  Sampling will start from scratch"
                     " since overwrite is enabled.\n");
        } else {
            return EXIT_ERROR;
        }
    }

    // make compressed model
    ArgModel model(c.model);
    if (!model.setup_maps(seq_region.chrom, seq_region.start, seq_region.end))
        return EXIT_ERROR;
    compress_model(&model, sites_mapping, c.compress_seq);

    // log original model
    log_model(model);

    // setup init ARG
    LocalTrees *trees = NULL;
    auto_ptr<LocalTrees> trees_ptr;
    if (c.arg_file != "") {
        // init ARG from file

        trees = new LocalTrees();
        trees_ptr = auto_ptr<LocalTrees>(trees);
        vector<string> seqnames;
        if (!read_init_arg(c.arg_file.c_str(), &c.model, trees, seqnames)) {
            printError("could not read ARG");
            return EXIT_ERROR;
        }

        if (!trees->set_seqids(seqnames, sequences.names)) {
            printError("input ARG's sequence names do not match input"
                       " sequences");
            return EXIT_ERROR;
        }

        printLog(LOG_LOW, "read input ARG (chrom=%s, start=%d, end=%d,"
                 " nseqs=%d)\n",
                 trees->chrom.c_str(), trees->start_coord, trees->end_coord,
                 trees->get_num_leaves());

        // check ARG matches sites/sequences
        if (trees->start_coord != seq_region.start ||
            trees->end_coord != seq_region.end) {
            printError("trees range does not match sites: tree(start=%d,"
                       " end=%d), sites(start=%d, end=%d) [compressed"
                       " coordinates]",
                       trees->start_coord, trees->end_coord,
                       seq_region.start, seq_region.end);
            return EXIT_ERROR;
        }

        // compress input ARG if compression is requested
        if (sites_mapping)
            compress_local_trees(trees, sites_mapping, true);

    } else {
        // create new init ARG
        trees = new LocalTrees(seq_region_compress.start,
                               seq_region_compress.end);
        trees->chrom = seq_region.chrom;
        trees_ptr = auto_ptr<LocalTrees>(trees);
    }


    // check for region sample
    if (c.resample_region_str != "") {
        if (!parse_region(c.resample_region_str.c_str(),
                          &c.resample_region[0], &c.resample_region[1]))
            {
                printError("--resample-region is not specified as 'start-end'");
                return EXIT_ERROR;
            }
        c.resample_region[0] -= 1; // convert to 0-index

        if (sites_mapping) {
            c.resample_region[0] = sites_mapping->compress(c.resample_region[0], -1);
            c.resample_region[1] = sites_mapping->compress(c.resample_region[1], 1);
        }
    }


    // init stats file
    string stats_filename = c.out_prefix + c.mcmcmc_prefix + STATS_SUFFIX;
    const char *stats_mode = (c.resume ? "a" : "w");
    if (!(c.stats_file = fopen(stats_filename.c_str(), stats_mode))) {
        printError("could not open stats file '%s'", stats_filename.c_str());
        return EXIT_ERROR;
    }

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);

    // sample ARG
    printLog(LOG_LOW, "\n");
    sample_arg(&model, &sequences, trees, sites_mapping, &c);

    // final log message
    maxrss = get_max_memory_usage() / 1000.0;
    printTimerLog(timer, LOG_LOW, "sampling time: ");
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);
    printLog(LOG_LOW, "FINISH\n");

    // clean up
    fclose(c.stats_file);

#ifdef ARGWEAVER_MPI
    MPI_Finalize();
#endif

    return 0;
}
