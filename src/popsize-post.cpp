// C/C++ includes
#include <time.h>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>

// arghmm includes
#include "argweaver/compress.h"
#include "argweaver/ConfigParam.h"
#include "argweaver/emit.h"
#include "argweaver/fs.h"
#include "argweaver/logging.h"
#include "argweaver/mem.h"
#include "argweaver/parsing.h"
#include "argweaver/sample_arg.h"
#include "argweaver/sequences.h"
#include "argweaver/total_prob.h"
#include "argweaver/track.h"
#include "argweaver/est_popsize.h"
#include "argweaver/mcmcmc.h"


using namespace argweaver;

// version info
#define VERSION_TEXT "0.1"
#define VERSION_INFO  "\
Popsize-post " VERSION_TEXT " \n\
Melissa Hubisz\n\
Sample popsize from arg-sample ARGs\n\
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
	popsize = 10000.0;
	ntimes = 20;
	maxtime = 200e3;
	delta = 0.01;

        // note that population size is independent of mutation and
        // and recombination rate conditional on ARG so these are really just
        // dummy variables needed to define the model
	mu = 2.5e-8;
	rho = 1.25e-8;
    }

    void make_parser()
    {
        config.clear();

        // input/output
        config.add(new ConfigParam<string>
                   ("-o", "--output", "<outfile>", &outfile,
                    "popsizes.txt",
		    "file where popsize samples are written"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<arg-sample output directory>", &arg_dir, "",
		    "directory where arg-sample *.smc.gz files are"));
	config.add(new ConfigParam<int>
		   ("", "--arg-start", "<first sample number>", &arg_start, 500,
		    "starting iteration of ARGs to use"));
	config.add(new ConfigParam<int>
		   ("", "--arg-step", "<arg step interval>", &arg_step, 20,
		    "sampling step for input ARGs"));
	config.add(new ConfigParam<int>
		   ("", "--mpi", "<number of MPI arg-sample runs to combine>",
		    &mpi, 0, "Will add \".<n>\" prefixes to n arg-sample outputs,"
		    " for n from  0 to <mpi-1>"));

        // model parameters
        config.add(new ConfigParamComment("Model parameters"));
        config.add(new ConfigParam<string>
                   ("-N", "--popsize", "<population size>", &popsize_str,
                    "10000",
                    "initial effective population size (default=1e4)"));
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
	config.add(new ConfigSwitch
                   ("", "--sample-popsize-const", &sample_popsize_const,
		    "sample popsize, keep constant across times", DEBUG_OPT));
	config.add(new ConfigParam<double>
                   ("", "--time-step", "<time>", &time_step, 0,
                    "linear time step in generations (optional)"));
	config.add(new ConfigParam<double>
		   ("", "--epsilon", "<val>", &epsilon,
		    0.01, "(for use with --sample-popsize) epsilon value for"
		    "Hamiltonian population size updates", DEBUG_OPT));
        config.add(new ConfigParam<double>
		   ("", "--pseudocount", "<val>", &pseudocount,
		    0.0, "gives weight to prior",
		    DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--init-popsize-random", &init_popsize_random,
                    "(for use with --sample-popsize). Initialize each"
                    " population size parameter to a random number sampled"
                    " uniformly in [5000, 50000]"));
	config.add(new ConfigParam<double>
		   ("", "--min-events", "<val>", &min_events,
		    1.0, "Minimum number of recomb events affecting a time"
		    " interval; time intervals with fewer events will be"
		    " combined with adjacent, more ancient time interval (default: 1)"));

        // sampling
        config.add(new ConfigParamComment("Sampling"));
        config.add(new ConfigParam<int>
                   ("-n", "--iters", "<# of iterations>", &niters, 1000,
                    "(default=1000)"));
        config.add(new ConfigSwitch
                   ("", "--resume", &resume, "resume a previous run"));
        config.add(new ConfigSwitch
                   ("", "--overwrite", &overwrite,
                    "force an overwrite of a previous run"));

        // misc
        config.add(new ConfigParamComment("Miscellaneous"));
        config.add(new ConfigParam<int>
                   ("", "--sample-step", "<sample step size>", &sample_step,
                    10, "number of iterations between steps (default=10)"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default=current time)"));

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
        return 0;
    }

    ConfigParser config;

    // input/output
    string outfile;
    string arg_dir;
    int mpi;

    // model parameters
    double popsize;
    string popsize_str;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    double delta;
    double time_step;
    string times_file;
    ArgModel model;
    bool init_popsize_random;
    string popsize_config_file;
    int sample_popsize_num;
    bool sample_popsize_const;
    double epsilon;
    double pseudocount;
    int arg_start, arg_step;
    double min_events;

    // search
    int niters;
    bool resume;
    bool overwrite;
    string resume_stage;
    int resume_iter;
    int resample_window;
    int resample_window_iters;

    // misc
    int sample_step;
    int randseed;

    // help/information
    bool quiet;
    int verbose;
    bool version;
    bool help;
    bool help_debug;

    //logging
    FILE *stats_file;
};



//=============================================================================
// logging

void set_up_logging(const Config &c, int level, const char *log_mode) {
    //    Logger *logger;
    setLogLevel(level);
    // log only to stdout
    //    logger = &g_logger;
}

// log the program version and start time
void log_intro(int level)
{
    time_t t = time(NULL);

    printLog(level, "popsize-post " VERSION_TEXT "\n");
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

    }

    printLog(LOG_LOW, "\n");
}



//=============================================================================
// statistics output
void print_stats_header(Config *config) {
    //    fprintf(config->stats_file, "stage\titer\tprior\tlikelihood\tjoint");
    fprintf(config->stats_file, "stage\titer");
    for (int i=0; i < config->model.ntimes-1; i++)
	fprintf(config->stats_file, "\tN%i", i);
    fprintf(config->stats_file, "\n");
}

void print_stats_popsizes(Config *config, int iter, ArgModel *model) {
    fprintf(config->stats_file, "popsize_mle\t%i", iter);
    for (int i=0; i < config->model.ntimes-1; i++)
	fprintf(config->stats_file, "\t%.1f", model->popsizes[2*i]);
    fprintf(config->stats_file, "\n");
}


/*void print_stats(FILE *stats_file, const char *stage, int iter,
                 ArgModel *model, const struct popsize_data &popsize_data_prior,
		 const struct popsize_data *popsize_data_post, int popsize_data_len,
                 const Config *config)
{

    double prior = calc_popsize_prob(model, popsize_data_prior);
    double likelihood = calc_avg_popsize_prob(model, popsize_data_post, popsize_data_len);
    double joint = prior + likelihood;

    // output stats
    fprintf(stats_file, "%s\t%d\t%f\t%f\t%f", stage, iter, likelihood);
    for (int i=0; i < config->model.ntimes-1; i++)
	fprintf(stats_file, "\t%f", model->popsizes[i]);
    fprintf(stats_file, "\n");
    fflush(stats_file);

    printLog(LOG_LOW, "\n"
             "prior:      %f\n"
             "likelihood: %f\n"
             "joint:      %f\n",
             prior, likelihood, joint);

	     }*/

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


bool read_arg(const char *arg_file, const ArgModel *model,
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


int main(int argc, char **argv)
{

    // parse command line arguments
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
        return ret;

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
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n", c.randseed);

    // setup model parameters
    if (c.times_file != "")
        c.model.set_times_from_file(c.times_file);
    else if (c.time_step)
        c.model.set_linear_times(c.time_step, c.ntimes);
    else
        c.model.set_log_times(c.maxtime, c.ntimes, c.delta);
    c.model.rho = c.rho;
    c.model.mu = c.mu;
    c.model.set_popsizes(c.popsize_str, c.model.ntimes);
    c.model.popsize_config.pseudocount = c.pseudocount;

    // log original model
    ArgModel model(c.model);
    log_model(model);

    // init stats file
    string stats_filename = c.outfile;
    const char *stats_mode = (c.resume ? "a" : "w");
    if (!(c.stats_file = fopen(stats_filename.c_str(), stats_mode))) {
        printError("could not open stats file '%s'", stats_filename.c_str());
        return EXIT_ERROR;
    }

    print_stats_header(&c);
    //    ExtendArray<struct popsize_data> suff_stats;
    for (int rep=c.arg_start; ; rep += c.arg_step) {
	struct popsize_data data;
	int num_read=0;
	for (int mpi=0; mpi==0 || mpi < c.mpi; mpi++) {
	    char file[10000];
	    if (!c.mpi)
		sprintf(file, "%s.%i.smc.gz", c.arg_dir.c_str(), rep);
	    else sprintf(file, "%s%i.%i.smc.gz", c.arg_dir.c_str(), mpi, rep);
	    LocalTrees *trees = new LocalTrees();
	    bool new_arg;
	    vector<string> seqnames;
	    new_arg = read_arg(file, &c.model, trees, seqnames);
	    if (new_arg) {
		printLog(LOG_LOW, "read input ARG from %s\n", file);
		// TODO: set pseudocount to zero first? Keep that separate?
		//	    popsize_sufficient_stats(&data, model, trees);
		//	    suff_stats.push(&data);
		popsize_sufficient_stats(&data, &model, trees, mpi != 0);
	    }
	    delete trees;
	    if (!new_arg) break;
	    num_read++;
	}
	if (num_read == 0) break;
	if (c.mpi > 0 && num_read != c.mpi) break;
	mle_popsize(&model, &data, c.min_events);
	print_stats_popsizes(&c, rep, &model);
	delete_popsize_data(&data);
    }

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);

    // final log message
    maxrss = get_max_memory_usage() / 1000.0;
    printTimerLog(timer, LOG_LOW, "sampling time: ");
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);
    printLog(LOG_LOW, "FINISH\n");

    // clean up
    fclose(c.stats_file);

    return 0;
}
