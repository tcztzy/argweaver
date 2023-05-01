// C/C++ includes
#ifdef ARGWEAVER_MPI
#include "mpi.h"
#endif
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
#include "argweaver/coal_records.h"
#include "argweaver/recomb.h"
#include "argweaver/model.h"
#include "argweaver/local_tree.h"

using namespace argweaver;

// version info
#define VERSION_TEXT "0.8.1"
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
const char *COAL_RECORDS_SUFFIX = ".cr";


const int EXIT_ERROR = 1;


class MigEvent {
public:
    MigEvent(char *name0, int from_pop, int to_pop, int time) :
        from_pop(from_pop), to_pop(to_pop), time(time) {
        name = string(name0);
    }
    string name;
    int from_pop;
    int to_pop;
    int time;
};


// parsing command-line options
class Config
{
public:

    Config()
    {
        make_parser();
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
        config.add(new ConfigParam<int>
                   ("-c", "--compress-seq", "<compression factor>",
                    &compress_seq, 1,
                    "return likelihoods of compressed sequence"));
        config.add(new ConfigSwitch
                   ("-C", "--compress-orig", &compress_orig,
                    "compute likelihoods with compression, using same compression"
                    " factor that was used in arg-sample command (inferred from"
                    " log file)"));
        config.add(new ConfigParam<string>
                   ("-o", "--output", "<outfile>", &outfile_name,
                    "likelihood.bed",
                    "outfile for likelihoods (bed format; default=likelihood.bed)"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<SMC file>", &arg_file, "",
                    "initial ARG file (*.smc) for resampling"));
        config.add(new ConfigParam<string>
                   ("", "--region", "<start>-<end>",
                    &region, "",
                    "return likelihood for this region only"));
        config.add(new ConfigParam<int>
                   ("", "--rep", "<MCMC_rep>",
                    &mcmc_rep, 0,
                    "MCMC rep (echoed in output file)"));
        config.add(new ConfigParam<string>
                   ("", "--region-bed", "<regions.bed>", &regions_bed_file, "",
                    "(alternative to --region) Bed file containing regions to"
                    " compute likelihood"));
        config.add(new ConfigParam<string>
                   ("", "--maskmap", "<sites mask>",
                    &maskmap, "",
                    "mask map file (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--subsites", "<subsites file>", &subsites_file,
                    "file listing NAMES from sites file (or sequences from"
                    " fasta) to keep; others will not be used"));
        config.add(new ConfigParam<string>
                   ("-l", "--log-file", "<log file>", &log_file, "",
                    "log file from arg-sample run used as input to read model"
                    "parameters"));

        // model parameters
        config.add(new ConfigParamComment("Model parameters"));
        config.add(new ConfigParam<string>
                   ("-M", "--mutmap", "<mutation rate map file>", &mutmap, "",
                    "mutation map file (optional)"));
        config.add(new ConfigParam<string>
                   ("-R", "--recombmap", "<recombination rate map file>",
                    &recombmap, "",
                    "recombination map file (optional)"));
        config.add(new ConfigParam<double>
                   ("-P", "--panmictic", "<Population size>",
                    &panmictic_popsize, -1,
                    "Ignore the population model and compute ARG likelihood"
                    " under panmactic population model"));
        config.add(new ConfigParam<string>
                   ("-m", "--migfile", "<mig event file>",
                    &migfile, "",
                    "If provided, count all occurences of migration events"
                    " named in this file. File is in the format:"
                    " <mig name> <from_pop> <to_pop> <time>"
                    " and may contain multiple events. Two counts will be produced"
                    " for each event- <event_1> is number that occur, and"
                    " <event_0> is number that do not occur"));


        // misc
        config.add(new ConfigParamComment("Miscellaneous"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default=current time)"));
        config.add(new ConfigSwitch
                   ("", "--overwrite", &overwrite,
                    "overwrite output file (default: append)"));

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

        // display version info
        if (version) {
            printf(VERSION_INFO);
            return EXIT_ERROR;
        }
        return 0;
    }

    ConfigParser config;

    // input/output
    string fasta_file;
    string sites_file;
    string subsites_file;
    string outfile_name;
    string arg_file;
    string region;
    string log_file;
    string regions_bed_file;
    string migfile;

    // model parameters
    string mutmap;
    string recombmap;
    string maskmap;
    ArgModel *model;
    int mcmc_rep;
    double panmictic_popsize;
    int compress_seq;
    bool compress_orig;

    // misc
    int randseed;

    // help/information
    bool quiet;
    int verbose;
    bool version;
    bool help;
    bool help_debug;
    bool overwrite;

    FILE *outfile;
};



int get_compress_from_logfile(string log_filename) {
    FILE *infile = fopen(log_filename.c_str(), "r");
    char *line;
    while (NULL != (line = fgetline(infile))) {
        chomp(line);
        vector<string>tokens;
        split(line, ' ', tokens);
        if (strcmp(tokens[0].c_str(), "command:")==0) {
            int compress=1;
            for (unsigned int i=1; i < tokens.size(); i++) {
                if ((strcmp(tokens[i].c_str(), "-c")==0 ||
                     strcmp(tokens[i].c_str(), "--compress-seq")==0)
                    && i < tokens.size()-1) {
                    compress = atoi(tokens[i+1].c_str());
                    delete [] line;
                    fclose(infile);
                    return compress;
                }
            }
            fprintf(stderr, "Warning: did not find compress in command\n");
            delete [] line;
            fclose(infile);
            return 1;
        }
        delete [] line;
    }
    fclose(infile);
    fprintf(stderr, "Error: did not find command in logfile\n");
    assert(false);
    return 1;
}



bool parse_region(const char *region, int *start, int *end)
{
    return sscanf(region, "%d-%d", start, end) == 2;
}

//=============================================================================
// logging

/*void set_up_logging(const Config &c, int level, const char *log_mode) {
    Logger *logger;
    setLogLevel(level);
    logger = &g_logger;
    }*/

// log the program version and start time
void log_intro(int level)
{
    time_t t = time(NULL);

    printLog(level, "arg-likelihood " VERSION_TEXT "\n");
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




//=============================================================================
// statistics output

bool read_init_arg(const char *arg_file, const ArgModel *model,
                   LocalTrees *trees, vector<string> &seqnames,
                   vector<int> *invisible_recomb_pos,
                   vector<Spr> *invisible_recombs)
{
    CompressStream stream(arg_file, "r");
    if (!stream.stream) {
        printError("cannot read '%s'", arg_file);
        return false;
    }
    return read_local_trees(stream.stream, model->times, model->ntimes,
                            trees, seqnames,
                            invisible_recomb_pos, invisible_recombs);
}



// region_start and region_end are 0-based, not compressed
void print_arg_likelihood(const ArgModel *model,
                          const Sequences *sequences,
                          const LocalTrees *trees,
                          const Config *c,
                          const Region *region,
                          const vector<int> &invisible_recomb_pos,
                          const vector<Spr> &invisible_recombs,
                          const SitesMapping *sites_mapping,
                          const vector<class MigEvent> &migevents) {
    int start, end;
    if (sites_mapping != NULL) {
        start = sites_mapping->compress(region->start);
        end = sites_mapping->compress(region->end-1)+1;
    } else {
        start = region->start;
        end = region->end;
    }
    double prior = calc_arg_prior(model, trees, NULL, NULL,
                                  start, end,
                                  invisible_recomb_pos, invisible_recombs);
    double prior2 = calc_arg_prior_recomb_integrate(model, trees,
                                                    start, end);
    double like = calc_arg_likelihood(model, sequences, trees,
                                      start, end);
    int noncompat = count_noncompat(trees, sequences, start, end);
    int nrecomb=0;
    int curr_end=trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end(); ++it) {
        int curr_start = curr_end;
        curr_end += it->blocklen;
        if (curr_end <= start) continue;
        if (curr_start >= end) break;
        if (curr_end != end) nrecomb++;
    }
    fprintf(c->outfile, "%s\t%i\t%i\t%i\t%f\t%f\t%f\t%i\t%i", region->chrom.c_str(),
            region->start, region->end, c->mcmc_rep, prior, prior2, like,
            nrecomb, noncompat);

    for (unsigned int i=0; i < migevents.size(); i++) {
        int count[2]={0,0};
        int migtime = migevents[i].time;
        int from_pop = migevents[i].from_pop;
        int to_pop = migevents[i].to_pop;
        int prev_time=0; // set prev_time to time interval immediately before migtime
        for (int j=0; j < model->ntimes; j++) {
            if (model->times[j] < migtime)
                prev_time = j;
            else break;
        }
        if (start <= trees->start_coord) {
            const LocalTree *tree = trees->begin()->tree;
            for (int j=0; j < tree->nnodes; j++) {
                if (tree->nodes[j].age > prev_time) continue;
                if (j != tree->root &&
                    tree->nodes[tree->nodes[j].parent].age <= prev_time)
                    continue;
                if (model->get_pop(tree->nodes[j].pop_path, prev_time) == from_pop) {
                    if (model->get_pop(tree->nodes[j].pop_path, prev_time+1) == to_pop)
                        count[1]++;
                    else count[0]++;
                }
            }
        }
        int curr_end = trees->start_coord;
        for (LocalTrees::const_iterator it=trees->begin(); it != trees->end(); ++it) {
            int curr_start = curr_end;
            curr_end += it->blocklen;
            if (curr_end <= start) continue;
            if (curr_start >= end) break;
            const Spr *spr = &it->spr;
            if (spr->recomb_time <= prev_time && spr->coal_time >= prev_time+1 &&
                model->get_pop(spr->pop_path, prev_time) == from_pop) {
                if (model->get_pop(spr->pop_path, prev_time+1) == to_pop)
                    count[1]++;
                else count[0]++;
            }
        }
        fprintf(c->outfile, "\t%i\t%i", count[1], count[0]);

    }
    fprintf(c->outfile, "\n");

}

void parse_status_file(string log_file, int mcmc_rep, ArgModel *model) {
    char stat_filename[log_file.length()+3];
    strcpy(stat_filename, log_file.c_str());
    assert(strcmp(&stat_filename[strlen(stat_filename)-4], ".log")==0);
    sprintf(&stat_filename[strlen(stat_filename)-4], ".stats");
    FILE *stats_file;
    if (!(stats_file = fopen(stat_filename, "r"))) {
        fprintf(stderr, "cannot open %s\n", stat_filename);
        exit(0);
    }
    char *line = NULL;
    line = fgetline(stats_file);
    if (!line) {
        printError("stats file is empty");
        exit(0);
    }
    chomp(line);
    vector<string> header;
    split(line, "\t", header);
    delete [] line;
    while ((line = fgetline(stats_file))) {
        vector<string> tokens;
        split(line, "\t", tokens);
        assert(tokens.size() == header.size());
        if (tokens[0] != "resample") continue;
        int iter;
        if (sscanf(tokens[1].c_str(), "%d", &iter) != 1) {
            printError("Error getting iter from stat file\n");
            assert(0);
        }
        if (iter != mcmc_rep) continue;
        model->init_params_from_statfile(header, line);
        fclose(stats_file);
        return;
    }
    printError("Did not find rep %i in stats file\n", mcmc_rep);
    exit(0);
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

    // setup logging
    //    set_up_logging(c, c.verbose, "w");

    // log intro
    log_prog_commands(LOG_LOW, argc, argv);
    Timer timer;


    // init random number generator
    // probably never used in this program
    if (c.randseed == 0)
        c.randseed = time(NULL);
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n", c.randseed);

    vector<class MigEvent> migevents;
    if (c.migfile != "") {
        FILE *infile = fopen(c.migfile.c_str(), "r");
        if (infile == NULL) {
            printError("Could not open %s", c.migfile.c_str());
            return EXIT_ERROR;
        }
        char tmp[1000];
        int from_pop, to_pop, mig_time;
        while (EOF != fscanf(infile, "%s %i %i %i",
                             tmp, &from_pop, &to_pop, &mig_time)) {
            migevents.push_back(MigEvent(tmp, from_pop, to_pop,
                                         mig_time));
        }
        fclose(infile);
        printf("Counting %i migevents", (int)migevents.size());
    }

    // read sequences
    Sites sites;
    Sequences sequences;
    Region seq_region;

    if (c.sites_file == "") {
        printError("Error: Require --sites");
        return EXIT_ERROR;
    }
    // read sites file
    // read sites
    CompressStream stream(c.sites_file.c_str());
    if (!stream.stream ||
        !read_sites(stream.stream, &sites)) {
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

    if (sites.get_num_sites() == 0) {
        fprintf(stderr, "Error: No sites\n");
        return EXIT_ERROR;
    }
    SitesMapping *sites_mapping = new SitesMapping();
    if (!find_compress_cols(&sites, c.compress_seq, sites_mapping)) {
        printError("unable to compress sequences at given compression level"
                   " (--compress-seq)");
        return EXIT_ERROR;
    }


    // compress sequences
    compress_sites(&sites, sites_mapping);
    make_sequences_from_sites(&sites, &sequences);

    // read in and compress mask
    if (c.maskmap != "") {
        TrackNullValue maskmap;
        CompressStream stream(c.maskmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &maskmap, seq_region)) {
            printError("cannot read mask map '%s'",
                       c.maskmap.c_str());
            return EXIT_ERROR;
        }

        // apply mask
        // TODO: when mask is compressed, only allow it to apply to
        // whole compressed regions
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
    if (c.log_file == "") {
        printError("--log-file is required\n");
        return EXIT_ERROR;
    }
    c.model = new ArgModel(c.log_file.c_str());
    parse_status_file(c.log_file, c.mcmc_rep, c.model);
    if (c.panmictic_popsize > 0) {
        if (c.model->pop_tree != NULL) {
            delete(c.model->pop_tree);
            c.model->pop_tree = NULL;
        }
        c.model->free_popsizes();
        c.model->set_popsizes(c.panmictic_popsize);
    }
    // log model has compressed rates so un-compress them; they
    // will later be re-compressed
    int old_compress = get_compress_from_logfile(c.log_file);
    if (c.compress_orig) {
        if (c.compress_seq != 1) {
            printError("ERROR: --compress-seq and --compress-orig should not"
                       " both be used");
            return EXIT_ERROR;
        }
        printLog(LOG_LOW, "Using compression factor %i\n", c.compress_seq);
        c.compress_seq = old_compress;
    }
    c.model->rho /= (double)old_compress;
    c.model->mu /= (double)old_compress;

    // read model parameter maps if given
    if (c.mutmap != "") {
        CompressStream stream(c.mutmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model->mutmap, seq_region)) {
            printError("cannot read mutation rate map '%s'", c.mutmap.c_str());
            return EXIT_ERROR;
        }
    }
    if (c.recombmap != "") {
        CompressStream stream(c.recombmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model->recombmap, seq_region)) {
            printError("cannot read recombination rate map '%s'",
                       c.recombmap.c_str());
            return EXIT_ERROR;
        }
    }


    compress_model(c.model, sites_mapping, c.compress_seq);

    c.model->log_model();

    // setup init ARG
    LocalTrees *trees = NULL;
    unique_ptr<LocalTrees> trees_ptr;
    if (c.arg_file == "") {
        printError("Error: --arg-file required\n");
        return EXIT_ERROR;
    }
    // init ARG from file
    trees = new LocalTrees();
    trees_ptr = unique_ptr<LocalTrees>(trees);
    vector<string> seqnames;
    vector<int> invisible_recomb_pos;
    vector<Spr> invisible_recombs;
    if (!read_init_arg(c.arg_file.c_str(), c.model, trees, seqnames,
                       &invisible_recomb_pos, &invisible_recombs)) {
        printError("could not read ARG");
        return EXIT_ERROR;
    }
    if (!trees->set_seqids(seqnames, sequences.names)) {
        printError("input ARG's sequence names do not match input"
                   " sequences");
        return EXIT_ERROR;
    }
    if (c.panmictic_popsize > 0)
        remove_population_paths(trees);
    if (c.region == "" && c.regions_bed_file == "") {
        char tmp[1000];
        sprintf(tmp, "%i-%i", trees->start_coord+1, trees->end_coord);
        c.region = string(tmp);
    }
    if (sites_mapping) {
        compress_local_trees(trees, sites_mapping);
        for (unsigned int i=0; i < invisible_recomb_pos.size(); i++)
            invisible_recomb_pos[i] = sites_mapping->compress(invisible_recomb_pos[i], 0,
                                                              i==0 ? 0 : invisible_recomb_pos[i-1]);
    }

    printLog(LOG_LOW, "read input ARG (chrom=%s, start=%d, end=%d,"
             " nseqs=%d)\n",
             trees->chrom.c_str(), trees->start_coord, trees->end_coord,
             trees->get_num_leaves());

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);


    if (!(c.outfile = fopen(c.outfile_name.c_str(),
                            c.overwrite ? "w" : "a"))) {
        printError("Could not open out file %s for writing\n",
                   c.outfile_name.c_str());
        return EXIT_ERROR;
    }

    if (c.overwrite) {
        fprintf(c.outfile, "#chrom\tstart\tend\trep\tprior\tprior2\tlikelihood\tnrecomb\tncompat");
        for (unsigned int i=0; i < migevents.size(); i++)
            fprintf(c.outfile, "\t%s_1\t%s_0", migevents[i].name.c_str(),
                    migevents[i].name.c_str());
        fprintf(c.outfile, "\n");
    }

    // get likelihod
    printLog(LOG_LOW, "\n");

    if (c.region != "") {
        int start, end;
        if (!parse_region(c.region.c_str(), &start, &end)) {
            printError("Error parsing region string %s\n", c.region.c_str());
            return EXIT_ERROR;
        }
        seq_region.start = start-1;
        seq_region.end = end;

        print_arg_likelihood(c.model, &sequences, trees, &c,
                             &seq_region,
                             invisible_recomb_pos, invisible_recombs,
                             sites_mapping, migevents);
    } else {
        if (c.regions_bed_file == "") {
            printError("Need to supply --region or --region-bed");
            return EXIT_ERROR;
        }
        FILE *bedfile;
        if (!(bedfile = fopen(c.regions_bed_file.c_str(), "r"))) {
            printError("Error opening bed file %s\n", c.regions_bed_file.c_str());
            return EXIT_ERROR;
        }
        char *line;
        while (NULL != (line = fgetline(bedfile))) {
            vector<string> tokens;
            if (line[0] == '#') {
                delete [] line;
                continue;
            }
            split(line, "\t", tokens);
            if (tokens.size() < 3) {
                printError("Should have at least three entries in each line of bed file");
                return EXIT_ERROR;
            }
            if (tokens[0] == sites.chrom) {
                seq_region.start = atoi(tokens[1].c_str());
                seq_region.end = atoi(tokens[2].c_str());
                print_arg_likelihood(c.model, &sequences, trees, &c,
                                     &seq_region,
                                     invisible_recomb_pos, invisible_recombs,
                                     sites_mapping, migevents);
            }
            delete [] line;
        }
        fclose(bedfile);
    }
    // TODO: call in loop for bed file

    // final log message
    maxrss = get_max_memory_usage() / 1000.0;
    printTimerLog(timer, LOG_LOW, "sampling time: ");
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);
    printLog(LOG_LOW, "FINISH\n");

    // clean up
    fclose(c.outfile);
    return 0;
}
