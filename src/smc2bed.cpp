#include "getopt.h"
#include <iostream>
#include <fstream>
#include <assert.h>

// argweaver includes
#include "argweaver/local_tree.h"
#include "argweaver/compress.h"
#include "argweaver/parsing.h"
#include "argweaver/model.h"

//using namespace spidir;
using namespace argweaver;

void print_usage() {
    printf("smc2bed: This program converts a single smc file into a\n"
           "  bed file. The bed file format is chrom,start,end,sample,tree.\n"
           "  The tree nodes are labelled with NHX-style comments indicating\n"
           "  the nodes and times of the recombination event which leads to\n"
           "  the next tree.\n\n"
           "This program is intended for use combining multiple SMC files\n"
           "  from different MCMC samples; the pipeline for doing this is to\n"
           "  run smc2bed on each file, piping the results to sort-bed, then\n"
           "  bgzip. The resulting file can be indexed using tabix.\n\n");
    printf("Usage: ./smc2bed [OPTIONS] <smc-file>\n"
           "  smc-file can be gzipped\n"
           " OPTIONS:\n"
           " --region START-END\n"
           "   Process only these coordinates (1-based)\n"
           " --sample <sample>\n"
           "   Give the sample number for this file; this is important\n"
           "   when combining multiple smc files."
           " --log-file <file.log>\n"
           "   Log file from arg-sample run; this is used as input to read model"
           "   parameters. If not provided, smc2bed will look for log file"
           "   in directory with smc file.\n");
}


bool guess_log_file(char *smc_file, char *log_file) {
    int len = strlen(smc_file);
    strcpy(log_file, smc_file);
    if (strcmp(&smc_file[len-7], ".smc.gz")==0) {
        int pos=len-8;
        while (pos >= 0 && smc_file[pos] != '.') pos--;
        if (pos < 0) return false;
        log_file[pos]='\0';
        char tmp_log_file[1000];
        strncpy(tmp_log_file, log_file, 999);
        sprintf(log_file, "%s%s", tmp_log_file, ".log");
        return true;
    }
    return false;
}


int main(int argc, char *argv[]) {
    char c;
    int region[2]={-1,-1};
    LocalTrees *trees=NULL;
    char *log_file = NULL;
    ArgModel *model;
    int sample=0, opt_idx;
    struct option long_opts[] = {
        {"region", 1, 0, 'r'},
        {"sample", 1, 0, 's'},
        {"log-file", 1, 0, 'l'},
        {"help", 0, 0, 'h'},
        {0,0,0,0}};
    while ((c = (char)getopt_long(argc, argv, "r:s:l:h", long_opts, &opt_idx))
           != -1) {
        switch (c) {
        case 'r':
            if (2 != (sscanf(optarg, "%d-%d", &region[0], &region[1]))) {
                fprintf(stderr, "error parsing region %s\n", optarg);
                return 1;
            }
            region[0]--;  //convert to 0-based
            break;
        case 's':
            sample = atoi(optarg);
            break;
        case 'l':
            log_file = optarg;
            break;
        case 'h':
            print_usage();
            return 0;
        case '?':
            fprintf(stderr, "unknown option. Try --help\n");
            return 1;
        }
    }
    if (optind != argc - 1) {
        fprintf(stderr, "Bad arguments. Try --help\n");
        return 1;
    }
    Logger *logger = new Logger(stderr, LOG_HIGH);
    g_logger.setChain(logger);

    if (log_file == NULL) {
        log_file = (char*)malloc((strlen(argv[optind])+10)*sizeof(char));
        if (!guess_log_file(argv[optind], log_file)) {
            fprintf(stderr, "Could not guess log file name, provide with -l");
            return 1;
        }
    }

    model = NULL;
    if (log_file != NULL) {
        model = new ArgModel(log_file);
    }

    CompressStream instream(argv[optind], "r");
    vector<string> seqnames;

    trees = new LocalTrees();
    if (!read_local_trees(instream.stream, model->times, model->ntimes,
                          trees, seqnames)) {
        fprintf(stderr, "Error parsing SMC file\n");
        return 1;
    }
    if (region[0] != -1) {
        LocalTrees *trees2 = partition_local_trees(trees, region[0], true);
        delete trees;
        trees = trees2;
    }
    if (region[1] != -1)
        partition_local_trees(trees, region[1], true);
    write_local_trees_as_bed(stdout, trees, seqnames,
                             model, sample);
    instream.close();
    return 0;
}
