// C/C++ includes
#include <time.h>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>
#include <vector>
#include <math.h>
#include <queue>
#include <set>
#include <map>

// argweaver includes
#include "argweaver/ConfigParam.h"
#include "argweaver/logging.h"
#include "argweaver/parsing.h"
#include "argweaver/track.h"
#include "argweaver/Tree.h"
#include "argweaver/tabix.h"
#include "argweaver/compress.h"
#include "argweaver/IntervalIterator.h"
#include "argweaver/model.h"
#include "argweaver/seq.h"
//#include "allele_age.h"


using namespace argweaver;
using namespace spidir;

/* Initial version: just have it output bedgraph with a stat for each
   line of input.
*/

#define VERSION_INFO "arg-summarize 0.3"
bool html;
int summarize=0;
int getNumSample=0;
int getMean=0;
int getStdev=0;
int getQuantiles=0;
vector <double> quantiles;
vector<string> node_dist_leaf1;
vector<string> node_dist_leaf2;
vector<string> min_coal_time_ind1;
vector<string> min_coal_time_ind2;
bool quiet;
vector<string> ind_dist_leaf1;
vector<string> ind_dist_leaf2;
set<string> cluster_group;

const int EXIT_ERROR = 1;
const int POPMODEL_OPT = 3;
const int EXPERIMENTAL_OPT = 4;


class MigStat {
public:
    MigStat(string name, int p0[2], double dt, const ArgModel *model,
            const string hap="") : name(name),hap(hap){
        bool found=false;
        // ensure that t[0] < t[1]
        p[0] = p0[0];
        p[1] = p0[1];
        for (int i=0; i < model->ntimes-1; i++)
            if (model->times[i] < dt && model->times[i+1] > dt) {
                t[0] = i;
                t[1] = i+1;
                found=true;
            }
        if (!found) {
            fprintf(stderr, "Error finding time intervals spanning mig time %f\n", dt);
            exit(1);
        }
    }
    string name;
    int p[2];
    int t[2];
    string hap;
};


/* class of miscellaenous data structures to be passed around arg-summarize
   functions */
class ArgSummarizeData {
public:
    ArgModel *model;
    vector< set<string> > group;
    vector<string> spr_leaf;

    //individuals to test coal into various groups; maps an individual
    // name to one or two haploid chromosome names
    map<string, set<string> > coalgroup_inds;

    //maps each haploid chromosome to a group identifier
    // does not have to map all chromosomes, only the ones with groups
    map<string, int> coalgroups;

    //maps integer from coalgroups to a group name
    vector<string> coalgroup_names;

    vector<MigStat> migstat;
};


class Config
{
public:
    Config()
    {
        sample_num=0;
        make_parser();
    }
    void make_parser() {
        config.clear();

        config.add(new ConfigParamComment("Input files"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg-file", "<file.bed.gz>", &argfile,
                    "Bed file containing args sampled by ARGweaver. Should"
                    " be created with smc2bed and sorted with sort-bed. If"
                    " using --region or --bedfile, also needs to be gzipped"
                    " and tabix'd"));
        config.add(new ConfigParam<string>
                   ("-r", "--region", "<chr:start-end>", &region,
                    "region to retrieve statistics from (1-based coords)"));
        config.add(new ConfigParam<string>
                   ("-b", "--bed-file", "<file.bed>", &bedfile,
                    "regions to retrieve statistics from (alternative to "
                    "region)"));
        config.add(new ConfigParam<string>
                   ("-s", "--subset", "<hap_list.txt>", &hapfile,
                    "file with list of leafs to keep (rest will be pruned)"));
        config.add(new ConfigParam<string>
                   ("", "--subset-inds", "<ind_list.txt>", &indfile,
                    "(Alternative to --subset) Subset ARG to only contain"
                    " individuals listed in file."
                    " Assume each individual is diploid with haploid names"
                    " ind_1 and ind_2"));
        config.add(new ConfigParam<string>
                   ("-f", "--snp-file", "<snp_file.bed>", &snpfile,
                    "Compute statistics for specific SNPs. Each output row will"
                    " contain information about each SNP (alleles and"
                    " frequencies), and statistics will be computed in two ways-"
                    " once across all MCMC samples, and once only for samples"
                    " where the site pattern for the SNP agrees with an infinite"
                    " sites model."
                    " The SNP file should have a header like this:\n"
                    " #NAMES NA06985_1       NA06985_2       NA06994_1 ...\n"
                    " and then each line should be tab-delimited with the format"
                    " chr,start,end,AAAACAAAAA, where the last column gives the"
                    " alleles for each haplotype in the order indicated in the"
                    " header. The file needs to be sorted and indexed with"
                    " tabix."));
        config.add(new ConfigParam<int>
                   ("-p", "--sample", "<sample_num>", &sample_num,
                    "If provided, only process trees in argfile with given MCMC"
			   " sample number"));
        config.add(new ConfigParam<string>
                   ("-l", "--log-file", "<file.log>", &logfile,
                    "Log file from arg-sample run, used to obtain model"
                    " parameters. Required in combination with several options."));
        config.add(new ConfigParam<string>
                  ("-m", "--mig-file", "<migfile.txt>", &migfile,
                   "Report statistics on migrations listed in this file. Format is:\n"
                   " statName pop1 pop2 t\n"
                   " will return a boolean statistic named statName indicating existance"
                   " of a lineage that is in pop1 at the time interval immediately"
                   " after time t and pop2 at the time interval before time t2"
                   " (indicating move from pop1 to pop2 looking backwards in time)",
		   POPMODEL_OPT));
        config.add(new ConfigParam<string>
                   ("", "--hap-mig-file", "<migfile.txt>", &hapmigfile,
                    "Similar to --mig-file, but migfile.txt has five columns:\n"
                    " statName pop1 pop2 t hap\n"
                    " where hap is one of the haploid lineages of the ARG, will"
                    " report whether the lineage whose descendent is hap goes from\n"
                    " pop1 to pop2 at time t (looking backwards)",
		    POPMODEL_OPT));

        config.add(new ConfigParamComment("Statistics to retrieve"));
        config.add(new ConfigSwitch
                   ("-E", "--tree", &rawtrees,
                    "output newick tree strings (cannot use summary options"
                    " with this)"));

        config.add(new ConfigSwitch
                   ("-T", "--tmrca", &tmrca,
                    "time to the most recent common ancestor"));
        config.add(new ConfigSwitch
                   ("-p", "--pi", &pi,
                    "average distance between two lineages in tree"));
        config.add(new ConfigSwitch
                   ("-B", "--branchlen", &branchlen, "total branch length"));
        config.add(new ConfigSwitch
                   ("-R", "--recomb", &recomb,
                    "recombination rate per generation/bp (NOTE this is not"
                    " particularly well estimated by the SMC)"));
        config.add(new ConfigSwitch
                   ("", "--recombs-per-time", &recombs_per_time,
                    "number of recombinations for each time interval"
                    " (requires --log-file)"));
        config.add(new ConfigSwitch
                   ("", "--invis-recombs-per-time", &invis_recombs_per_time,
                    "number of invisiblel recombinations for each time interval"
                    " (requires --log-file)", EXPERIMENTAL_OPT));
        config.add(new ConfigSwitch
                   ("", "--branchlen-per-time", &branchlen_per_time,
                   "total branchlen existing at each time interval"
                    " (does not count root branch; requires --log-file)"));
        config.add(new ConfigSwitch
                   ("-K", "--breaks", &breaks,
                    "recombinations per bp (not normalized by tree size"));
        config.add(new ConfigSwitch
                   ("-H", "--tmrca-half", &tmrca_half,
                    "time for half of samples to reach a common ancestor"));
        config.add(new ConfigSwitch
                   ("-F", "--rth", &rth,
                    "relative TMRCA Halftime (tmrca_half/tmrca)"));
        config.add(new ConfigSwitch
                   ("-P", "--popsize", &popsize,
                    "popsize estimated by coal rates in local tree"
                    " (not recommended for estimating population sizes)"));
        config.add(new ConfigSwitch
                   ("", "--max-coal-rate", &max_coal_rate,
                    "Maximum rate of coalescence in tree (measured in fraction"
                    " of branches coalescing per generation, take max over all pairs"
                    " of discrete time points", EXPERIMENTAL_OPT));
        config.add(new ConfigSwitch
                   ("-A", "--allele-age", &allele_age,
                    "(requires --snp). Compute allele age for SNPs. Only works"
                    " for pre-phased SNPs"));
        config.add(new ConfigSwitch
                   ("", "--min-allele-age", &min_allele_age,
                    "(requires --snp). Compute minimum possible age for SNP"));
        config.add(new ConfigParam<string>
                   ("-D", "--node-dist", "<leaf0;leaf1,leaf2;...>",
                    &node_dist,
                    "leaf branch length (for single leaf) OR distance between"
                    " pairs of leafs. Multiple statistics can be requested by"
                    " separating arguments with semi-colon. In this example"
                    " return length of leaf leaf0, then distance between leaf1 and leaf2"));
        config.add(new ConfigSwitch
                   ("", "--node-dist-all", &node_dist_all,
                    "return distance between all pairs of leaf nodes."
                    "Requires --subset option (though this may specify a file"
                    " with all leaves)"));
        /*	config.add(new ConfigParam<string>
		   ("-I", "--ind-dist", "<ind1_hap1,ind1_hap2;ind2_hap1,ind2_hap2,...>",
		    &ind_dist,
		    "return three statistics for each individual (identified"
		    " by two haplotypes): the minimum coalescence time to rest of"
		    " tree, the maximum coalescence time to rest of tree, and"
		    " a 0/1 indicating whether the two lineages coalesce with each"
		    " other before any other lineage")); */
        config.add(new ConfigParam<string>
                   ("", "--min-coal-time", "<ind1,ind2;...>",
                    &min_coal_time,
                    "return the most recent coalescence time between two"
                    " diploid individuals, whose lineages are named ind1_1,"
                    " ind1_2, ind2_1, ind2_2. May specify multiple pairs of"
                    " individuals, will return the min coal time for each pair"));
        config.add(new ConfigSwitch
                   ("", "--min-coal-time-all", &min_coal_time_all,
                    "Same as min-coal-time, but for all pairs of individuals."
                    " Each individual must have two linages named with the "
                    " convention ind_1, ind_2. Requires --subset-inds argument"
                    " (though this may specify a file with all individuals"));
        config.add(new ConfigSwitch
                   ("-Z", "--zero-len", &zero,
                    "number of branches of length zero"));
        config.add(new ConfigSwitch
                   ("-C", "--coalcounts", &coalcounts,
                    "number of coal events at each discretized time point"
                    " (requires --log-file)"));
        config.add(new ConfigSwitch
                   ("", "--coalcounts-cluster", &coalcounts_cluster,
                    "maximum number of coalescences clustered together at each time point"
                    " (requires --log-file)"));
        config.add(new ConfigParam<string>
                   ("-G", "--group", "<group_file1,groupfile2,...>", &groupfile,
                    "Output boolean indicating whether individuals in each group"
                    " file cluster together in the local tree."));
        /*        config.add(new ConfigParam<string>
                   ("", "--group-no-root", "<group_noroot_file>",
                    &group_noroot_file,
                    "Like --group, but consider grouping in unrooted tree"));*/
        /*        config.add(new ConfigParam<string>
                   ("", "--coal-groups", "<coal_group_file>", &coalgroup_file,
                    "(For use with --coal-group-inds); Output boolean indicating"
                    " whether individuals in coal-group-inds coalesce with"
                    " each group in coal_group_file. Coal_group_file should"
                    " be a two column file with haploid chromosome name and"
                    " group name"));
        config.add(new ConfigParam<string>
                   ("", "--coal-group-inds", "<coal_group_inds>",
                    &coalgroup_inds_file,
                    "(For use with --coal-groups) File containing individuals"
                    " to test coalescence with groups defined in coal_group_file"
                    " coal_group_inds file should contain one individual per"
                    " line, which could be one or two haploid chromosome"
                    " names"));
        config.add(new ConfigParam<string>
                   ("", "--spr-leaf", "<leaf1,leaf2,...>",
                    &spr_leaf_names,
                    "For each leaf, indicate binary whether the SPR operation"
                    "at the end of segment affects the branch coming from that"
                    "leaf"));*/
        config.add(new ConfigParam<string>
                   ("", "--cluster-test", "<group_file>", &cluster_group_file,
                    "Return maximum likelihood ratio statistic"
                    " relating to whether lineages in group_file are randomly"
                    " disperesed between the children, as well as the time of the"
                    " most significant node (maximum taken over nodes)",
                    EXPERIMENTAL_OPT));
        /*                    "Return minimum number of pruning operations required to "
                    " reduce each local tree to the lineages named in <group_file>."
                    " (Lower numbers indicate clustering)"));*/
        config.add(new ConfigSwitch
                   ("-N", "--numsample", &numsample,
                    "number of MCMC samples covering each region"));

        config.add(new ConfigParamComment("Summary options (if not given,"
                                          " statistics will be output for each"
                                          " MCMC sample)"));
        config.add(new ConfigSwitch
                   ("-M", "--mean", &mean,
                    "return mean across all MCMC samples"));
        config.add(new ConfigSwitch
                   ("-S", "--stdev", &stdev,
                    "return standard deviation across all MCMC samples"));
        config.add(new ConfigParam<string>
                   ("-Q", "--quantile", "<q1,q2,q3,...>", &quantile,
                    "return the requested quantiles for each samples"));

        config.add(new ConfigParamComment("Misceallaneous"));
        config.add(new ConfigParam<int>
                   ("-u", "--burnin", "<num>", &burnin, 0,
                    "Discard results from iterations < burnin before computing statistics"));
        config.add(new ConfigSwitch
                   ("-n", "--no-header", &noheader, "Do not output header"));
        config.add(new ConfigParam<string>
                   ("-t", "--tabix-dir", "<tabix dir>", &tabix_dir,
                    "Specify the directory of the tabix executable"));
        config.add(new ConfigSwitch
                   ("", "--html", &html,
                    "output HTML instead of plain text (useful with --tree;"
                    " makes link to GIF image)", EXPERIMENTAL_OPT));
        config.add(new ConfigSwitch
                   ("-q", "--quiet", &quiet, "Proceed quietly"));
        config.add(new ConfigSwitch
                   ("-v", "--version", &version, "display version information"));
        config.add(new ConfigSwitch
                   ("-h", "--help", &help, "display help information"));
	config.add(new ConfigSwitch
		   ("-d", "--help-popmodel", &help_popmodel, "display help information specific to ARGweaver-D (population model)"));
        config.add(new ConfigSwitch
                   ("", "--help-advanced", &help_advanced,
                    "display help information for experimental/advanced options"));
    }

    int parse_args(int argc, char **argv) {
        if (!config.parse(argc, (const char**) argv)) {
            if (argc < 2)
                config.printHelp();
            return EXIT_ERROR;
        }
        if (help) {
            config.printHelp();
            return EXIT_ERROR;
        }
	if (help_popmodel) {
            config.printHelp(stderr, POPMODEL_OPT);
            return EXIT_ERROR;
        }
        if (help_advanced) {
            config.printHelp(stderr, EXPERIMENTAL_OPT);
            return EXIT_ERROR;
        }
        if (version) {
            printf(VERSION_INFO);
            return EXIT_ERROR;
        }
        return 0;
    }
    ConfigParser config;

    string argfile;
    string region;
    string bedfile;
    string hapfile;
    string indfile;
    string snpfile;
    string logfile;
    string groupfile;
    string coalgroup_inds_file;
    string coalgroup_file;
    string spr_leaf_names;
    string migfile;
    string hapmigfile;
    int sample_num;

    bool rawtrees;
    bool html;
    bool tmrca;
    bool pi;
    bool branchlen;
    bool recomb;
    bool recombs_per_time;
    bool invis_recombs_per_time;
    bool branchlen_per_time;
    bool breaks;
    bool tmrca_half;
    bool rth;
    bool popsize;
    bool max_coal_rate;
    bool allele_age;
    bool min_allele_age;
    string node_dist;
    bool node_dist_all;
    string min_coal_time;
    bool min_coal_time_all;
    string cluster_group_file;
    string ind_dist;
    bool zero;
    bool coalcounts;
    bool coalcounts_cluster;
    bool numsample;

    bool mean;
    bool stdev;
    string quantile;

    int burnin;
    bool noheader;
    string tabix_dir;
    bool quiet;
    bool version;
    bool help;
    bool help_popmodel;
    bool help_advanced;
};

void checkResults(IntervalIterator<vector<double> > *results) {
    Interval<vector<double> > summary=results->next();
    vector<vector <double> > scores;
    while (summary.start != summary.end) {
        scores = summary.get_scores();
        if (scores.size() > 0) {
            if (html) printf("<tr><td>\n");
            printf("%s\t", summary.chrom.c_str());
            if (html) printf("</td><td>");
            printf("%i\t", summary.start);
            if (html) printf("</td><td>");
            printf("%i", summary.end);
            vector<double> tmpScore(scores.size());
            int numscore = scores[0].size();
            assert(numscore > 0);
            for (int i=0; i < numscore; i++) {
                int have_mean = 0;
                double meanval=-1;
                for (unsigned int j=0; j < scores.size(); j++)
                    tmpScore[j] = scores[j][i];
                if (i==0 && getNumSample > 0) {
                    if (html) printf("</td><td>");
                    printf("\t%i", (int)scores.size());
                }
                for (int j=1; j <= summarize; j++) {
                    if (getMean==j) {
                        meanval = compute_mean(tmpScore);
                        have_mean=1;
                        if (html) printf("</td><td>");
                        printf("\t%g", meanval);
                    } else if (getStdev==j) {
                        if (!have_mean)
                            meanval = compute_mean(tmpScore);
                        if (html) printf("</td><td>");
                        printf("\t%g", compute_stdev(tmpScore, meanval));
                    } else if (getQuantiles==j) {
                        vector<double> q =
                            compute_quantiles(tmpScore, quantiles);
                        for (unsigned int k=0; k < quantiles.size(); k++) {
                        if (html) printf("</td><td>");
                            printf("\t%g", q[k]);
                        }
                    }
                }
            }
            printf("\n");
            if (html) printf("</td></tr>\n");
        }
        summary = results->next();
    }
}

class BedLine {
public:
    BedLine(char *chr, int start, int end, int sample, char *nwk,
            SprPruned *trees=NULL) :
        start(start), end(end), sample(sample),
        trees(trees) {
        chrom = new char[strlen(chr)+1];
        strcpy(chrom, chr);
        if (nwk != NULL) {
            newick = (char*)malloc((strlen(nwk)+1)*sizeof(char));
            strcpy(newick, nwk);
        } else newick=NULL;
    };
    ~BedLine() {
        //      if (orig_tree != NULL) delete orig_tree;
        //      if (pruned_tree != NULL) delete pruned_tree;
        delete [] chrom;
        if (newick != NULL)
            free(newick);
        //        delete trees;
    }
    char *chrom;
    int start;
    int end;
    int sample;
    SprPruned *trees;
    char *newick;
    vector<double> stats;
    char derAllele, otherAllele;
    int derFreq, otherFreq;
    int infSites;
};


void scoreBedLine(BedLine *line, vector<string> &statname,
                  ArgSummarizeData &data,
                  double allele_age=-1, double min_allele_age = -1,
                  int infsites=-1) {
    Tree * tree = (line->trees->pruned_tree != NULL ?
                   line->trees->pruned_tree :
                   line->trees->orig_tree);
    double bl=-1.0;
    int node_dist_idx=0;
    int min_coal_time_idx=0;
    const ArgModel *model = data.model;
    int ind_dist_idx=0;
    if (line->stats.size() == statname.size()) return;
    line->stats.resize(statname.size());
    for (unsigned int i=0; i < statname.size(); i++) {
        if (statname[i] == "tmrca")
            line->stats[i] = tree->tmrca();
        else if (statname[i]=="tmrca_half")
            line->stats[i] = tree->tmrca_half();
        else if (statname[i]=="pi")
            line->stats[i] = tree->avg_pairwise_distance();
        else if (statname[i]=="branchlen") {
            if (bl < 0) {
                line->stats[i] = tree->total_branchlength();
                bl=line->stats[i];
            }
        }
        else if (statname[i]=="rth")
            line->stats[i] = tree->rth();
        else if (statname[i]=="popsize")
            line->stats[i] = tree->popsize();
        else if (statname[i]=="recomb") {
            if (bl < 0) bl = tree->total_branchlength();
            line->stats[i] = 1.0/(bl*(double)(line->end - line->start));
        }
        else if (statname[i]=="breaks") {
            line->stats[i] = 1.0/((double)(line->end - line->start));
        }
        else if (statname[i]=="zero_len") {
            line->stats[i] = tree->num_zero_branches();
        }
        else if (statname[i]=="max_coal_rate") {
            line->stats[i] = tree->maxCoalRate(model);
        }
        else if (statname[i]=="tree") {
            if (line->trees->pruned_tree != NULL) {
                string tmp =
                    line->trees->pruned_tree->format_newick(false, true, 1,
                                                     &line->trees->pruned_spr);
                //pruned tree will be fewer characters than whole tree
                sprintf(line->newick, "%s", tmp.c_str());
            }
        }
        else if (statname[i]=="allele_age")
            line->stats[i] = allele_age;
        else if (statname[i]=="min_allele_age")
            line->stats[i] = min_allele_age;
        else if (statname[i]=="inf_sites")
            line->stats[i] = (double)infsites;
        else if (statname[i].substr(0, 9)=="node_dist") {
	    if (node_dist_leaf2[node_dist_idx].empty()) {
		line->stats[i] = tree->branch_len(node_dist_leaf1[node_dist_idx]);
	    } else {
		line->stats[i] =
		    tree->distBetweenLeaves(node_dist_leaf1[node_dist_idx],
					    node_dist_leaf2[node_dist_idx]);
	    }
            node_dist_idx++;
        }
        else if (statname[i].substr(0,13)=="min_coal_time") {
            line->stats[i] = tree->minCoalBetweenInds(min_coal_time_ind1[min_coal_time_idx],
                                                      min_coal_time_ind2[min_coal_time_idx]);
            min_coal_time_idx++;
        }
        else if (statname[i].substr(0, 8)=="recombs.") {
            const NodeSpr *nodespr = (line->trees->pruned_tree != NULL ?
                                      &(line->trees->pruned_spr) :
                                      &(line->trees->orig_spr));
            for (int j=0; j < model->ntimes; j++) {
                assert(i+j < statname.size() &&
                       statname[i+j].substr(0, 8)=="recombs.");
                line->stats[i+j] = 0;
            }
            if (nodespr->recomb_node != NULL) {
                int t = model->discretize_time(nodespr->recomb_time);
                line->stats[i + t] = 1;
            }
            i += model->ntimes - 1;
        }
        else if (statname[i].substr(0, 14)=="invis-recombs.") {
            const NodeSpr *nodespr = (line->trees->pruned_tree != NULL ?
                                      &(line->trees->pruned_spr) :
                                      &(line->trees->orig_spr));
            for (int j=0; j < model->ntimes; j++) {
                assert(i+j < statname.size() &&
                       statname[i+j].substr(0, 14)=="invis-recombs.");
                line->stats[i+j] = 0;
            }
            if (nodespr->is_invisible()) {
                int t = model->discretize_time(nodespr->recomb_time);
                line->stats[i + t] = 1;
            }
            i += model->ntimes - 1;
        }
        else if (statname[i].substr(0, 10)=="branchlen.") {
            for (int j=0; j < model->ntimes; j++) {
                assert(i+j < statname.size() &&
                       statname[i+j].substr(0, 10)=="branchlen.");
                line->stats[i+j] = 0;
            }
            for (int j=0; j < tree->nnodes; j++) {
                if (tree->nodes[j] == tree->root) continue;
                int age1 = model->discretize_time(tree->nodes[j]->age);
                int age2 = model->discretize_time(tree->nodes[j]->parent->age);
                for (int k=age1; k < age2; k++) {
                    line->stats[i + k] += (model->times[k + 1] - model->times[k]);
                }
            }
            i += model->ntimes - 1;
        }
	else if (statname[i].substr(0, 8)=="ind_dist") {
	    string h1 = ind_dist_leaf1[ind_dist_idx];
	    string h2 = ind_dist_leaf2[ind_dist_idx];
	    Node *parent = tree->are_sisters(h1, h2);
	    if (parent != NULL) {
		line->stats[i] = line->stats[i+1] = tree->branch_len(h1) + parent->dist;
		line->stats[i+2] = 1;
	    } else {
		double d1 = tree->branch_len(h1);
		double d2 = tree->branch_len(h2);
		line->stats[i] = min(d1, d2);
		line->stats[i+1] = max(d1, d2);
		line->stats[i+2] = 0;
	    }
	    i += 2;
	    ind_dist_idx++;
	}
        else if (statname[i].substr(0, 11)=="coalcounts.") {
            vector<double>coal_counts = tree->coalCounts(model->times, model->ntimes);
            for (unsigned int j=0; j < coal_counts.size(); j++) {
                assert(i+j < statname.size() &&
                       statname[i+j].substr(0,11)=="coalcounts.");
                line->stats[i+j] = coal_counts[j];
            }
            i += coal_counts.size() - 1;
        } else if (statname[i].substr(0, 19)=="coalcounts-cluster.") {
            vector<int>coal_counts = tree->coalCountsCluster(model->times, model->ntimes);
            for (unsigned int j=0; j < coal_counts.size(); j++) {
                assert(i+j < statname.size() &&
                       statname[i+j].substr(0,19)=="coalcounts-cluster.");
                line->stats[i+j] = (double)coal_counts[j];
            }
            i += coal_counts.size() - 1;
        }
        else if (statname[i].substr(0, 5)=="group") {
            for (unsigned int j=0; j < data.group.size(); j++)
                line->stats[i+j] = (int)tree->isGroup(data.group[j]);
            i += data.group.size()-1;
        }
        else if (statname[i].substr(0, 9)=="spr_leaf-") {
            const NodeSpr *nodespr = (line->trees->pruned_tree != NULL ?
                                      &(line->trees->pruned_spr) :
                                      &(line->trees->orig_spr));
            for (unsigned int j=0; j < data.spr_leaf.size(); j++) {
                line->stats[i+j] =
                    ( (nodespr->coal_node==NULL ||
                       !nodespr->coal_node->longname.compare(data.spr_leaf[j]))
                      ||
                      (nodespr->recomb_node==NULL ||
                       !nodespr->recomb_node->longname.compare(data.spr_leaf[j]))
                      );
            }
            i += data.spr_leaf.size() - 1;
        }
        else if (statname[i].substr(0, 5)=="coal-") {
            for (map<string,set<string> >::iterator it=
                     data.coalgroup_inds.begin();
                 it != data.coalgroup_inds.end(); ++it) {
                string hap1, hap2="";
                set<string> inds = it->second;
                set<string>::iterator it3 = inds.begin();
                hap1 = *it3;
                it3++;
                if (it3 != inds.end()) {
                    hap2 = *it3;
                    it3++;
                    assert(it3 == inds.end());
                }
                vector<double> tmpstats =
                    tree->coalGroup(hap1, hap2,
                                    data.coalgroups,
                                    data.coalgroup_names.size());
                for (unsigned int j=0; j < tmpstats.size(); j++)
                    line->stats[i+j] = tmpstats[j];
                i += tmpstats.size() - 1;
            }
        }
        else if (data.migstat.size() > 0) {
            for (unsigned int j=0; j < data.migstat.size(); j++) {
                if (statname[i] == data.migstat[j].name) {
                    line->stats[i] = (int)tree->haveMig(data.migstat[j].p,
                                                        data.migstat[j].t,
                                                        data.model,
                                                        data.migstat[j].hap);
                }
            }
        }
        else if (statname[i] == "cluster_stat") {
            line->stats[i] = tree->cluster_test(cluster_group, &line->stats[i+1]);
            //            line->stats[i] = (double)tree->num_prune_to_group(cluster_group);
            i++;
        }
        else {
            fprintf(stderr, "Error: unknown stat %s\n", statname[i].c_str());
            exit(1);
        }
    }
}


struct CompareBedLineSample
{
    bool operator()(const BedLine *l1, const BedLine *l2) const
    {
        return l1->sample < l2->sample;
    }
};


struct CompareBedLineEnd
{
    bool operator()(const BedLine *l1, const BedLine *l2) const
    {
        //could generalize and sort by start, but for this purpose we are
        //only comparing entries with same start
        assert(l1->start == l2->start);
        if (l1->end == l2->end)
            return (l1->sample < l2->sample);
        return l1->end < l2->end;
    }
};

void processNextBedLine(BedLine *line,
                        IntervalIterator<vector<double> > *results,
                        vector<string> &statname,
                        char *region_chrom, int region_start, int region_end,
                        ArgSummarizeData &data) {
    static int counter=0;
    static list<BedLine*> bedlist;

    if (line != NULL) {
        if (line->stats.size() == 0)
            scoreBedLine(line, statname, data);
        if (region_chrom != NULL) {
            assert(strcmp(region_chrom, line->chrom)==0);
            if (line->end > region_end) line->end = region_end;
            if (line->start < region_start) line->start = region_start;
            assert(line->start < line->end);
        }
    }
    if (!summarize) {
        // this little bit of code ensures that output is sorted. The
        // summarizeRegion code should ensure that it comes here sorted by
        // start coordinate, but not necessarily by end coordinate (at least
        // not when subsetting by individuals).
        // so, this stores BedLine elements just long enough until a new start
        // coordinate is encountered, then sorts them all by end coordinate (and
        // sample), then outputs them all.  Need to call this function one last
        // time with line==NULL to output the final set of rows.
        if (bedlist.size() > 0 &&
            (line==NULL || bedlist.front()->start < line->start)) {
            bedlist.sort(CompareBedLineEnd());
            for (list<BedLine*>::iterator it=bedlist.begin();
                 it != bedlist.end(); ++it) {
                BedLine *l = *it;
                if (html) printf("<tr><td>");
                printf("%s\t", l->chrom);
                if (html) printf("</td><td>");
                printf("%i\t", l->start);
                if (html) printf("</td><td>");
                printf("%i\t", l->end);
                if (html) printf("</td><td>");
                printf("%i", l->sample);
                for (unsigned int i=0; i < statname.size(); i++) {
                    if (statname[i]=="tree") {
                        if (!html) {
                            printf("\t%s", l->newick);
                        } else {
                            printf("</td><td nowrap>");
                            printf("\t<a href=\"http://mhubisz.genome-mirror.cshl.edu/cgi-bin/phyloGif?phyloGif_width=240&phyloGif_height=512&phyloGif_branchLengths=on&phyloGif_underscores=on&phyloGif_tree=");
                           for (unsigned int i=0; i < strlen(l->newick); i++) {
                              if (l->newick[i]=='(') {
                                printf("%%28");
                              } else if (l->newick[i]==':') {
                                printf("%%3A");
                              } else if (l->newick[i]==',') {
                                printf("%%2C");
                              } else if (l->newick[i]==')') {
                                printf("%%29");
                              } else if (l->newick[i]==';') {
                                printf("%%3B");
                              } else if (l->newick[i]=='&') {
                                printf("%%26");
                              } else if (l->newick[i]=='[') {
                                printf("%%5B");
                              } else if (l->newick[i]==']') {
                                printf("%%5D");
                              } else printf("%c", l->newick[i]);
                           }
                           printf("%%0D%%0A\">%s</a>", l->newick);
                        }
                    } else {
                        if (html) printf("</td><td>");
                        printf("\t%g", l->stats[i]);
                    }
                }
                printf("\n");
                if (html) printf("</td></tr>\n");
                delete l;
            }
            bedlist.clear();
        }
        if (line != NULL) bedlist.push_back(line);
    } else {
        if (line != NULL) {
            results->append(line->chrom, line->start, line->end, line->stats);
            counter++;
            if (counter%100==0) {
                checkResults(results);
            }
            delete line;
        }
    }
}

class SnpStream {
public:
    SnpStream(TabixStream *snp_in) : snp_in(snp_in) {
        char tmp[1000], c;
        string str;
        // first string should be NAMES or #NAMES
        assert(1==fscanf(snp_in->stream, "%s", tmp));
        done=0;
        while ('\n' != (c=fgetc(snp_in->stream)) && c!=EOF) {
            assert(c=='\t');
            assert(1==fscanf(snp_in->stream, "%s", tmp));
            str = string(tmp);
            inds.push_back(str);
        }
        firstLine=true;
    }

    int readNext() {
        int tmpStart;
        char a;
	vector<set<string> > allele_inds;
	int count[4]={0,0,0,0};
        bool alreadyRead=false;
        if (done) return 1;
        if (firstLine) {
            if (EOF == fscanf(snp_in->stream, "%s", chr)) {
                done=1;
                return 1;
            }
            isBed=!(strcmp(chr, "REGION")==0 || strcmp(chr, "#REGION")==0);
            firstLine=false;
            if (!isBed) {
                char c;
                while ('\n'!=(c=fgetc(snp_in->stream)) && c!=EOF);
                if (c==EOF) {
                    done=1;
                    return 1;
                } else {
                    return readNext();
               }
            } else {
                assert(2==fscanf(snp_in->stream, "%i %i", &tmpStart, &coord));
                alreadyRead=true;
            }
        }
        if (isBed) {
            if (!alreadyRead) {
                if (EOF==fscanf(snp_in->stream, "%s %i %i", chr, &tmpStart, &coord)) {
                    done=1;
                    return 1;
                }
            }
            assert(tmpStart==coord-1);
        } else {
            if (EOF == fscanf(snp_in->stream, "%i", &coord)) {
                done = 1;
                return 1;
            }
        }
	for (int i=0; i < 4; i++)
	    allele_inds.push_back(set<string>());
        assert('\t' == fgetc(snp_in->stream));
	allele1=allele2='N';
        for (unsigned int i=0; i < inds.size(); i++) {
            a=fgetc(snp_in->stream);
            if (a=='N') continue;
            a = toupper(a);
            assert(a=='A' || a=='C' || a=='G' || a=='T');
	    int aval = dna2int[(int)a];
	    allele_inds[aval].insert(inds[i]);
	    count[aval]++;
        }
        //make sure that allele1 is always minor allele
	// if there are more than 2 alleles, then allele2 is
	// most common and allele1 is second-most-common.
	// Others are treated as Ns
	int allele2_val = 0;
	for (int i=1; i < 4; i++)
	    if (count[i] > count[allele2_val])
		allele2_val = i;
	int allele1_val = (allele2_val == 0 ? 1 : 0);
	for (int i=1; i < 4; i++) {
	    if (i == allele2_val) continue;
	    if (count[i] > count[allele1_val])
		allele1_val = i;
	}
        // go to next line if this line is invariant
        if (count[allele1_val] == 0) return readNext();
        moreThanTwoAlleles=false;
        for (int i=0; i < 4; i++) {
            if (i == allele1_val || i == allele2_val) continue;
            if (count[i] > 0) moreThanTwoAlleles=true;
        }
	allele1 = int2dna[allele1_val];
	allele2 = int2dna[allele2_val];
	allele1_inds = allele_inds[allele1_val];
	allele2_inds = allele_inds[allele2_val];
        return 0;
    }


    void scoreAlleleAge(BedLine *l, vector<string> statname,
                        ArgSummarizeData &data) {
        int num_derived, total;
        assert(l->start < coord);
        assert(l->end >= coord);
        Tree *t;
        if (l->trees->pruned_tree != NULL)
            t = l->trees->pruned_tree;
        else t = l->trees->orig_tree;

        set<string> prune;
        set<string> derived_in_tree;
        for (map<string,int>::iterator it=t->nodename_map.begin();
             it != t->nodename_map.end(); ++it) {
            if (t->nodes[it->second]->nchildren != 0) continue;
            if (allele1_inds.find(it->first) != allele1_inds.end())
                derived_in_tree.insert(it->first);
            else if (allele2_inds.find(it->first) == allele2_inds.end())
                prune.insert(it->first);
        }
        Tree *t2=NULL;
        if (prune.size() > 0) {
            t2 = t->copy();
            t2->prune(prune);
            t = t2;
        }
        set<Node*> derived;
        for (set<string>::iterator it=derived_in_tree.begin();
             it != derived_in_tree.end(); ++it) {
            map<string,int>::iterator it2 = t->nodename_map.find(*it);
            assert(it2 != t->nodename_map.end());
            derived.insert(t->nodes[it2->second]);
        }
        num_derived = (int)derived.size();
        total = (t->nnodes+1)/2;

        set<Node*>lca = t->lca(derived);
        set<Node*>lca2;
        int major_is_derived=0;
        if (lca.size() > 1) {
            set<Node*> derived2;
            for (map<string,int>::iterator it=t->nodename_map.begin();
                 it != t->nodename_map.end(); ++it) {
                if (t->nodes[it->second]->nchildren != 0) continue;
                if (derived.find(t->nodes[it->second]) == derived.end()) {
                    derived2.insert(t->nodes[it->second]);
                }
            }
            set <Node*>lca2 = t->lca(derived2);
            if (lca2.size() < lca.size()) {
                major_is_derived=1;
                lca = lca2;
            }
        }
        double age=0.0;
        double minage=0.0;
        if (!moreThanTwoAlleles) {
            for (set<Node*>::iterator it4=lca.begin(); it4 != lca.end(); ++it4) {
                Node *n = *it4;
                assert(n != t->root);
                double tempage = n->age + (n->parent->age - n->age)/2;  //midpoint
                minage = n->age;
                if (tempage > age) age = tempage;
            }
        }
        if (moreThanTwoAlleles || num_derived == 0 || total-num_derived == 0)
            age = -1;
        scoreBedLine(l, statname, data, age, minage, lca.size()==1);
        l->derAllele = (major_is_derived ? allele2 : allele1);
        l->otherAllele = (major_is_derived ? allele1 : allele2);
        l->derFreq = (major_is_derived ? total-num_derived : num_derived);
        l->otherFreq = (major_is_derived ? num_derived : total - num_derived);
        l->infSites = (lca.size() == 1 && !moreThanTwoAlleles);

        if (t2 != NULL) delete t2;
    }

    TabixStream *snp_in;
    vector<string> inds;
    set<string> allele1_inds;
    set<string> allele2_inds;
    char allele1, allele2;  //minor allele, major allele
    char chr[100];
    int coord;  //1-based
    int done;
    bool isBed;
    bool firstLine;
    bool moreThanTwoAlleles;
};



void print_summaries(vector<double> &stat) {
    double meanval=0;
    int have_mean=0;
    for (int j=1; j <= summarize; j++) {
        if (getMean==j) {
            if (html) printf("</td><td>");
            if (stat.size() > 0) {
                meanval = compute_mean(stat);
                have_mean=1;
                printf("\t%g", meanval);
            } else printf("\tNA");
        } else if (getStdev==j) {
            if (html) printf("</td><td>");
            if (stat.size() > 1) {
                if (!have_mean)
                    meanval = compute_mean(stat);
                printf("\t%g", compute_stdev(stat, meanval));
            } else printf("\tNA");
        } else if (getQuantiles==j) {
            if (html) printf("</td><td>");
            if (stat.size() > 0) {
                vector<double> q = compute_quantiles(stat, quantiles);
                for (unsigned int k=0; k < quantiles.size(); k++) {
                    printf("\t%g", q[k]);
                }
            } else {
                for (unsigned int k=0; k < quantiles.size(); k++) {
                    if (html) printf("</td><td>");
                    printf("\tNA");
                }
            }
        }
    }
}


int summarizeRegionBySnp(Config *config, const char *region,
                         set<string> inds, vector<string> statname,
                         ArgSummarizeData &data) {
    TabixStream snp_infile(config->snpfile, region, config->tabix_dir);
    TabixStream infile(config->argfile, region, config->tabix_dir);
    vector<string> token;
    map<int,BedLine*> last_entry;
    map<int,BedLine*>::iterator it;
    char chrom[1000], c;
    int start, end, sample;
    BedLine *l=NULL;
    const ArgModel *model = data.model;

    if (snp_infile.stream == NULL) return 1;
    if (infile.stream == NULL) return 1;
    while (EOF != (c=fgetc(infile.stream))) {
        ungetc(c, infile.stream);
        if (c != '#') break;
        while ('\n' != (c=fgetc(infile.stream))) {
            if (c==EOF) return 0;
        }
    }
    SnpStream snpStream = SnpStream(&snp_infile);
    char *newick;
    while (1) {
        if (EOF==fscanf(infile.stream, "%s %i %i %i",
                        chrom, &start, &end, &sample)) return 0;
        assert('\t' == fgetc(infile.stream));
        newick = fgetline(infile.stream);
        if (sample >= config->burnin && (config->sample_num == 0 ||
                                         config->sample_num == sample))
            break;
        delete [] newick;
    }
    chomp(newick);

    while (1) {
        list<BedLine*> bedlist;
        bedlist.clear();
        snpStream.readNext();
        if (snpStream.done) break;
        // first check already-parsed BedLines and score any that overlap SNP
        for (it=last_entry.begin(); it != last_entry.end(); it++) {
            l = it->second;
            if (l->start < snpStream.coord && l->end >= snpStream.coord) {
                snpStream.scoreAlleleAge(l, statname, data);
                bedlist.push_back(l);
            }
        }
        //now look through bed file until we get to one that starts after SNP
        while (start != -1 && snpStream.coord > start) {
            it = last_entry.find(sample);
            if (it == last_entry.end() ||
                it->second->trees->orig_spr.recomb_node == NULL) {
                SprPruned *trees;
                if (it != last_entry.end()) {
                    l = it->second;
                    delete l->trees;
                    delete &*l;
                }
                trees = new SprPruned(newick, inds, model);
                l = new BedLine(chrom, start, end, sample, newick, trees);
                last_entry[sample] = l;
            } else {
                l = it->second;
                l->trees->update(newick, model);
                free(l->newick);
                l->newick = (char*)malloc((strlen(newick)+1)*sizeof(char));
                strcpy(l->newick, newick);
                l->start = start;
                l->end = end;
            }
            if (snpStream.coord <= end) {
                snpStream.scoreAlleleAge(l, statname, data);
                bedlist.push_back(l);
            }
            while (1) {
                if (4 != fscanf(infile.stream, "%s %i %i %i",
                                chrom, &start, &end, &sample)) {
                    start = -1;
                    break;
                } else {
                    assert('\t' == fgetc(infile.stream));
                    delete [] newick;
                    newick = fgetline(infile.stream);
                    chomp(newick);
                    if (sample >= config->burnin && (config->sample_num == 0 || config->sample_num==sample))
                       break;
                }
            }
        }
        if (bedlist.size() > 0) {
            if (summarize == 0) {
                bedlist.sort(CompareBedLineSample());
                for (list<BedLine*>::iterator it=bedlist.begin();
                     it != bedlist.end(); ++it) {
                    BedLine *l = *it;
                    if (html) printf("<tr><td>");
                    printf("%s\t", l->chrom);
                    if (html) printf("</td><td>");
                    printf("%i\t", snpStream.coord-1);
                    if (html) printf("</td><td>");
                    printf("%i\t", snpStream.coord);
                    if (html) printf("</td><td>");
                    printf("%i\t", l->sample);
                    if (html) printf("</td><td>");
                    printf("%c\t", l->derAllele);
                    if (html) printf("</td><td>");
                    printf("%c\t", l->otherAllele);
                    if (html) printf("</td><td>");
                    printf("%i\t", l->derFreq);
                    if (html) printf("</td><td>");
                    printf("%i", l->otherFreq);
                    for (unsigned int i=0; i < statname.size(); i++) {
                        if (statname[i]=="tree") {
                            if (html) printf("</td><td>");
                            printf("\t%s", l->newick);
                        } else if (statname[i]=="infSites") {
                            if (html) printf("</td><td>");
                            printf("\t%i", (int)(l->stats[i]==1));
                        } else {
                            if (html) printf("</td><td>");
                            printf("\t%g", l->stats[i]);
                        }
                    }
                    if (html) printf("</td></tr>");
                    printf("\n");
                    l->stats.clear();
                }
            } else {
                //now output three versions- one for all samples,
                //one for same derived allele, one for infinite sites
                BedLine* first = *(bedlist.begin());
                int same=0, diff=0, infsites=0,
                    derFreq, otherFreq;
                char derAllele, otherAllele;
                for (list<BedLine*>::iterator it=bedlist.begin();
                     it != bedlist.end(); ++it) {
                    BedLine *l = *it;
                    if (l->derAllele == first->derAllele) same++; else diff++;
                    infsites += l->infSites;
                }
                if (same >= diff) {
                    derAllele = first->derAllele;
                    otherAllele = first->otherAllele;
                    derFreq = first->derFreq;
                    otherFreq = first->otherFreq;
                } else {
                    derAllele=first->otherAllele;
                    otherAllele = first->derAllele;
                    derFreq = first->otherFreq;
                    otherFreq = first->derFreq;
                }
                if (html) printf("<tr><td>");
                printf("%s\t", l->chrom);
                if (html) printf("</td><td>");
                printf("%i\t", snpStream.coord-1);
                if (html) printf("</td><td>");
                printf("%i\t", snpStream.coord);
                if (html) printf("</td><td>");
                printf("%c\t", derAllele);
                if (html) printf("</td><td>");
                printf("%c\t", otherAllele);
                if (html) printf("</td><td>");
                printf("%i\t", derFreq);
                if (html) printf("</td><td>");
                printf("%i\t", otherFreq);
                if (html) printf("</td><td>");
                printf("%i\t", (int)bedlist.size());
                if (html) printf("</td><td>");
                printf("%i", infsites);
                if (html) printf("</td><td>");
                for (unsigned int i=0; i < statname.size(); i++) {
                    if (statname[i] != "inf_sites") {
                        // first compute stats across all
                        vector<double> stat;
                        for (list<BedLine*>::iterator it=bedlist.begin();
                             it != bedlist.end(); ++it) {
                            BedLine *l = *it;
                            stat.push_back(l->stats[i]);
                        }
                        print_summaries(stat);

                        stat.clear();
                        //now stats for infinite sites set
                        for (list<BedLine*>::iterator it=bedlist.begin();
                             it != bedlist.end(); ++it) {
                            BedLine *l = *it;
                            if (l->infSites)
                                stat.push_back(l->stats[i]);
                            l->stats.clear();
                        }
                        print_summaries(stat);
                    }
                }
                printf("\n");
                if (html) printf("</td></tr>\n");
            }
        }
    }
    delete [] newick;

    for (map<int,BedLine*>::iterator it=last_entry.begin();
         it != last_entry.end(); ++it) {
        BedLine *l = it->second;
        delete(l->trees);
        delete l;
    }
    return 0;
}


int summarizeRegionNoSnp(Config *config, const char *region,
                         set<string> inds, vector<string>statname,
                         ArgSummarizeData &data) {
    TabixStream *infile;
    char c;
    char *region_chrom = NULL;
    char chrom[1000];
    vector<string> token;
    int region_start=-1, region_end=-1, start, end, sample;
    IntervalIterator<vector<double> > results;
    queue<BedLine*> bedlineQueue;
    map<int,BedLine*> bedlineMap;
    map<int,SprPruned*> trees;
    map<int,SprPruned*>::iterator it;
    const ArgModel *model = data.model;
    /*
      Class BedLine contains chr,start,end, newick tree, parsed tree.
      parsed tree may be NULL if not parsing trees but otherwise will
      be populated, either by parsing the newick or an SPR operation
      on previous tree.
      Parsed tree has recomb_node, recomb_time, coal_node, coal_time set
      (recomb_node==NULL => no recomb. Only happens in full tree at end
      of regions analyzed by arg-sample)

      Queue bedlineQueue contains pointers to this class, will be output to
      results in order (first in, first out).
      bedlineMap<int,bedlineQueue> maps samples to pointers of the most recently
      read instance of bedline for each sample. It points to the same objects
      as bedlineQueue (not copies).

      For each line of file:
      Read chr, start, end, sample, tree string.
      Look up bedlineMap<sample> = lastSample
      If (lastSample == NULL) {
         parse tree. Make new bedline object, add it to bedlineMap and end
        of bedlineQueue.
      } else if (lastSample->recomb_node != NULL) {
        apply SPR to lastSample->tree to create new parsed tree. Use this tree
        to create new bedline object, add it to bedlineMap and end of
        bedlineQueue.
      } else { //lastSample does not end in recomb
        assert that lastSample->end==start and lastSample->chr==chr
        update lastSample->end=end
        determine if tree string ends in recomb, update recomb_node,
        recomb_time, etc if so.
        (This is a tricky part esp if there is pruning involved).
      }
      while (first element of queue ends in recomb) {
        compute statistics for first element of queue
        add to intervaliterator
        if bedlineMap<sample>==(first element in queue),
          set bedlineMap<sample>=NULL
        pop from queue
      }

      After reading all lines:
      go through queue and dump everything to intervalIterator...

    */

    infile = new TabixStream(config->argfile, region, config->tabix_dir);
    if (infile->stream == NULL) return 1;

    //parse region to get region_chrom, region_start, region_end.
    // these are only needed to truncate results which fall outside
    // of the boundaries (tabix returns anything that overlaps)
    if (region != NULL) {
        split(region, "[:-]", token);
        if (token.size() != 3) {
            fprintf(stderr,
                    "Error: bad region format (%s); should be chr:start-end\n",
                    region);
            return 1;
        }
        region_chrom = new char[token[0].size()+1];
        //remove commas from integer coordinates in case they are
        // copied from browser
        token[1].erase(std::remove(token[1].begin(), token[1].end(), ','),
                       token[1].end());
        token[2].erase(std::remove(token[2].begin(), token[2].end(), ','),
                       token[2].end());
        strcpy(region_chrom, token[0].c_str());
        region_start = atoi(token[1].c_str())-1;
        region_end = atoi(token[2].c_str());
    }
    while (EOF != (c=fgetc(infile->stream))) {
        ungetc(c, infile->stream);
        if (c!='#') break;
        while ('\n' != (c=fgetc(infile->stream))) {
            if (c==EOF) return 0;
        }
    }

    while (EOF != fscanf(infile->stream, "%s %i %i %i",
                         chrom, &start, &end, &sample)) {
        assert('\t'==fgetc(infile->stream));
        char* newick = fgetline(infile->stream);
        if (config->sample_num != 0 && sample != config->sample_num) {
            delete [] newick;
            continue;
        }
        chomp(newick);
        if (sample < config->burnin) {
            delete [] newick;
            continue;
        }
        it = trees.find(sample);
        if (it == trees.end())   //first tree from this sample
            trees[sample] = new SprPruned(newick, inds, model);
        else trees[sample]->update(newick, model);

        map<int,BedLine*>::iterator it3 = bedlineMap.find(sample);
        BedLine *currline;
        if (it3 == bedlineMap.end()) {
            currline = new BedLine(chrom, start, end, sample, newick,
                                   trees[sample]);
            bedlineMap[sample] = currline;
            bedlineQueue.push(currline);
        } else {
            currline = it3->second;
            assert(strcmp(currline->chrom, chrom)==0);
            assert(currline->end == start);
            currline->end = end;
        }

        //assume orig_spr.recomb_node == NULL is a rare occurrence that happens
        // at the boundaries of regions analyzed by arg-sample; treat these as
        // recombination events
        if (trees[sample]->orig_spr.recomb_node == NULL ||
            trees[sample]->pruned_tree == NULL ||
            trees[sample]->pruned_spr.recomb_node != NULL) {
            scoreBedLine(currline, statname, data);
            bedlineMap.erase(sample);
        }

        while (bedlineQueue.size() > 0) {
            BedLine *firstline = bedlineQueue.front();
            if (firstline->stats.size() == statname.size()) {
                processNextBedLine(firstline, &results, statname,
                                   region_chrom, region_start, region_end,
                                   data);
                bedlineQueue.pop();
            } else break;
        }
        delete [] newick;
    }
    infile->close();
    delete infile;

    while (bedlineQueue.size() > 0) {
        BedLine *firstline = bedlineQueue.front();
        processNextBedLine(firstline, &results, statname,
                           region_chrom, region_start, region_end, data);

        //        delete firstline;
        bedlineQueue.pop();
    }

    if (summarize) {
        results.finish();
        checkResults(&results);
    } else {
        processNextBedLine(NULL, &results, statname, region_chrom,
                           region_start, region_end, data);
    }

    it = trees.begin();
    while (it != trees.end()) {
        delete it->second;
        advance(it, 1);
    }
    if (region_chrom != NULL) delete[] region_chrom;
    return 0;
}

int summarizeRegion(Config *config, const char *region,
                    set<string> inds, vector<string>statname,
                    ArgSummarizeData &data) {
    if (config->snpfile.empty())
        return summarizeRegionNoSnp(config, region, inds, statname, data);
    else
        return summarizeRegionBySnp(config, region, inds, statname, data);
}


int main(int argc, char *argv[]) {
    Config c;
    int ret = c.parse_args(argc, argv);
    ArgSummarizeData data;
    if (ret)
        return ret;

    quiet = c.quiet;
    if (c.argfile.empty()) {
        fprintf(stderr, "Error: must specify argfile\n");
        return 1;
    }
    if (!c.logfile.empty()) {
        data.model = new ArgModel(c.logfile.c_str());
    } else data.model = NULL;
    if (c.html) {
        html=true;
        printf("<html>\n");
        //printf("<link rel=\"stylesheet\" type=\"text/css\" href=\"ARGweaver.css\" />\n");
        if (c.rawtrees) {
            printf("<frameset cols=\"66%%,*\">\n");
        }
    }
    set<string> haps;
    if (!c.hapfile.empty()) {
        ifstream in(c.hapfile.c_str());
        string line;
        if (in.is_open()) {
            while ( getline(in, line) ) {
                haps.insert(line);
            }
            in.close();
        } else {
            fprintf(stderr, "Error opening %s.\n", c.hapfile.c_str());
            return 1;
        }
    }
    set<string> inds;
    if (!c.indfile.empty()) {
        ifstream in(c.indfile.c_str());
        string line;
        if (in.is_open()) {
            while ( getline(in, line) ) {
                inds.insert(line);
                haps.insert(line + string("_1"));
                haps.insert(line + string("_2"));
            }
            in.close();
        } else {
            fprintf(stderr, "Error opening %s.\n", c.indfile.c_str());
            return 1;
        }
    }

    vector<string> statname;
    if (c.tmrca)
        statname.push_back(string("tmrca"));
    if (c.pi)
        statname.push_back(string("pi"));
    if (c.branchlen)
        statname.push_back(string("branchlen"));
    if (c.recomb)
        statname.push_back(string("recomb"));
    if (c.breaks)
        statname.push_back(string("breaks"));
    if (c.tmrca_half)
        statname.push_back(string("tmrca_half"));
    if (c.rth)
        statname.push_back(string("rth"));
    if (c.popsize)
        statname.push_back(string("popsize"));
    if (c.max_coal_rate) {
        if (data.model == NULL) {
            fprintf(stderr, "--log-file required for --max-coalrate\n");
            return 1;
        }
        statname.push_back(string("max_coal_rate"));
    }
    if (c.allele_age) {
        statname.push_back(string("allele_age"));
        statname.push_back(string("inf_sites"));
    }
    if (c.min_allele_age)
        statname.push_back(string("min_allele_age"));
    if (c.zero)
        statname.push_back(string("zero_len"));

    vector<string> groupfiles;
    if (!c.groupfile.empty()) {
        split(c.groupfile.c_str(), ',', groupfiles);
        for (unsigned int i=0; i < groupfiles.size(); i++) {
            char tmp[100];
            sprintf(tmp, "%i", i);
            statname.push_back(string("group") + string(tmp));
        }
    }
    if (!(c.coalgroup_file.empty() && c.coalgroup_inds_file.empty())) {
        if (c.coalgroup_file.empty() || c.coalgroup_inds_file.empty()) {
            fprintf(stderr, "Error: --coal-groups and --coal-group-inds"
                    " should be used together\n");
            return 1;
        }
        FILE *infile = fopen(c.coalgroup_inds_file.c_str(), "r");
        char tmp[1000], tempc;
        if (infile == NULL) {
            fprintf(stderr, "Error opening %s.\n",
                    c.coalgroup_inds_file.c_str());
            return 1;
        }
        while (EOF != fscanf(infile, "%s", tmp)) {
            string ind = string(tmp);
            //read in first haploid chrom for ind
            assert(1 == fscanf(infile, "%s", tmp));
            set<string> indset;
            indset.insert(string(tmp));

            // check to see if there's a second haploid chrom for ind
            while ('\n' != (tempc=fgetc(infile)) && tempc!=EOF)
                if (!isspace(tempc)) break;
            if (!isspace(tempc)) {
                ungetc(tempc, infile);
                assert(1==fscanf(infile, "%s", tmp));
                indset.insert(string(tmp));
            }
            //make sure there are no more
            while ('\n' != (tempc=fgetc(infile)) && tempc!=EOF) {
                if (!isspace(tempc)) {
                    fprintf(stderr, "Error: expected only two haploids per"
                            " line in group_indfile\n");
                    return 1;
                }
            }
            data.coalgroup_inds[ind] = indset;
            if (tempc==EOF) break;
        }
        fclose(infile);
        infile = fopen(c.coalgroup_file.c_str(), "r");
        if (infile == NULL) {
            fprintf(stderr, "Error opening %s.\n",
                    c.coalgroup_file.c_str());
            return 1;
        }
        char group[1000], hap[1000];
        while (EOF != fscanf(infile, "%s %s", hap, group)) {
            int group_num;
            for (group_num=0;
                 group_num < (int)data.coalgroup_names.size();
                 group_num++)
                if (strcmp(data.coalgroup_names[group_num].c_str(),
                           group)==0) break;
            if (group_num == (int)data.coalgroup_names.size())
                data.coalgroup_names.push_back(string(group));
            data.coalgroups[string(hap)] = group_num;
        }
        fclose(infile);

        for (map<string,set<string> >::iterator it=data.coalgroup_inds.begin();
             it != data.coalgroup_inds.end(); ++it) {
            for (unsigned int j=0; j < data.coalgroup_names.size(); j++) {
                string tmpname = "coal-" + it->first + "-" +
                    data.coalgroup_names[j];
                statname.push_back(tmpname);
            }
            statname.push_back("coal-" + it->first + "-mixedgroup");
        }
    }
    if (c.coalcounts) {
        if (data.model == NULL) {
            fprintf(stderr, "Error: --log-file required with --coalcounts\n");
            return 1;
        }
        for (int i=0; i < data.model->ntimes; i++) {
            char tmp[1000];
            sprintf(tmp, "coalcounts.%i", i);
            statname.push_back(string(tmp));
        }
    }
    if (c.coalcounts_cluster) {
        if (data.model == NULL) {
            fprintf(stderr, "Error: --log-file required with --coalcounts-cluster\n");
            return 1;
        }
        for (int i=0; i < data.model->ntimes; i++) {
            char tmp[1000];
            sprintf(tmp, "coalcounts-cluster.%i", i);
            statname.push_back(string(tmp));
        }
    }
    if (c.recombs_per_time) {
        if (data.model == NULL) {
            fprintf(stderr, "Error: --log-file required with --recombs-per-time\n");
            return 1;
        }
        for (int i=0; i < data.model->ntimes; i++) {
            char tmp[1000];
            sprintf(tmp, "recombs.%i", i);
            statname.push_back(string(tmp));
        }
    }
    if (c.invis_recombs_per_time) {
        if (data.model == NULL) {
            fprintf(stderr, "Error: --log-file required with --invis-recombs-per-time\n");
            return 1;
        }
        for (int i=0; i < data.model->ntimes; i++) {
            char tmp[1000];
            sprintf(tmp, "invis-recombs.%i", i);
            statname.push_back(string(tmp));
        }
    }
    if (c.branchlen_per_time) {
        if (data.model == NULL) {
            fprintf(stderr, "Error: --log-file required with --branchlen-per-time\n");
            return 1;
        }
        for (int i=0; i < data.model->ntimes; i++) {
            char tmp[1000];
            sprintf(tmp, "branchlen.%i", i);
            statname.push_back(string(tmp));
        }
    }
    if (!c.node_dist.empty()) {
        vector<string> tokens1, tokens2;
        split(c.node_dist.c_str(), ';', tokens1);
        for (unsigned int i=0; i < tokens1.size(); i++) {
            split(tokens1[i].c_str(), ',', tokens2);
	    if (tokens2.size() == 1) {
		statname.push_back(string("node_dist-") + tokens2[0]);
		node_dist_leaf1.push_back(tokens2[0]);
		node_dist_leaf2.push_back("");
	    } else {
		if (tokens2.size() != 2) {
		    fprintf(stderr, "Bad format to --node-dist argument; expect"
			    " one or two leaf names separated by comma");
		    return 1;
		}
		statname.push_back(string("node_dist-") + tokens2[0]
				   + string(",") + tokens2[1]);
		node_dist_leaf1.push_back(tokens2[0]);
		node_dist_leaf2.push_back(tokens2[1]);
	    }
        }
    }
    if (c.node_dist_all) {
        if (haps.size() == 0) {
            fprintf(stderr, "Must use --subset or --subset-inds argument with"
                    " --node-dist-all\n");
            return 1;
        }
        for (set<string>::iterator it=haps.begin(); it != haps.end(); it++) {
            for (set<string>::iterator it2 = it; it2 != haps.end(); it2++) {
                if (it != it2) {
                    statname.push_back(string("node_dist-") + *it + string(",") + *it2);
                    node_dist_leaf1.push_back(*it);
                    node_dist_leaf2.push_back(*it2);
                }
            }
        }
    }
    if (!c.min_coal_time.empty()) {
        vector<string> tokens1, tokens2;
        split(c.min_coal_time.c_str(), ';', tokens1);
        for (unsigned int i=0; i < tokens1.size(); i++) {
            split(tokens1[i].c_str(), ',', tokens2);
            if (tokens2.size() != 2) {
                fprintf(stderr, "--min-coal-time; bad argument\n");
                return 1;
            }
            statname.push_back(string("min_coal_time-") + tokens1[i]);
            min_coal_time_ind1.push_back(tokens2[0]);
            min_coal_time_ind2.push_back(tokens2[1]);
        }
    }
    if (c.min_coal_time_all) {
        if (inds.size() == 0) {
            fprintf(stderr, "Must use --subset-inds argument with"
                    " --min-coal-time-all\n");
            return 1;
        }
        for (set<string>::iterator it=inds.begin(); it != inds.end(); it++) {
            for (set<string>::iterator it2=it; it2 != inds.end(); it2++) {
                if (it != it2) {
                    statname.push_back(string("min_coal_time-") + *it + string(",") + *it2);
                    min_coal_time_ind1.push_back(*it);
                    min_coal_time_ind2.push_back(*it2);
                }
            }
        }
    }
    if (!c.spr_leaf_names.empty()) {
        split(c.spr_leaf_names.c_str(), ',', data.spr_leaf);
        for (unsigned int i=0; i < data.spr_leaf.size(); i++) {
            statname.push_back(string("spr_leaf-") + data.spr_leaf[i]);
        }
    }
    if (!c.cluster_group_file.empty()) {
        FILE *infile = fopen(c.cluster_group_file.c_str(), "r");
        char leaf_name[10000];
        if (infile == NULL) {
            fprintf(stderr, "Error opening %s\n", c.cluster_group_file.c_str());
            return 0;
        }
        while (EOF != fscanf(infile, "%s", leaf_name)) {
            cluster_group.insert(string(leaf_name));
        }
        fclose(infile);
        statname.push_back(string("cluster_stat"));
        statname.push_back(string("cluster_time"));
    }

    if (!c.ind_dist.empty()) {
	vector<string> tokens1, tokens2;
	split(c.ind_dist.c_str(), ';', tokens1);
	for (unsigned int i=0; i < tokens1.size(); i++) {
	    split(tokens1[i].c_str(), ',', tokens2);
	    if (tokens2.size() != 2) {
		fprintf(stderr, "Bad format to --ind-dist argument; expect"
			" two leaf names separated by comma");
		return 1;
	    }
	    statname.push_back(string("ind_dist-") + tokens2[0]
			       + string(",") + tokens2[1] +
			       string("-min"));
	    statname.push_back(string("ind_dist-") + tokens2[0]
			       + string(",") + tokens2[1] +
			       string("-max"));
	    statname.push_back(string("ind_dist-") + tokens2[0]
			       + string(",") + tokens2[1] +
			       string("-homo"));
	    ind_dist_leaf1.push_back(tokens2[0]);
	    ind_dist_leaf2.push_back(tokens2[1]);
	}
    }
    if (!c.migfile.empty()) {
        if (data.model == NULL) {
            fprintf(stderr, "--log-file required with --mig-file\n");
            exit(1);
        }
        FILE *infile = fopen(c.migfile.c_str(), "r");
        if (infile == NULL) {
            fprintf(stderr, "Error opening %s\n", c.migfile.c_str());
            exit(1);
        }
        char migname[1000];
        int p[2];
        double dt;
        while (EOF != fscanf(infile, "%s %i %i %lf", migname, &p[0], &p[1],
                             &dt)) {
            data.migstat.push_back(MigStat(string(migname), p, dt, data.model));
            statname.push_back(string(migname));
        }
        fclose(infile);
    }
    if (!c.hapmigfile.empty()) {
        if (data.model == NULL) {
            fprintf(stderr, "--log-file required with --hap-mig-file\n");
            exit(1);
        }
        FILE *infile = fopen(c.hapmigfile.c_str(), "r");
        if (infile == NULL) {
            fprintf(stderr, "Error opening %s\n", c.migfile.c_str());
            exit(1);
        }
        char migname[1000], hap[1000];
        int p[2];
        double dt;
        while (EOF != fscanf(infile, "%s %i %i %lf %s", migname, &p[0], &p[1],
                             &dt, hap)) {
            data.migstat.push_back(MigStat(string(migname), p, dt, data.model, string(hap)));
            statname.push_back(string(migname));
        }
        fclose(infile);

    }
    if (c.rawtrees)
        statname.push_back(string("tree"));

    if (c.numsample)
        getNumSample=++summarize;
    if (c.mean)
        getMean=++summarize;
    if (c.stdev)
        getStdev=++summarize;;
    if (!c.quantile.empty()) {
        getQuantiles=++summarize;
        vector<string> tokens;
        split(c.quantile.c_str(), ',', tokens);
        for (unsigned int i=0; i < tokens.size(); i++) {
            double q=atof(tokens[i].c_str());
            //              fprintf(stderr, "getting quantile %lf\n",q);
            quantiles.push_back(q);
        }
    }

    if ((!c.region.empty()) && (!c.bedfile.empty())) {
        fprintf(stderr, "Error: --bed and --region cannot be used together.\n");
        return 1;
    }
    if (statname.size() == 0) {
        fprintf(stderr, "Error: need to specify a tree statistic\n");
        return 1;
    }

    if (summarize && statname.size()==0) {
        fprintf(stderr,
                "Error: need to specify a tree statistic (e.g., --tmrca,"
                " --pi, --popsize, --recomb, --allele-age, etc)\n");
        return 1;
    }
    if (summarize && c.rawtrees) {
        fprintf(stderr, "Error: --trees not compatible with summary statistics"
                " (--mean, --quantile, --stdev, --numsample)\n");
        return 1;
    }
    if ((c.recomb || c.breaks) && !c.snpfile.empty()) {
        fprintf(stderr, "Error: cannot use --recomb or --breaks with"
                " --allele-age\n");
        return 1;
    }
    if ((c.allele_age || c.min_allele_age) && c.snpfile.empty()) {
        fprintf(stderr, "Error: need to specify snp file with --snp to use"
                " --allele-age or --min-allele-age\n");
        return 1;
    }

    if (!c.noheader) {
        printf("## %s\n", VERSION_INFO);
        if (html) printf("<br>\n");
        printf("##");
        for (int i=0; i < argc; i++) printf(" %s", argv[i]);
        printf("\n");
        if (html) printf("<br><table><tr><td>\n");
        printf("##chrom\t");
        if (html) printf("</td><td>\n");
        printf("chromStart\t");
        if (html) printf("</td><td>\n");
        printf("chromEnd");
        if (summarize==0) {
            if (html) printf("</td><td>\n");
            printf("\tMCMC_sample");
        }
        if (!c.snpfile.empty()) {
            if (html) printf("</td><td>\n");
            printf("\tderAllele");
            if (html) printf("</td><td>\n");
            printf("\tancAllele");
            if (html) printf("</td><td>\n");
            printf("\tderFreq");
            if (html) printf("</td><td>\n");
            printf("\tancFreq");
        }
        if (c.snpfile.empty() && getNumSample > 0) {
            if (html) printf("</td><td>\n");
            printf("\tnumsample");
        }
        if (summarize && !c.snpfile.empty()) {
            if (html) printf("</td><td>\n");
            printf("\tnumsample-all");
            if (html) printf("</td><td>\n");
            printf("\tnumsample-infsites");
        }
        vector<string> stattype;
        if (c.snpfile.empty()) {
            stattype.push_back("");
        } else {
            stattype.push_back("-all");
            //    stattype.push_back("-derConsensus");
            stattype.push_back("-infsites");
        }

        for (unsigned int j=0; j < statname.size(); j++) {
            if (summarize==0) {
                if (html) printf("</td><td>\n");
                printf("\t%s", statname[j].c_str());
            }
            if (statname[j] != "inf_sites") {
                for (unsigned int k=0; k < stattype.size(); k++) {
                    for (int i=1; i <= summarize; i++) {
                        if (getMean==i) {
                            if (html) printf("</td><td>\n");
                            printf("\t%s%s_mean",
                                   statname[j].c_str(), stattype[k].c_str());
                        } else if (getStdev==i) {
                            if (html) printf("</td><td>\n");
                            printf("\t%s%s_stdev",
                                   statname[j].c_str(), stattype[k].c_str());
                        } else if (getQuantiles==i) {
                            for (unsigned int l=0; l < quantiles.size(); l++) {
                                if (html) printf("</td><td>\n");
                                printf("\t%s%s_quantile_%.3f",
                                       statname[j].c_str(),
                                       stattype[k].c_str(),
                                       quantiles[l]);
                            }
                        }
                    }
                }
            }
        }
        printf("\n");
        if (html) printf("</td></tr>\n");
    } else if (html) printf("<table valing=\"top\">\n");



    for (unsigned int i=0; i < groupfiles.size(); i++) {
        ifstream in(groupfiles[i].c_str());
        string line;
        data.group.resize(groupfiles.size());
        if (in.is_open()) {
            while ( getline(in, line) ) {
                data.group[i].insert(line);
            }
            in.close();
        } else {
            fprintf(stderr, "Error opening %s.\n", groupfiles[i].c_str());
            return 1;
        }
    }

    /*    if (!c.group_noroot_file.empty()) {
        ifstream in(c.group_noroot_file.c_str());
        string line;
        if (in.is_open()) {
            while ( getline(in, line) ) {
                data.group_noroot.insert(line);
            }
            in.close();
        } else {
            fprintf(stderr, "Error opening %s.\n", c.group_noroot_file.c_str());
            return 1;
        }
        }*/

    if (c.bedfile.empty()) {
        summarizeRegion(&c, c.region.empty() ? NULL : c.region.c_str(),
                        haps, statname, data);
    } else {
        CompressStream bedstream(c.bedfile.c_str());
        char *line;
        vector<string> token;
        char *regionStr;
        if (!bedstream.stream) {
            fprintf(stderr, "error reading %s\n", c.bedfile.c_str());
            return 1;
        }
        while ((line = fgetline(bedstream.stream))) {
            split(line, '\t', token);
            if (token.size() < 3) {
                fprintf(stderr, "expected at least 3 files in %s\n",
                        c.bedfile.c_str());
                return 1;
            }
            regionStr = new char[token[0].size()+token[1].size()+
                                 token[2].size()+3];
            int start = atoi(token[1].c_str());
            int end = atoi(token[2].c_str());
            sprintf(regionStr, "%s:%i-%i", token[0].c_str(), start+1, end);
            summarizeRegion(&c, regionStr, haps, statname, data);
            delete [] regionStr;
        }
        bedstream.close();
    }
    if (html) printf("</table>\n</html>\n");

    return 0;
}
