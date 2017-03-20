//=============================================================================
// sample full ARGs
//

#ifndef ARGWEAVER_SAMPLE_ARG_H
#define ARGWEAVER_SAMPLE_ARG_H

// c++ includes
#include <vector>

// arghmm includes
#include "local_tree.h"
#include "model.h"
#include "sequences.h"


namespace argweaver {

using namespace std;

void sample_arg_seq(const ArgModel *model, Sequences *sequences,
                    LocalTrees *trees, bool random=false, int num_buildup=1);

void resample_arg(const ArgModel *model, Sequences *sequences,
                  LocalTrees *trees);

void resample_arg_all(const ArgModel *model, Sequences *sequences,
                      LocalTrees *trees, double prob_path_switch);

void resample_arg_leaf(const ArgModel *model, Sequences *sequences,
                       LocalTrees *trees, int node);

void resample_arg_random_leaf(const ArgModel *model, Sequences *sequences,
			      LocalTrees *trees);

bool resample_arg_mcmc(const ArgModel *model, Sequences *sequences,
                       LocalTrees *trees, double heat=1.0);

void resample_arg_mcmc_all(const ArgModel *model, Sequences *sequences,
                           LocalTrees *trees, bool do_leaf,
                           int window, int niters, double heat=1.0,
                           bool no_resample_mig=false);

void resample_arg_climb(const ArgModel *model, Sequences *sequences,
                        LocalTrees *trees, double recomb_preference);

void remax_arg(const ArgModel *model, const Sequences *sequences,
               LocalTrees *trees, int nremove=1);

/*void cond_sample_arg_seq(const ArgModel *model, Sequences *sequences,
                         LocalTrees *trees,
                         LocalTree *start_tree, LocalTree *end_tree,
                         const vector<int> &full_seqids);

void sample_arg_seq_region(const ArgModel *model, const Sequences *sequences,
                           LocalTrees *trees, int region_start, int region_end);
*/

double resample_arg_region(
    const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int region_start, int region_end, int niters,
    bool open_ended=true, double heat=1.0, int time_interval=-1, int hap=-1);

double resample_arg_cut(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    int window_start=-1, int window_end=-1);

double resample_arg_cut(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    int window, int step, int niters);

double resample_arg_regions(
    const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int window, int niters=1,
    double heat=1.0, int time_interval=-1, int hap=-1);

int resample_arg_by_time_and_hap(
    const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int time_interval, int hap);

void resample_migrates(ArgModel *model,
                       const LocalTrees *trees,
                       vector<Spr> &invisible_recombs);

} // namespace argweaver

#endif // ARGWEAVER_SAMPLE_ARG_H
