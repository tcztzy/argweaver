// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "local_tree.h"
#include "sequences.h"
#include "trans.h"
#include "total_prob.h"

using namespace std;

namespace argweaver {


double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees, int start_coord, int end_coord)
{
    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();
    if (start_coord < trees->start_coord)
        start_coord = trees->start_coord;
    if (end_coord < 0 || end_coord > trees->end_coord)
        end_coord = trees->end_coord;

    // special case for truck genealogies
    if (trees->nnodes < 3)
        return lnl += log(.25) * (end_coord - start_coord);

    // get sequences for trees
    char *seqs[nseqs];
    for (int j=0; j<nseqs; j++)
        seqs[j] = sequences->seqs[trees->seqids[j]];

    int end = trees->start_coord;
    int mu_idx = 0, rho_idx = 0;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        end = start + it->blocklen;
        if (end <= start_coord) continue;
        if (start >= end_coord) break;
        if (start < start_coord) start = start_coord;
        if (end > end_coord) end = end_coord;
        LocalTree *tree = it->tree;
        ArgModel local_model;

        //note: this is approximate, uses mu/rho from center of block
        model->get_local_model((start+end)/2, local_model, &mu_idx, &rho_idx);
        lnl += likelihood_tree(tree, &local_model, seqs, nseqs, start, end);
    }

    return lnl;
}


    // TODO: This fills in compressed sites with A's... should
    // take mask into account!
// NOTE: trees should be uncompressed and sequences compressed
//start_coord and end_coord uncompressed, 0 based
double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees,
                           const SitesMapping* sites_mapping,
                           const TrackNullValue *maskmap_uncompressed,
                           int start_coord, int end_coord)
{
    if (!sites_mapping)
        return calc_arg_likelihood(model, sequences, trees, start_coord, end_coord);

    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();
    const char default_char = 'A';

    if (start_coord < trees->start_coord)
        start_coord = trees->start_coord;
    if (end_coord < 0 || end_coord > trees->end_coord)
        end_coord = trees->end_coord;

    // special case for truck genealogies
    if (trees->nnodes < 3)
        return lnl += log(.25) * (end_coord - start_coord);

    int end = trees->start_coord;
    int mu_idx = 0;
    int rho_idx = 0;
    int mask_pos=0;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        end = start + it->blocklen;
        if (end <= start_coord) continue;
        if (start >= end_coord) break;
        if (start < start_coord)
            start = start_coord;
        if (end > end_coord)
            end = end_coord;
        int blocklen = end - start;
        LocalTree *tree = it->tree;

        // get sequences for trees
        char *seqs[nseqs];
        char *matrix = new char [blocklen*nseqs];
        for (int j=0; j<nseqs; j++)
            seqs[j] = &matrix[j*blocklen];

        // find first site within this block
        unsigned int i2 = 0;

        // copy sites into new alignment
        for (int i=start; i<end; i++) {
            while (i2 < sites_mapping->all_sites.size() &&
                   sites_mapping->all_sites[i2] < i)
                i2++;
            if (i2 < sites_mapping->all_sites.size() &&
                i == sites_mapping->all_sites[i2]) {
                // copy site
                for (int j=0; j<nseqs; j++)
                    seqs[j][i-start] = sequences->seqs[trees->seqids[j]][i2];
            } else {
                // copy non-variant site
                char c=default_char;
                if (maskmap_uncompressed->find(i, &mask_pos))
                    c='N';
                for (int j=0; j<nseqs; j++)
                    seqs[j][i-start] = c;
            }
        }

        ArgModel local_model;
        model->get_local_model((start+end)/2, local_model,
                               &mu_idx, &rho_idx);
        lnl += likelihood_tree(tree, &local_model, seqs, nseqs, 0, end-start);

        delete [] matrix;
    }

    return lnl;
}

//=============================================================================
// ARG prior


// The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
// with population size 'n'
double log_prob_coal_counts(int a, int b, double t, double n)
{
    double lnC1 = 0.0, lnC2=0.0, C3=1.0;

    for (int y=0; y<b; y++)
        lnC1 += log((b+y)*(a-y)/double(a+y));
    assert(!isnan(lnC1));
    double s = -b*(b-1)*t/2.0/n;

    for (int k=b+1; k<a+1; k++) {
        const double k1 = double(k - 1);
        double val;
        lnC2 += log(double(b+k1)*(a-k1)/(a+k1)/(k-b));
        C3 *= -1.0;
        val = -k*k1*t/2.0/n + lnC2 + log((2.0*k-1)/double(k1+b));
        //        s += exp(-k*k1*t/2.0/n) * (2*k-1) / double(k1+b) * C2;
        if (1) {
            double loga = s;
            double logc = val;
            if (logc > loga) {
                double tmp = logc;
                logc = loga;
                loga = tmp;
            }
            s = loga + log(1.0 + C3*exp(logc - loga));
        }
        /*        if (s > val)
            s += log(1.0 + C3*exp(val-s));
            else s = val + log(1.0 + C3*exp(s-val));*/
        /*        double s2 = log(exp(s) + C3*exp(val));
        assert(!isnan(s2));
        s = s2;*/
    }

    for (int i=2; i<=b; i++)
        s -= log(i);
    assert(!isnan(s));
    assert(!isinf(s));
    return s + lnC1;
}


// The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
// with population size 'n'
double prob_coal_counts(int a, int b, double t, double n)
{
    double C = 1.0;

    for (int y=0; y<b; y++)
        C *= (b+y)*(a-y)/double(a+y);

    double s = exp(-b*(b-1)*t/2.0/n) * C;

    for (int k=b+1; k<a+1; k++) {
        const double k1 = double(k - 1);
        C *= double(b+k1)*(a-k1)/(a+k1)/(b-k);
        s += exp(-k*k1*t/2.0/n) * (2*k-1) / double(k1+b) * C;
    }

    for (int i=1; i<=b; i++)
        s /= i;

    return s;
}

double calc_tree_prior(const ArgModel *model, const LocalTree *tree,
                       LineageCounts &lineages)
{
    lineages.count(tree, model->pop_tree);
    double lnl = 0.0;

    for (int pop=0; pop < model->num_pops(); pop++) {
        for (int i=0; i<model->ntimes-1; i++) {
            int a;
            if (i == 0) {
                assert((lineages.ncoals_pop[pop][0] +
                        lineages.nbranches_pop[pop][0])%2 == 1);
                a = (lineages.ncoals_pop[pop][0] +
                     lineages.nbranches_pop[pop][0])/2;
                assert( a >= 0);
            } else {
                a = lineages.nbranches_pop[pop][2*i-1];
            }
            if (a <= 1) break;
            int b = lineages.nbranches_pop[pop][2*i];
            double t = model->coal_time_steps[2*i];
            if (i > 0) t += model->coal_time_steps[2*i-1];
            lnl += log_prob_coal_counts(a, b, t, 2.0 * model->popsizes[pop][2*i]);
            assert(!isnan(lnl));
            assert(!isinf(lnl));
        }
    }

    if (model->pop_tree != NULL) {
        for (int i=0; i < tree->nnodes; i++) {
            int start_time = tree->nodes[i].age;
            int end_time = ( tree->nodes[i].parent == -1 ? model->ntimes - 1 :
                             tree->nodes[tree->nodes[i].parent].age );
            if (start_time != end_time)
                lnl += log(model->path_prob(tree->nodes[i].pop_path,
                                            start_time, end_time));
        }
    }

    // TODO: top prior

    return lnl;
}


void calc_coal_rates_spr(const ArgModel *model, const LocalTree *tree,
                         const Spr &spr, const LineageCounts &lineages,
                         double *coal_rates)
{
    int broken_age = tree->nodes[tree->nodes[spr.recomb_node].parent].age;

    for (int i=spr.recomb_time*2; i<=2*spr.coal_time; i++) {
        int pop_time = (i+1)/2;
        int spr_pop = model->get_pop(spr.pop_path, pop_time);
        int recomb_parent_pop =
            model->get_pop(tree->nodes[spr.recomb_node].pop_path, pop_time);
        int nbranches = lineages.nbranches_pop[spr_pop][i]
            - int((!model->smc_prime) && i/2 < broken_age && spr_pop == recomb_parent_pop);
        coal_rates[i] = model->coal_time_steps[i] * nbranches /
            (2.0 * model->popsizes[spr_pop][i]);
	if (coal_rates[i] < 0) {
	    assert(0);
	}
    }
}


double calc_spr_prob(const ArgModel *model, const LocalTree *tree,
                     const Spr &spr, LineageCounts &lineages,
                     double treelen, double **num_coal, double **num_nocoal,
                     double coal_weight, int lineages_counted)
{
    if (spr.recomb_node == tree->root) {
        assert(spr.coal_node == tree->root);
        assert(model->pop_tree != NULL);
        return (model->path_prob(spr.pop_path, spr.recomb_time, spr.coal_time));
    }
    assert(spr.recomb_node != tree->root);
    const LocalNode *nodes = tree->nodes;
    const int root_age = nodes[tree->root].age;

    // get tree length, if it is not already given
    if (treelen < 0)
        treelen = get_treelen(tree, model->times, model->ntimes, false);

    // get lineage counts
    if (!lineages_counted) {
        lineages.count(tree, model->pop_tree);
        lineages.nrecombs[root_age]--;
    }

    double lnl = 0.0;

    // probability of recombination location in tree (from eq 5)
    int k = spr.recomb_time;
    if (model->smc_prime) {
        int pop = model->get_pop(spr.pop_path, k);
        int nrecomb = lineages.ncoals_pop[pop][k];
        double blen_above=0, blen_below=0;
        if (k == root_age) {
            nrecomb--;
        } else {
            assert(k < root_age);
            blen_above = model->coal_time_steps[2*k]
                * (double)lineages.nbranches_pop[pop][2*k];
        }
        if (k > 0) {
            blen_below = model->coal_time_steps[2*k-1]
                * (double)lineages.nbranches_pop[pop][2*k-1];
        }
        lnl += log((blen_above + blen_below) /
                   (treelen * (double)nrecomb));

    } else {
        lnl += log(lineages.nbranches[k] * model->time_steps[k] /
                   (lineages.nrecombs[k] *
                    (treelen + model->time_steps[nodes[tree->root].age])));
    }

    // probability of re-coalescence
    double coal_rates[2*model->ntimes];
    calc_coal_rates_spr(model, tree, spr, lineages, coal_rates);
    int j = spr.coal_time;
    int broken_age = nodes[nodes[spr.recomb_node].parent].age;

    // probability of recoalescence on chosen branch
    int recomb_parent_pop = model->get_pop(nodes[spr.recomb_node].pop_path,
                                           spr.coal_time);
    int spr_pop = model->get_pop(spr.pop_path, spr.coal_time);
    int ncoals_j = lineages.ncoals_pop[spr_pop][j]
        - int((!model->smc_prime) && j <= broken_age && spr_pop == recomb_parent_pop)
        - int((!model->smc_prime) && j == broken_age && spr_pop == recomb_parent_pop);
    lnl -= log(ncoals_j);
    //printf("     ncoals_j = %i\n", ncoals_j);

    // probability of recoalescing in chosen time interval
    if (j < model->ntimes - 2)
        lnl += log(1.0 - exp(- coal_rates[2*j] -
                             (j>k ? coal_rates[2*j-1] : 0.0)));

    // probability of not coalescing before time interval j
    for (int m=2*k; m<2*j-1; m++)
        lnl -= coal_rates[m];

    // probability of spr population path
    lnl += log(model->path_prob(spr.pop_path, k, j));

    // add to counts of coal/no_coal events, if not null
    if (num_coal != NULL) {
	if (k == j)
            num_coal[spr_pop][2*j] += coal_weight;
	else {
            // Note: these heuristics don't make a lot of sense but are
            // only used for reporting counts to monitor MCMC estimate of popsize
            num_coal[spr_pop][2*j] += coal_weight / 2.0;
            num_coal[spr_pop][2*j-1] += coal_weight / 2.0;
            num_nocoal[spr_pop][2*j-1] += coal_weight / 2.0;
            for (int i=2*k; i < 2*j-1; i++) {
                int pop = model->get_pop(spr.pop_path, (k+1)/2);
                num_nocoal[pop][i] += coal_weight;
            }
	}
    }
    assert(!isinf(lnl));
    assert(!isnan(lnl));
    return lnl;
}


double calc_self_recomb_prob(const ArgModel *model, LocalTree *tree,
                             LineageCounts &lineages, double treelen) {
    double lnprob=-INFINITY;
    for (int i=0; i < tree->nnodes; i++) {
        if (i == tree->root) continue;
        int parent_age = tree->nodes[tree->nodes[i].parent].age;
        for (int j = tree->nodes[i].age; j <= parent_age; j++) {
            for (int k=j; k <= parent_age; k++) {
                Spr spr(i,j,i,k,tree->nodes[i].pop_path);
                lnprob = logadd(lnprob,
                                calc_spr_prob(model, tree, spr,
                                              lineages, treelen,
                                              NULL, NULL, 1.0, true));
            }
        }
    }
    return lnprob;
}


// calculate the probability of an ARG given the model parameters
double calc_arg_prior(const ArgModel *model, const LocalTrees *trees,
		      double **num_coal, double **num_nocoal,
                      int start_coord, int end_coord)
{
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes, model->num_pops());

    if (num_coal != NULL) {
	assert(num_nocoal != NULL);
        for (int pop=0; pop < model->num_pops(); pop++) {
            for (int i=0; i < 2*model->ntimes-1; i++)
                num_coal[pop][i] = num_nocoal[pop][i] = 0;
        }
    }
    if (start_coord < trees->start_coord)
        start_coord = trees->start_coord;
    if (end_coord < 0 || end_coord > trees->end_coord)
        end_coord = trees->end_coord;

    // first tree prior
    //    lnl += calc_tree_prior(model, trees->front().tree, lineages);
    //    printLog(LOG_MEDIUM, "tree_prior: %f\n", lnl);

    int end = trees->start_coord;
    int mu_idx = 0, rho_idx = 0;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
        int start=end;
        end += it->blocklen;
        if (end <= start_coord) {++it; continue;}
        if (start >= end_coord) break;
        if (start < start_coord)
            start = start_coord;
        if (end > end_coord)
            end = end_coord;
        int blocklen = start - end;
        LocalTree *tree = it->tree;
        double treelen = get_treelen(tree, model->times, model->ntimes, false);
        ArgModel local_model;
        model->get_local_model((start+end)/2, local_model, &mu_idx, &rho_idx);
        lineages.count(tree, model->pop_tree);

        // not sure what this is for but it is only used for non-SMC' calcs
        lineages.nrecombs[tree->nodes[tree->root].age]--;

        // calculate probability P(blocklen | T_{i-1})
        double recomb_rate = max(local_model.rho * treelen, local_model.rho);
        double self_recomb_prob = exp(calc_self_recomb_prob(&local_model, tree, lineages, treelen));
        assert(self_recomb_prob >= 0 && self_recomb_prob < 1.0);
        recomb_rate = max(local_model.rho * treelen, local_model.rho)*(1.0-self_recomb_prob);
        if (end < end_coord) {
            // not last block
            // probability of recombining after blocklen
	    lnl += log(recomb_rate) - recomb_rate * blocklen;

            // get SPR move information
            ++it;
            const Spr *spr = &it->spr;
            lnl += calc_spr_prob(&local_model, tree, *spr, lineages, treelen,
				 num_coal, num_nocoal, 1.0, true);

        } else {
            // last block
            // probability of not recombining after blocklen
            lnl += - recomb_rate * blocklen;
            ++it;
        }
    }
    assert(!isnan(lnl));
    assert(!isinf(lnl));
    return lnl;
}

double calc_arg_prior_recomb_integrate(const ArgModel *model,
                                       const LocalTrees *trees,
                                       double **num_coal, double **num_nocoal) {
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes, model->num_pops());

    if (num_coal != NULL) {
	assert(num_nocoal != NULL);
	for (int i=0; i < 2*model->ntimes-1; i++)
	    num_coal[i] = num_nocoal[i] = 0;
    }

    // first tree prior- ignore for now
    lnl += calc_tree_prior(model, trees->front().tree, lineages);
    //    printLog(LOG_MEDIUM, "tree_prior: %f\n", lnl);

    int rho_idx = 0;
    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
        end += it->blocklen;
        LocalTree *tree = it->tree;
        int blocklen = it->blocklen;
        double treelen = get_treelen(tree, model->times, model->ntimes, false);

        // calculate probability P(blocklen | T_{i-1})
        double recomb_rate = max(model->get_local_rho(trees->start_coord,
                                                      &rho_idx)
                                 * treelen, model->rho);

        //for single site, probability of no recomb
        double pr_no_recomb = exp(-recomb_rate);
        double pr_recomb = 1.0 - pr_no_recomb;
        const int root_age = tree->nodes[tree->root].age;
        lineages.count(tree, model->pop_tree);
        lineages.nrecombs[root_age]--;
        if (end >= trees->end_coord)
            blocklen++;
        if (blocklen > 1) {

            double recomb_sum = 0.0;
            for (int i=0; i < tree->nnodes; i++) {
                const LocalNode *node = &(tree->nodes[i]);
                if (node->parent != -1) {
                    int maxage = tree->nodes[node->parent].age;
                    assert(node->age <= maxage);
                    for (int age=node->age; age <= maxage; age++) {
                        Spr spr(i, age, node->parent, maxage, node->pop_path);
                        double val =
                            exp(calc_spr_prob(model, tree, spr, lineages,
                                              treelen, NULL, NULL, 0, true));
                        recomb_sum += val;
                        // printf("%i\t%i\t%f\t%f\n", i, age, val, recomb_sum);

                    }
                }
            }
            lnl += ((double)blocklen-1.0) * log(pr_no_recomb +
                                                pr_recomb * recomb_sum);
        }

        if (end < trees->end_coord) {
            // not last block, add probability of any recomb that results in
            // same topology as sampled SPR

            // get SPR move information
            ++it;
            const Spr *real_spr = &it->spr;
            int node = real_spr->recomb_node;
            int parent = tree->nodes[node].parent;
            int max_age = min(tree->nodes[parent].age,
                              real_spr->coal_time);
            assert(tree->nodes[node].age <= max_age);
            assert(real_spr->recomb_time >= tree->nodes[node].age &&
                   real_spr->recomb_time <= max_age);
            double recomb_sum = 0.0;
            int target_path = model->consistent_path(tree->nodes[node].pop_path,
                                                     real_spr->pop_path,
                                                     tree->nodes[node].age,
                                                     real_spr->recomb_time,
                                                     real_spr->coal_time);
            for (int age=tree->nodes[node].age; age <= max_age; age++) {
                Spr spr(node, age, real_spr->coal_node, real_spr->coal_time,
                        target_path);
                double val = exp(calc_spr_prob(model, tree, spr, lineages,
                                               treelen, num_coal, num_nocoal,
                                               age == real_spr->recomb_time
                                               ? 1.0 : 0.0, true));
                recomb_sum += val;
            }
            int sib = tree->nodes[parent].child[0] == node ?
                tree->nodes[parent].child[1] : tree->nodes[parent].child[0];
            bool possible_no_recomb = false;
            if (real_spr->coal_node == sib) {
                if (model->paths_equal(real_spr->pop_path,
                                       tree->nodes[node].pop_path,
                                       real_spr->recomb_time,
                                       real_spr->coal_time)) {
                    if (real_spr->coal_time == tree->nodes[parent].age)
                        possible_no_recomb = true;
                    for (int age=tree->nodes[sib].age; age <= max_age; age++) {
                        Spr spr(sib, age, node, real_spr->coal_time,
                                tree->nodes[sib].pop_path);
                        recomb_sum += exp(calc_spr_prob(model, tree, spr, lineages,
                                                        treelen, num_coal, num_nocoal,
                                                        0, true));
                    }
                }
            } else if (real_spr->coal_node == parent) {
                int path = model->consistent_path(tree->nodes[node].pop_path,
                                                  tree->nodes[parent].pop_path,
                                                  tree->nodes[node].age,
                                                  tree->nodes[parent].age,
                                                  real_spr->coal_time);
                if (model->paths_equal(path, real_spr->pop_path,
                                       real_spr->recomb_time,
                                       real_spr->coal_time)) {
                    if (real_spr->coal_time == tree->nodes[parent].age)
                        possible_no_recomb = true;
                    path = model->consistent_path(tree->nodes[sib].pop_path,
                                                  tree->nodes[parent].pop_path,
                                                  tree->nodes[sib].age,
                                                  tree->nodes[parent].age,
                                                  real_spr->coal_time);
                    for (int age=tree->nodes[sib].age; age <= max_age; age++) {
                        Spr spr(sib, age, parent, real_spr->coal_time, path);
                        recomb_sum += exp(calc_spr_prob(model, tree, spr, lineages,
                                                        treelen, num_coal, num_nocoal,
                                                        0, true));
                    }
                }
            }
            lnl += log(pr_recomb * recomb_sum +
                       (possible_no_recomb ? pr_no_recomb : 0.0));
	    if (isinf(lnl))
		assert(0);
        } else {
            ++it;
        }
    }
    assert(!isnan(lnl));
    assert(!isinf(lnl));
    return lnl;
}



// calculate the probability of the sequences given an ARG
double calc_arg_joint_prob(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees)
{
    return calc_arg_likelihood(model, sequences, trees) +
        calc_arg_prior(model, trees, NULL, NULL);
}



//=============================================================================
// C interface
/*
extern "C" {

double arghmm_likelihood(LocalTrees *trees,
                         double *times, int ntimes,
                         double mu,
                         char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    return calc_arg_likelihood(&model, &sequences, trees);
}


double arghmm_likelihood_parsimony(LocalTrees *trees,
                                   double *times, int ntimes,
                                   double mu,
                                   char **seqs, int nseqs, int seqlen)
{

    // setup model, local trees, sequences
//    ArgModel model(ntimes, times, NULL, 0.0, mu);
//    Sequences sequences(seqs, nseqs, seqlen);
//    return calc_arg_likelihood_parsimony(&model, &sequences, trees);

    abort();
    return 0.0;
}


double arghmm_prior_prob(LocalTrees *trees,
                         double *times, int ntimes, double *popsizes,
                         double rho)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    return calc_arg_prior(&model, trees, NULL, NULL);
}


    // NOT CALLED
double arghmm_tree_prior_prob(LocalTrees *trees,
                              double *times, int ntimes, double *popsizes)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, 0.0, 0.0);
    LineageCounts lineages(ntimes, model->num_pops());

    //printf("%f\n", calc_tree_prior_approx(&model, trees->front().tree,
    //                                      lineages));

    return calc_tree_prior(&model, trees->front().tree, lineages);
}


double arghmm_joint_prob(LocalTrees *trees,
                         double *times, int ntimes, double *popsizes,
                         double mu, double rho,
                         char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    return calc_arg_joint_prob(&model, &sequences, trees);
}



} // extern "C"

*/
} // namespace argweaver

