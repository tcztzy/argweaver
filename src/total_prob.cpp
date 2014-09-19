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
                           const LocalTrees *trees)
{
    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();

    // special case for truck genealogies
    if (trees->nnodes < 3)
        return lnl += log(.25) * sequences->length();

    // get sequences for trees
    char *seqs[nseqs];
    for (int j=0; j<nseqs; j++)
        seqs[j] = sequences->seqs[trees->seqids[j]];

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        end = start + it->blocklen;
        LocalTree *tree = it->tree;

        lnl += likelihood_tree(tree, model, seqs, nseqs, start, end);
    }

    return lnl;
}


// NOTE: trees should be uncompressed and sequences compressed
double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees,
                           const SitesMapping* sites_mapping)
{
    if (!sites_mapping)
        return calc_arg_likelihood(model, sequences, trees);

    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();
    const char default_char = 'A';

    // special case for truck genealogies
    if (trees->nnodes < 3)
        return lnl += log(.25) * sequences->length();

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        int blocklen = it->blocklen;
        end = start + blocklen;
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
                for (int j=0; j<nseqs; j++)
                    seqs[j][i-start] = default_char;
            }
        }

        lnl += likelihood_tree(tree, model, seqs, nseqs, 0, end-start);

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
    lineages.count(tree);
    int nleaves = tree->get_num_leaves();
    double lnl = 0.0;

    for (int i=0; i<model->ntimes-1; i++) {
        int a = (i == 0 ? nleaves : lineages.nbranches[i-1]);
        int b = lineages.nbranches[i];
	double t = model->coal_time_steps[2*i];
	if (i > 0) t += model->coal_time_steps[2*i-1];
        lnl += log_prob_coal_counts(a, b, t, 2.0 * model->popsizes[2*i]);
        assert(!isnan(lnl));
        assert(!isinf(lnl));
    }

    // TODO: top prior

    return lnl;
}


double calc_tree_prior_approx(const ArgModel *model, const LocalTree *tree,
                              LineageCounts &lineages)
{
    lineages.count(tree);
    double lnl = 0.0;

    for (int i=0; i<tree->nnodes; i++) {
        if (tree->nodes[i].is_leaf())
            continue;
        int time = tree->nodes[i].age;

        // remove lineage counts
        for (int j=0; j<time; j++)
            lineages.nbranches[j]--;

        lnl += log(calc_state_priors(time, lineages.nbranches, lineages.ncoals,
                                     0, model->popsizes, model->coal_time_steps,
                                     model->ntimes));
        lnl += log(lineages.ncoals[time]);
    }

    if (isnan(lnl))
        lnl = 0.0;

    return lnl;
}


void calc_coal_rates_full_tree(const ArgModel *model, const LocalTree *tree,
                               const Spr &spr, const LineageCounts &lineages,
                               double *coal_rates)
{
    int broken_age = tree->nodes[tree->nodes[spr.recomb_node].parent].age;

    for (int i=0; i<2*model->ntimes-1; i++) {
        int nbranches = lineages.nbranches[i/2] - int(i/2 < broken_age);
        coal_rates[i] = model->coal_time_steps[i] * nbranches /
            (2.0 * model->popsizes[i]);
	if (coal_rates[i] < 0) {
	    assert(0);
	}
    }
}


double calc_spr_prob(const ArgModel *model, const LocalTree *tree,
                     const Spr &spr, LineageCounts &lineages,
                     double treelen, double *num_coal, double *num_nocoal,
                     double coal_weight, int lineages_counted)
{
    assert(spr.recomb_node != tree->root);
    const LocalNode *nodes = tree->nodes;
    const int root_age = nodes[tree->root].age;

    // get tree length, if it is not already given
    if (treelen < 0)
        treelen = get_treelen(tree, model->times, model->ntimes, false);
    const double treelen_b = treelen + model->time_steps[nodes[tree->root].age];

    // get lineage counts
    if (!lineages_counted) {
        lineages.count(tree);
        lineages.nrecombs[root_age]--;
    }

    double lnl = 0.0;

    // probability of recombination location in tree (from eq 5)
    int k = spr.recomb_time;
    lnl += log(lineages.nbranches[k] * model->time_steps[k] /
                   (lineages.nrecombs[k] * treelen_b));

    // probability of re-coalescence
    double coal_rates[2*model->ntimes];
    calc_coal_rates_full_tree(model, tree, spr, lineages, coal_rates);
    int j = spr.coal_time;
    int broken_age = nodes[nodes[spr.recomb_node].parent].age;

    // probability of recoalescence on chosen branch
    int ncoals_j = lineages.ncoals[j]
        - int(j <= broken_age) - int(j == broken_age);
    lnl -= log(ncoals_j);
    //printf("     ncoals_j = %i\n", ncoals_j);

    // probability of recoalescing in chosen time interval
    if (j < model->ntimes - 2)
        lnl += log(1.0 - exp(- coal_rates[2*j] -
                             (j>k ? coal_rates[2*j-1] : 0.0)));

    // probability of not coalescing before time interval j
    for (int m=2*k; m<2*j-1; m++)
        lnl -= coal_rates[m];

    // add to counts of coal/no_coal events, if not null
    if (num_coal != NULL) {
	if (k == j)
            num_coal[2*j] += coal_weight;
	else {
            num_coal[2*j] += coal_weight / 2.0;
            num_coal[2*j-1] += coal_weight / 2.0;
            for (int i=2*k; i < 2*j-1; i++)
                num_nocoal[i] += coal_weight;
            num_nocoal[2*j-1] += coal_weight / 2.0;
	}
    }
    assert(!isinf(lnl));
    assert(!isnan(lnl));
    return lnl;
}



// calculate the probability of an ARG given the model parameters
double calc_arg_prior(const ArgModel *model, const LocalTrees *trees,
		      double *num_coal, double *num_nocoal)
{
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes);

    if (num_coal != NULL) {
	assert(num_nocoal != NULL);
	for (int i=0; i < 2*model->ntimes-1; i++)
	    num_coal[i] = num_nocoal[i] = 0;
    }

    // first tree prior
    //    lnl += calc_tree_prior(model, trees->front().tree, lineages);
    //    printLog(LOG_MEDIUM, "tree_prior: %f\n", lnl);

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
        end += it->blocklen;
        LocalTree *tree = it->tree;
        int blocklen = it->blocklen;
        double treelen = get_treelen(tree, model->times, model->ntimes, false);

        // calculate probability P(blocklen | T_{i-1})
        double recomb_rate = max(model->get_local_rho(trees->start_coord)
                                 * treelen, model->rho);

        if (end < trees->end_coord) {
            // not last block
            // probability of recombining after blocklen
       //     lnl += log(recomb_rate) - recomb_rate * blocklen;
            lnl += -recomb_rate * (blocklen - 1) + log(1.0 - exp(-recomb_rate));

            // get SPR move information
            ++it;
            const Spr *spr = &it->spr;
            lnl += calc_spr_prob(model, tree, *spr, lineages, treelen,
				 num_coal, num_nocoal);

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
                                       double *num_coal, double *num_nocoal) {
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes);

    if (num_coal != NULL) {
	assert(num_nocoal != NULL);
	for (int i=0; i < 2*model->ntimes-1; i++)
	    num_coal[i] = num_nocoal[i] = 0;
    }

    // first tree prior- ignore for now
    //    lnl += calc_tree_prior(model, trees->front().tree, lineages);
    //    printLog(LOG_MEDIUM, "tree_prior: %f\n", lnl);

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
        end += it->blocklen;
        LocalTree *tree = it->tree;
        int blocklen = it->blocklen;
        double treelen = get_treelen(tree, model->times, model->ntimes, false);

        // calculate probability P(blocklen | T_{i-1})
        double recomb_rate = max(model->get_local_rho(trees->start_coord)
                                 * treelen, model->rho);

        //for single site, probability of no recomb
        double pr_no_recomb = exp(-recomb_rate);
        double pr_recomb = 1.0 - pr_no_recomb;
        const int root_age = tree->nodes[tree->root].age;
        lineages.count(tree);
        lineages.nrecombs[root_age]--;
        if (end >= trees->end_coord)
            blocklen++;
        if (blocklen > 1) {

            if (0) {
               // write_newick_tree(stdout, tree, NULL, model->times, -1, true);
              //  printf("\n");
                double tempsum=0.0;
                for (int i=0; i < tree->nnodes; i++) { // i is recomb node
                  //  printf("Recomb node %i\n", i);
                    if (tree->nodes[i].parent == -1) continue;
                    int i_start = tree->nodes[i].age;
                    int i_end = tree->nodes[tree->nodes[i].parent].age;
                    for (int recomb_time=i_start; recomb_time <= i_end;
                         recomb_time++) {
                        for (int j=0; j < tree->nnodes; j++) { //j is coal node
                            if (i==j) continue;
                            bool is_sib= (j!=tree->root &&
                                          ( tree->nodes[j].parent ==
                                            tree->nodes[i].parent) );
                            int j_start = tree->nodes[j].age;
                            int j_end = tree->nodes[j].parent == -1
                                ? model->ntimes-2 :
                                tree->nodes[tree->nodes[j].parent].age;
                            if (is_sib) j_end--;
                            for (int coal_time = j_start; coal_time <= j_end; coal_time++) {
                                if (coal_time < recomb_time) continue;
                                Spr spr(i, recomb_time, j, coal_time);
                                double val =
                                    exp(calc_spr_prob(model, tree, spr, lineages,
                                                      treelen, NULL, NULL, 0,
                                                      true));
                                tempsum += val;
                            }
                        }
                    }
                }
                if (fabs(tempsum-1.0) > 0.00001)
                  printf("tempsum=%f\n", tempsum);
            }

            double recomb_sum = 0.0;
            for (int i=0; i < tree->nnodes; i++) {
                const LocalNode *node = &(tree->nodes[i]);
                if (node->parent != -1) {
                    int maxage = tree->nodes[node->parent].age;
                    assert(node->age <= maxage);
                    for (int age=node->age; age <= maxage; age++) {
                        Spr spr(i, age, node->parent, maxage);
                        double val =
                            exp(calc_spr_prob(model, tree, spr, lineages,
                                              treelen, NULL, NULL, 0, true));
                        recomb_sum += val;
                        //                        printf("%i\t%i\t%f\t%f\n", i, age, val, recomb_sum);

                    }
                }
            }
            //lnl += ((double)blocklen-1.0) * log(pr_no_recomb);
            lnl += ((double)blocklen-1.0) * log(pr_no_recomb +
                                                pr_recomb * recomb_sum);
            //printf("%e\t%e\n", (blocklen-1.0)*log(pr_no_recomb + pr_recomb * recomb_sum),
              //      log(recomb_rate) - recomb_rate * blocklen);
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
            for (int age=tree->nodes[node].age; age <= max_age; age++) {
                Spr spr(node, age, real_spr->coal_node, real_spr->coal_time);
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
                if (real_spr->coal_time == tree->nodes[parent].age)
                    possible_no_recomb = true;
                for (int age=tree->nodes[sib].age; age <= max_age; age++) {
                    Spr spr(sib, age, node, real_spr->coal_time);
                    recomb_sum += exp(calc_spr_prob(model, tree, spr, lineages,
                                                    treelen, num_coal, num_nocoal,
                                                    0, true));
                }
            } else if (real_spr->coal_node == parent) {
                if (real_spr->coal_time == tree->nodes[parent].age)
                    possible_no_recomb = true;
              for (int age=tree->nodes[sib].age; age <= max_age; age++) {
                 Spr spr(sib, age, parent, real_spr->coal_time);
                 recomb_sum += exp(calc_spr_prob(model, tree, spr, lineages,
                                                 treelen, num_coal, num_nocoal,
                                                 0, true));
              }
            }
            lnl += log(pr_recomb * recomb_sum +
                       (possible_no_recomb ? pr_no_recomb : 0.0));
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
    /*
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    return calc_arg_likelihood_parsimony(&model, &sequences, trees);
    */
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


double arghmm_tree_prior_prob(LocalTrees *trees,
                              double *times, int ntimes, double *popsizes)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, 0.0, 0.0);
    LineageCounts lineages(ntimes);

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


} // namespace argweaver

