//=============================================================================
// transitions

#include "matrices.h"
#include "total_prob.h"
#include "thread.h"
#include "trans.h"


namespace argweaver {

void TransMatrix::initialize(const ArgModel *model, int nstates)
{
    ntimes = model->ntimes;
    npaths = model->num_pop_paths();
    smc_prime = model->smc_prime;
    int data_len=0;
    if (smc_prime) {
        data_len = npaths * ntimes + 2 * ntimes;
    } else {
        data_len = ntimes*2 + npaths*ntimes*4 + npaths*npaths*ntimes*3;
    }
    data_alloc = new double [data_len];

    C1_prime = new MultiArray(3, npaths, npaths, 2*ntimes);
    Q1_prime = new MultiArray(3, npaths, npaths, 2*ntimes);
    if (smc_prime) {
        B0_prime = new MultiArray(2, npaths, ntimes);
        B1_prime = new MultiArray(3, npaths, npaths, ntimes);
        B2_prime = new MultiArray(3, npaths, npaths, ntimes);
        E0_prime = new MultiArray(2, npaths, ntimes);
        E1_prime = new MultiArray(3, npaths, npaths, ntimes);
        E2_prime = new MultiArray(3, npaths, npaths, ntimes);
        F0_prime = new MultiArray(2, npaths, ntimes);
        F1_prime = new MultiArray(3, npaths, npaths, ntimes);
        F2_prime = new MultiArray(3, npaths, npaths, ntimes);
        C0_prime = new MultiArray(2, npaths, 2*ntimes);
        Q0_prime = new MultiArray(2, npaths, 2*ntimes);
        G0_prime = new MultiArray(2, npaths, ntimes);
        G1_prime = new MultiArray(3, npaths, npaths, ntimes);
        G2_prime = new MultiArray(3, npaths, npaths, ntimes);
        self_recomb = new double[nstates];
    } else {
        E = new double* [npaths];
        lnB = new double** [npaths];
        lnE2 = new double** [npaths];
        lnNegG1 = new double** [npaths];
        G2 = new double* [npaths];
        G3 = new double* [npaths];
    }
    path_prob = new double* [npaths];
    D = &data_alloc[0];
    norecombs = &data_alloc[ntimes];
    int idx=ntimes*2;
    for (int i=0; i < npaths; i++) {
        if (!smc_prime) {
            E[i] = &data_alloc[idx]; idx += ntimes;
            G2[i] = &data_alloc[idx]; idx += ntimes;
            G3[i] = &data_alloc[idx]; idx += ntimes;
            lnB[i] = new double* [npaths];
            lnE2[i] = new double* [npaths];
            lnNegG1[i] = new double* [npaths];
            for (int j=0; j < npaths; j++) {
                lnB[i][j] = &data_alloc[idx]; idx += ntimes;
                lnE2[i][j] = &data_alloc[idx]; idx += ntimes;
                lnNegG1[i][j] = &data_alloc[idx]; idx += ntimes;
            }
        }
        path_prob[i] = &data_alloc[idx]; idx += ntimes;
    }
    assert(idx == data_len);

    if (npaths > 1) {
        paths_equal = new bool*** [npaths];
        for (int i=0; i < npaths; i++) {
            paths_equal[i] = new bool** [npaths];
            for (int j=0; j < npaths; j++) {
                paths_equal[i][j] = new bool *[ntimes];
                for (int k=0; k < ntimes; k++) {
                    paths_equal[i][j][k] = new bool[ntimes];
                    for (int l=k; l < ntimes; l++) {
                        paths_equal[i][j][k][l] =
                            model->paths_equal(i,j,k,l);
                    }
                }
            }
        }
        paths_equal_max = new int** [npaths];
        for (int i=0; i < npaths; i++) {
            paths_equal_max[i] = new int* [npaths];
            for (int j=0; j < npaths; j++) {
                paths_equal_max[i][j] = new int[ntimes];
                for (int k=0; k < ntimes; k++) {
                    paths_equal_max[i][j][k] = -1;
                    for (int l=k; l < ntimes; l++) {
                        if (paths_equal[i][j][k][l])
                            paths_equal_max[i][j][k] = l;
                        else break;
                    }
                }
            }
        }
    } else {
        paths_equal = NULL;
        paths_equal_max = NULL;
    }
}

void calc_coal_rates_partial_tree(const ArgModel *model, const LocalTree *tree,
                                  const LineageCounts *lineages,
                                  MultiArray *coal_rates,
                                  MultiArray *coal_rates_noprime,
                                  bool *do_path, int minage)
{
    if (model->smc_prime) {
        for (int path=0; path < model->num_pop_paths(); path++) {
            if (do_path != NULL && do_path[path] == false) continue;
            for (int i=0; i < 2*model->ntimes - 1; i++) {
                int pop = model->get_pop(path, (i+1)/2);
                int nbranch = lineages->nbranches_pop[pop][i];
                coal_rates_noprime->set(model->coal_time_steps[i] *
                                        nbranch / (2.0 * model->popsizes[pop][i]),
                                        path, i);
                coal_rates_noprime->set(1.0, path, 2 * model->ntimes - 1);
            }
            coal_rates_noprime->set(1.0, path, 2 * model->ntimes - 1);
            for (int path2=0; path2 < model->num_pop_paths(); path2++) {
                if (do_path != NULL && do_path[path2] == false) continue;
                for (int i=0; i < 2*model->ntimes - 1; i++) {
                    int pop = model->get_pop(path, (i+1)/2);
                    int pop2 = model->get_pop(path2, (i+1)/2);
                    int nbranch = lineages->nbranches_pop[pop][i];
                    if (pop == pop2) nbranch++;
                    coal_rates->set(model->coal_time_steps[i] *
                                    nbranch / (2.0 * model->popsizes[pop][i]),
                                    path, path2, i);
                }
                coal_rates->set(1.0, path, path2, 2*model->ntimes - 1);
            }
        }
    } else {
        for (int path=0; path < model->num_pop_paths(); path++) {
            if (do_path != NULL && do_path[path] == false) continue;
            for (int path2=0; path2 < model->num_pop_paths(); path2++) {
                if (do_path != NULL && do_path[path2] == false) continue;
                for (int i=0; i<2*model->ntimes-1; i++) {
                    int pop = model->get_pop(path, (i+1)/2);
                    int pop2 = model->get_pop(path2, (i+1)/2);
                    int nbranch = lineages->nbranches_pop[pop][i];
                    if (pop != pop2) nbranch--;
                    else if (i < minage*2) nbranch--;  //same pop
                    if (nbranch < 0) nbranch = 0;
                    //                assert(nbranch >= 0);
                    coal_rates->set(model->coal_time_steps[i] *
                                    nbranch / (2.0 * model->popsizes[pop][i]),
                                    path, path2, i);
                }
                coal_rates->set(1.0, path, path2, 2 * model->ntimes - 1);
            }
        }
    }
}

// self recombination probability
// recomb along a branch that starts at time "minage"
// and coalesces at time d (d <= end time of branch)
// given that the branch being threaded coalesces at time a
double TransMatrix::self_recomb_prob(int a, int path_a,
                                     int minage, int d, int path_d) {
    double prob = 0.0;
    if (minage < d) {
        double pathprob = path_prob[path_d][d];
        double cterm;
        if (d <= a)
            cterm = exp(-C1_prime->get(path_d, path_a, 2*d-2));
        else
            cterm = exp(-C0_prime->get(path_d, 2*d-2)
                        +C0_prime->get(path_d, 2*a-1)
                        -C1_prime->get(path_d, path_a, 2*a-1));
        double bterm;
        if (d-1 > a) {
            bterm = B0_prime->get(path_d, d-1)
                - B0_prime->get(path_d, a)
                + B1_prime->get(path_d, path_a, a)
                + B2_prime->get(path_d, path_a, a-1);
        } else if (d-1 == a) {
            bterm = B1_prime->get(path_d, path_a, a)
                + B2_prime->get(path_d, path_a, a-1);
        } else {
            bterm = B2_prime->get(path_d, path_a, d-1);
        }
        if (minage > 0) {
            if (minage-1 > a) {
                bterm  -= ( B0_prime->get(path_d, minage-1)
                            - B0_prime->get(path_d, a)
                            + B1_prime->get(path_d, path_a, a)
                            + B2_prime->get(path_d, path_a, a-1));
            } else if (minage-1 == a) {
                bterm -= ( B1_prime->get(path_d, path_a,  a)
                           + B2_prime->get(path_d, path_a, a-1));
            } else {
                bterm -= B2_prime->get(path_d, path_a, minage-1);
            }
        }
        double eterm;
        if (d > a) {
            eterm = E0_prime->get(path_d, d);
        } else if (d == a) {
            eterm = E1_prime->get(path_d, path_a, d);
        } else {
            eterm = E2_prime->get(path_d, path_a, d);
        }

        prob = pathprob * D[a] * cterm * bterm * eterm;
        assert(!isnan(prob));
    }
    double term2;
    if (d > a) {
        term2 = G0_prime->get(path_d, d) * F0_prime->get(path_d, d);
    } else if (d == a) {
        term2 = G1_prime->get(path_d, path_a, d)
            * F1_prime->get(path_d, path_a, d);
    } else {
        term2 = G2_prime->get(path_d, path_a, d)
            * F2_prime->get(path_d, path_a, d);
    }
    prob += D[a] * term2;
    assert(!isnan(prob));
    return prob;
}



void TransMatrix::assert_transmat(const LocalTree *tree,
                                  const ArgModel *model,
                                  const States &states,
                                  bool internal, int minage) const {
    /*    const double *times = model->times;
    int root_age_index, root_age, minage;
    double subtree_age = times[tree->nodes[subtree_root].age];
    double treelen0;
    int minage;
    if (internal) {
        root_age_index = tree->nodes[maintree_root].age;
        root_age = times[root_age_index];
        minage = max(minage, tree->nodes[subtree_root].age);
        subtree_age = times[tree->nodes[subtree_root].age];
        treelen0 = get_treelen_internal(tree, times, ntimes);
    } else {
        root_age_index = tree->nodes[tree->root].age;
        root_age = times[root_age_index];
        treelen0 = get_treelen(tree, times, ntimes, false);
        minage = minage0;
        } */
    for (unsigned int s1=0; s1 < states.size(); s1++) {
        const State state1 = states[s1];
        int a = state1.time;
        /*        double treelen = treelen0 + times[a] - minage;
        if (a > root_age_index)
        treelen += (a - times[root_age_index]);
        double pr_recomb = 1.0 - exp(-rho * treelen);*/
        double total_prob=0.0;
        //        double total_prob_recomb=0.0;
        for (unsigned int s2=0; s2 < states.size(); s2++) {
            const State state2 = states[s2];
            double p1=0, p2=0;
            p1 = get_time(a, state2.time, -1,
                          state1.pop_path,
                          state2.pop_path,
                          -1,
                          minage, false, -1);
            if (state2.node == state1.node) {
                p2 = get_time(a, state2.time, tree->nodes[state2.node].age,
                              state1.pop_path, state2.pop_path,
                              tree->nodes[state2.node].pop_path,
                              minage, true, s1) - p1;
            }
            total_prob += p1 + p2;
        }
        //        printf("total_prob %i = %e\n", s1, total_prob);
    }

}

// Calculate transition probability within a local block
void TransMatrix::calc_transition_probs_smcPrime(const LocalTree *tree,
                                                 const ArgModel *model,
                                                 const States &states,
                                                 const LineageCounts *lineages,
                                                 bool internal0, int minage0)
{
    assert(model->smc_prime);
    internal = internal0;
    minage = minage0;
    // get model parameters
    const int ntimes = model->ntimes;
    const double *times = model->times;
    const double rho = model->rho;
    int* const* ncoals_pop = lineages->ncoals_pop;
    int* const* nbranches_pop = lineages->nbranches_pop;
    const int num_paths = model->num_pop_paths();
    bool have_pop_path[num_paths];
    const int subtree_root = internal ? tree->nodes[tree->root].child[0] : -1;
    const int maintree_root = internal ? tree->nodes[tree->root].child[1] : -1;
    if (num_paths > 1) {
        for (int i=0; i < num_paths; i++)
            have_pop_path[i]=false;
        for (unsigned int i=0; i < states.size(); i++)
            have_pop_path[states[i].pop_path] = true;
        for (int i=0; i < tree->nnodes; i++) {
            if (i != subtree_root)
                have_pop_path[tree->nodes[i].pop_path] = true;
        }
    } else have_pop_path[0] = true;

    // determine tree information: root, root age, tree length
    int root_age_index;
    double root_age;
    double treelen;
    if (internal) {
        const double subtree_age = times[tree->nodes[subtree_root].age];
        root_age_index = tree->nodes[maintree_root].age;
        root_age = times[root_age_index];
        // NOTE: subtree_age is discounted in advance to offset +times[b]
        treelen = get_treelen_internal(tree, times, ntimes) - subtree_age;
        minage = max(minage, tree->nodes[subtree_root].age);
    } else {
        root_age_index = tree->nodes[tree->root].age;
        root_age = times[root_age_index];
        treelen = get_treelen(tree, times, ntimes, false);
    }

    calc_coal_rates_partial_tree(model, tree, lineages,
                                 Q1_prime, Q0_prime,
                                 have_pop_path, minage);

    // compute cumulative coalescent rates
    for (int i=0; i < num_paths; i++) {
        if (have_pop_path[i]) {
            for (int j=0; j < num_paths; j++) {
                if (have_pop_path[j]) {
                    for (int b=0; b < 2*ntimes - 1; b++)
                        C1_prime->set(C1_prime->get(i, j, b-1)
                                      + Q1_prime->get(i, j, b),
                                      i, j, b);
                }
            }
        }
    }
    for (int i=0; i < num_paths; i++) {
        if (have_pop_path[i]) {
            for (int j=0; j < 2*ntimes - 1; j++) {
                C0_prime->set(C0_prime->get(i, j-1) + Q0_prime->get(i, j),
                              i, j);
            }
        }
    }

    // calculate transition matrix terms
    for (int b=0; b<ntimes-1; b++) {
        // get tree length
        double treelen2 = treelen + times[b];
        if (b > root_age_index){
            // add wrapped branch
            treelen2 += times[b] - root_age;
        }
        norecombs[b] = exp(-max(rho * treelen2, rho));
        D[b] = ( treelen2 == 0 ? 0.0 : (1.0 - exp( -rho * treelen2)) / treelen2 );
    }

    for (int path=0; path < num_paths; path++) {
        if (! have_pop_path[path]) continue;
        for (int b=0; b < ntimes-1; b++) {
            double curr_path_prob = model->path_prob(path, 0, b);
            path_prob[path][b] = curr_path_prob;
            int pop = model->get_pop(path, b);
            int ncoal = ncoals_pop[pop][b];
            if (b >= ntimes - 2) {
                E0_prime->set(1.0/(double)ncoal, path, b);
                F0_prime->set(1.0/(double)ncoal, path, b);
            } else {
                E0_prime->set((1.0 - exp(-Q0_prime->get(path, 2*b)
                                         -Q0_prime->get(path, 2*b-1)))
                              / (double)ncoal, path, b);
                F0_prime->set((1.0 - exp(-Q0_prime->get(path, 2*b)))
                              / (double)ncoal, path, b);
            }

            int nrecomb = ncoals_pop[pop][b];
            int nbranch_above=0, nbranch_below=0;
            if (b < root_age_index) {
                if (b > 0) nbranch_below = nbranches_pop[pop][2*b-1];
                nbranch_above = nbranches_pop[pop][2*b];
            } else if (b == root_age_index) {
                if (b > 0) nbranch_below = nbranches_pop[pop][2*b-1];
                nbranch_above = 0;
                nrecomb--;
            } else  { // b > root_age_index
                nrecomb=1; // branchlens are zero; so recomb prob is zero;
                // set this to 1 to avoid by-zero division
            }
            double total_blen_b = (b == 0 ? 0.0 :
                                   model->coal_time_steps[2*b-1] * nbranch_below)
                + model->coal_time_steps[2*b] * nbranch_above;
            B0_prime->set(b == 0 ? 0 :
                          B0_prime->get(path, b-1), path, b);
            if (total_blen_b > 0.0)
                B0_prime->addVal(total_blen_b / (double)nrecomb *
                                 exp(C0_prime->get(path, 2*b-1))
                                 / curr_path_prob,
                                 path, b);
            G0_prime->set(total_blen_b == 0.0 ? 0 : total_blen_b / (double)nrecomb,
                          path, b);

            for (int path2=0; path2 < num_paths; path2++) {
                if (! have_pop_path[path2]) continue;
                int ncoal1 = ncoal;
                int ncoal2 = ncoal;
                int pop2 = model->get_pop(path2, b);
                if (pop2 == pop) {
                    ncoal1 += 2;
                    ncoal2 += 1;
                }
                if (b >= ntimes - 2) {
                    E1_prime->set(1.0/(double)ncoal1, path, path2, b);
                    E2_prime->set(1.0/(double)ncoal2, path, path2, b);
                    F1_prime->set(1.0/(double)ncoal1, path, path2, b);
                    F2_prime->set(1.0/(double)ncoal2, path, path2, b);
                } else {
                    E1_prime->set(( 1.0 - exp(-Q1_prime->get(path, path2, 2*b-1)
                                              -Q0_prime->get(path, 2*b)))
                                  / ((double)ncoal1), path, path2, b);
                    E2_prime->set(
                                  ( 1.0 - exp(-Q1_prime->get(path, path2, 2*b-1)
                                              -Q1_prime->get(path, path2, 2*b)))
                                  / ((double)ncoal2), path, path2, b);
                    F1_prime->set(( 1.0 - exp(-Q0_prime->get(path, 2*b) ))
                                  / ((double)ncoal1), path, path2, b);
                    F2_prime->set(
                                  ( 1.0 - exp(-Q1_prime->get(path, path2, 2*b)))
                                  / ((double)ncoal2), path, path2, b);
                }
                int nrecomb1 = ncoals_pop[pop][b];
                int nrecomb2 = nrecomb1;
                int nbranch_below1, nbranch_below2;
                int nbranch_above1, nbranch_above2;
                double blen1, blen2;
                nbranch_below1 = nbranch_below2 =
                    ( b == 0 ? 0 : nbranches_pop[pop][2*b-1] );
                nbranch_above1 = nbranch_above2 = nbranches_pop[pop][2*b];
                if (pop == pop2) {
                    if (b > minage) {
                        nbranch_below1++;
                        nbranch_below2++;
                    }
                    if (b >= minage) {
                        nbranch_above2++;
                        nrecomb1 += 2;
                        nrecomb2++;
                    }
                }
                if (b >= root_age_index) {
                    nrecomb1--;
                    nbranch_above1=0;
                }
                blen1 = ( b == 0 ? 0 :
                          model->coal_time_steps[2*b-1]*(double)nbranch_below1 )
                    + model->coal_time_steps[2*b] * (double)nbranch_above1;
                blen2 = ( b == 0 ? 0 :
                          model->coal_time_steps[2*b-1]*(double)nbranch_below2 )
                    + model->coal_time_steps[2*b] * (double)nbranch_above2;

                // B1 is not a sum since it applies to case where previous
                // branch coalesces at time b and not to consecutive intervals
                // of times
                B1_prime->set(blen1 == 0.0 ? 0 : blen1 / (double)nrecomb1
                              * exp(C1_prime->get(path, path2, 2*b-1))
                              / curr_path_prob,
                              path, path2, b);
                B2_prime->set(b == 0 ? 0 :
                              B2_prime->get(path, path2, b-1),
                              path, path2, b);
                if (blen2 > 0)
                    B2_prime->addVal( blen2 / (double)nrecomb2
                                      * exp(C1_prime->get(path, path2, 2*b-1))
                                      / curr_path_prob,
                                      path, path2, b);
                G1_prime->set(blen1 == 0.0 ? 0 : blen1 / (double)nrecomb1,
                              path, path2, b);
                G2_prime->set(blen2 == 0.0 ? 0 : blen2 / (double)nrecomb2,
                              path, path2, b);
            }
        }
    }
    // now need to calculate self_recombs[node][time] for all nodes, times
    // not implemented efficiently yet
    for (unsigned int s=0; s < states.size(); s++) {
        int path_a = states[s].pop_path;
        int a = states[s].time;
        int node = states[s].node;
        double prob=0.0;
        // first consider new branch which goes from minage to a
        for (int d=minage; d <= a; d++) {
            prob += self_recomb_prob(a, path_a, minage, d,
                                     path_a);

        }
        if (internal && node == maintree_root) {
            assert(a >= root_age_index);
            for (int d=root_age_index; d <= a; d++)
                prob += self_recomb_prob(a, path_a, root_age_index, d,
                                         tree->nodes[maintree_root].pop_path);
        } else if ((!internal) && node == tree->root) {
            assert(a >= root_age_index);
            for (int d = root_age_index; d <= a; d++)
                prob += self_recomb_prob(a, path_a, root_age_index, d,
                                         tree->nodes[tree->root].pop_path);
        }
        for (int i=0; i < tree->nnodes; i++) {
            if (i == tree->root) continue;
            if (internal && i == subtree_root) continue;
            if (internal && i == maintree_root) continue;
            int age = tree->nodes[i].age;
            int parent_age = tree->nodes[tree->nodes[i].parent].age;
            int path_d = tree->nodes[i].pop_path;
            if (i == node) {
                for (int d=age; d <= a; d++)
                    prob += self_recomb_prob(a, path_a, age, d, path_d);
                for (int d=a; d <=parent_age; d++) {
                    prob += self_recomb_prob(a, path_a, a, d, path_d);
                }
            } else {
                for (int d=age; d <= parent_age; d++)
                    prob += self_recomb_prob(a, path_a, age, d, path_d);
            }
        }
        self_recomb[s] = prob;
    }
    if (false) {
        assert_transmat(tree, model, states, lineages, minage0);
    }
}


// Calculate transition probability within a local block
void TransMatrix::calc_transition_probs(const LocalTree *tree,
                                        const ArgModel *model,
                                        const States &states,
                                        const LineageCounts *lineages,
                                        bool internal0, int minage0)
{
    if (model->smc_prime) {
        calc_transition_probs_smcPrime(tree, model, states, lineages,
                                       internal0, minage0);
        return;
    }
    assert(!model->smc_prime);
    internal = internal0;
    minage = minage0;
    // get model parameters
    const int ntimes = model->ntimes;
    const double *times = model->times;
    const double *time_steps = model->time_steps;
    const double rho = model->rho;
    const int* nbranches = lineages->nbranches;
    const int* nrecombs = lineages->nrecombs;
    int* const* ncoals_pop = lineages->ncoals_pop;
    const int num_paths = model->num_pop_paths();
    bool have_pop_path[num_paths];
    const int subtree_root = internal ? tree->nodes[tree->root].child[0] : -1;
    const int maintree_root = internal ? tree->nodes[tree->root].child[1] : -1;
    if (num_paths > 1) {
        for (int i=0; i < num_paths; i++)
            have_pop_path[i]=false;
        for (unsigned int i=0; i < states.size(); i++)
            have_pop_path[states[i].pop_path] = true;
        for (int i=0; i < tree->nnodes; i++) {
            if (i != subtree_root)
                have_pop_path[tree->nodes[i].pop_path] = true;
        }
    } else have_pop_path[0] = true;

    // determine tree information: root, root age, tree length
    int root_age_index;
    double root_age;
    double treelen;
    if (internal) {
        const double subtree_age = times[tree->nodes[subtree_root].age];
        root_age_index = tree->nodes[maintree_root].age;
        root_age = times[root_age_index];
        // NOTE: subtree_age is discounted in advance to offset +times[b]
        treelen = get_treelen_internal(tree, times, ntimes) - subtree_age;
        minage = max(minage, tree->nodes[subtree_root].age);
    } else {
        root_age_index = tree->nodes[tree->root].age;
        root_age = times[root_age_index];
        treelen = get_treelen(tree, times, ntimes, false);
    }

    calc_coal_rates_partial_tree(model, tree, lineages,
                                 Q1_prime, Q0_prime,
                                 have_pop_path, minage);

    // compute cumulative coalescent rates
    for (int i=0; i < num_paths; i++) {
        for (int j=0; j < num_paths; j++) {
            if (have_pop_path[i] && have_pop_path[j]) {
                for (int b=0; b < 2*ntimes - 1; b++)
                    C1_prime->set(C1_prime->get(i, j, b-1)
                                  + Q1_prime->get(i, j, b),
                                  i, j, b);
            }
        }
    }

    // calculate transition matrix terms
    for (int b=0; b<ntimes-1; b++) {
        // get tree length
        double treelen2 = treelen + times[b];
        double treelen2_b;
        if (b > root_age_index){
            // add wrapped branch
            treelen2 += times[b] - root_age;

            // add basal branch
            treelen2_b = treelen2 + time_steps[b];
        } else {
            // add basal branch
            treelen2_b = treelen2 + time_steps[root_age_index];
        }
        norecombs[b] = exp(-max(rho * treelen2, rho));
        D[b] = (1.0 - exp( -rho * treelen2)) / treelen2_b;
    }

    for (int path=0; path < num_paths; path++) {
        if (! have_pop_path[path]) continue;
        for (int b=0; b < ntimes-1; b++) {
            double curr_path_prob = model->path_prob(path, 0, b);
            path_prob[path][b] = curr_path_prob;
            int pop = model->get_pop(path, b);
            int ncoal = ncoals_pop[pop][b];
            int nrecomb = nrecombs[b];
            int nbranch = nbranches[b];
            if (b < minage) {
                ncoal--;
                nrecomb--;
                nbranch--;
            }
            for (int path2=0; path2 < num_paths; path2++) {
                if (! have_pop_path[path2]) continue;
                double term = C1_prime->get(path, path2, 2*b-1) + log(
                    time_steps[b] * (nbranch + 1.0) /
                    (nrecomb + 1.0) / curr_path_prob);
                if (b == 0)
                    lnB[path][path2][b] = term;
                else
                    lnB[path][path2][b] =
                        logadd(lnB[path][path2][b-1], term);
                lnE2[path][path2][b] = -C1_prime->get(path, path2, 2*b-2) +
                    (b < ntimes - 2 ?
                     log(1 - exp(-Q1_prime->get(path, path2, 2*b)
                                 -Q1_prime->get(path, path2, 2*b-1))) : 0.0);
                assert(!isnan(exp(2.0*lnE2[path][path2][b])));
                lnNegG1[path][path2][b] =
                    C1_prime->get(path, path2, 2*b-1) +
                    log( - time_steps[b] / curr_path_prob * (
                        (nbranch / (nrecomb + 1.0 + int(b < root_age_index)))
                        - (nbranch + 1.0) / (nrecomb + 1.0)));
            }
            G2[path][b] = (b<ntimes-2 ?
                           1.0 - exp(-Q1_prime->get(path, path, 2*b))
                           : 1.0) *
                time_steps[b] *
                (nbranch + 1.0) / (nrecomb + 1.0) / curr_path_prob;
            G3[path][b] = (b<ntimes-2 ?
                           1.0 - exp(-Q1_prime->get(path, path, 2*b))
                           : 1.0) *
                time_steps[b] *
                (nbranch / (nrecomb + 1.0 + int(b < root_age_index)) /
                 curr_path_prob);
            E[path][b] = 1.0 / ncoal;
        }
    }
    if (false) {
        assert_transmat(tree, model, states, lineages, minage0);
    }
}

double TransMatrix::get_time(int a, int b, int c,
                             int path_a, int path_b, int path_c,
                             int minage, bool same_node, int state_a) const
{
    if (a < minage || b < minage)
        return 0.0;

    const int p = ( npaths == 1 ? ntimes :
                    paths_equal_max[path_a][path_b][minage] );
    if (p == -1) return 0.0;

    if (!smc_prime) {
        double prob;
        double term1 = D[a] * E[path_b][b] * path_prob[path_b][b];
        double minage_term = 0.0;
        if (minage > 0) {
            minage_term = exp(lnE2[path_b][path_b][b] +
                              lnB[path_b][path_b][minage-1]);
        }
        if (p < a && p < b) {
            prob = term1 * (exp(lnE2[path_b][path_b][b] +
                                lnB[path_b][path_b][p])
                            - minage_term);
        } else if (a <= p && a < b) {
            prob = term1 * (exp(lnE2[path_b][path_b][b] +
                                lnB[path_b][path_b][a]) -
                            exp(lnE2[path_b][path_b][b] +
                                lnNegG1[path_b][path_b][a])
                            - minage_term);
        } else if (a == b) {
            prob = term1 * ((b > 0 ? exp(lnE2[path_b][path_b][b] +
                                         lnB[path_b][path_b][b-1]) : 0.0) +
                            G3[path_b][b] - minage_term);
        } else { // b < a
            prob = term1 * ((b > 0 ? exp(lnE2[path_b][path_b][b] +
                                         lnB[path_b][path_b][b-1]) : 0.0)
                            + G2[path_b][b] - minage_term);
        }
        if (isnan(prob)) {
            assert(false);
        }
        if (! same_node) return prob;

        // now add same_node term
        // norecomb case
        if (a == b &&
            (npaths == 1 || paths_equal[path_a][path_b][minage][a]))
            prob += norecombs[a];

        if (npaths > 1 && a < b &&
            !(paths_equal[path_b][path_c][a][b] &&
              paths_equal[path_a][path_b][minage][a])) {
            if (isnan(prob))
                assert(false);
            return prob;
        }
        if (npaths > 1 && b <= a &&
            !(paths_equal[path_a][path_c][b][a] &&
              paths_equal[path_b][path_a][minage][b]))
            return prob;

        term1 = D[a] * E[path_c][b] * path_prob[path_c][b];
        minage_term = 0.0;
        if (c > 0)
            minage_term = exp(lnE2[path_c][path_b][b]
                              + lnB[path_c][path_b][c-1]);

        if (a < b) {
            prob += term1 * (exp(lnE2[path_c][path_b][b] +
                                 lnB[path_c][path_b][a]) -
                             exp(lnE2[path_c][path_b][b] +
                                 lnNegG1[path_c][path_b][a])
                             - minage_term);
        } else if (a == b) {
            prob += term1 * ((b > 0 ? exp(lnE2[path_c][path_b][b] +
                                          lnB[path_c][path_b][b-1]) : 0.0) +
                             G3[path_c][b] - minage_term);
        } else { // b < a
            prob += term1 * ((b > 0 ? exp(lnE2[path_c][path_b][b] +
                                          lnB[path_c][path_b][b-1]) : 0.0)
                             + G2[path_c][b] - minage_term);
        }
        if (isnan(prob))
            assert(0);
        return prob;
    } else {  //smc prime calculations
        double term1=0, term2=0, minage_term=0;
        term1 = D[a]*path_prob[path_a][b];
        if (a < b) {
            term1 *= E0_prime->get(path_b, b)
                * exp(-C0_prime->get(path_b, 2*b-2)
                      + C0_prime->get(path_b, 2*a-1)
                      - C1_prime->get(path_b, path_a, 2*a-1));
        } else if (a == b) {
            term1 *= E1_prime->get(path_b, path_a, b)
                * exp(-C1_prime->get(path_b, path_a, 2*b-2));
        } else {
            term1 *= E2_prime->get(path_b, path_a, b)
                * exp(-C1_prime->get(path_b, path_a, 2*b-2));
        }
        int k_max = b-1;
        if (a < k_max) k_max = a;
        if (p < k_max) k_max = p;
        double b_term;
        if (k_max < a) {
            b_term = B2_prime->get(path_b, path_a, k_max);
        } else {
            assert(k_max == a);
            b_term = B2_prime->get(path_b, path_a, a-1) +
                B1_prime->get(path_b, path_a, a);
        }
        //            term1 *= b_term;
        if (b == a && b <= p) {
            term2 = D[a] * G1_prime->get(path_b, path_a, b)
                * F1_prime->get(path_b, path_a, b);
        } else if (b < a && b <= p) {
            term2 = D[a] * G2_prime->get(path_b, path_a, b)
                * F2_prime->get(path_b, path_a, b);
        }

        if (minage > 0) {
            assert(minage <= b);
            assert(minage <= a);
            assert(minage <= p);
            minage_term = B2_prime->get(path_b, path_a, minage-1);
        }

        double prob = term1 * (b_term - minage_term) + term2;
        assert(prob >=0 && prob <= 1);

        if (!same_node) return prob;  // must be recombination on threaded branch

        if (npaths > 1 && a < b &&
            !(paths_equal[path_b][path_c][a][b] &&
              paths_equal[path_a][path_b][minage][a])) {
            if (isnan(prob))
                assert(false);
            return prob;
        }
        if (npaths > 1 && b <= a &&
            !(paths_equal[path_a][path_c][b][a] &&
              paths_equal[path_b][path_a][minage][b]))
            return prob;

        assert(state_a >= 0);
        // now add same_node term
        // no recomb or self_recomb case
        if (a < b) {
            prob += D[a] * path_prob[path_c][b]
                * E0_prime->get(path_c, b)
                * exp(-C0_prime->get(path_c, 2*b-2)
                      +C0_prime->get(path_c, 2*a-1)
                      -C1_prime->get(path_c, path_a, 2*a-1))
                * (B2_prime->get(path_c, path_a, a-1)
                   + B1_prime->get(path_c, path_a, a)
                   - B2_prime->get(path_c, path_a, c-1));
        } else if (a == b) {
            prob *= 2.0;  // because could coal to parent or sister branch
            prob += norecombs[a] + self_recomb[state_a];

            prob += 2.0 * (( D[a] * path_prob[path_c][b]
                             * E1_prime->get(path_c, path_a, b)
                             * exp(-C1_prime->get(path_c, path_a, 2*b-2))
                             * ( B2_prime->get(path_c, path_a, b-1)
                                 - B2_prime->get(path_c, path_a, c-1)))
                           + ( D[a] * G1_prime->get(path_c, path_a, b) *
                               F1_prime->get(path_c, path_a, b)));
        } else if (a > b) {
            prob += D[a] * path_prob[path_c][b]
                * E2_prime->get(path_c, path_a, b)
                * exp(-C1_prime->get(path_c, path_a, 2*b-2))
                * ( B2_prime->get(path_c, path_a, b-1)
                    - B2_prime->get(path_c, path_a, c-1))
                + D[a] * G2_prime->get(path_c, path_a, b)
                * F2_prime->get(path_c, path_a, b);
        }
        if (isnan(prob))
            assert(0);
        assert(prob >=0 && prob <= 1);
        return prob;
    }
}

// Copies transition probability to dense probability matrix
void get_transition_probs(const LocalTree *tree, const ArgModel *model,
                          const States &states, const LineageCounts *lineages,
                          const TransMatrix *matrix, double **transprob)
{
    // get tree and model information
    const int nstates = states.size();

    // calculate full state transition matrix
    for (int i=0; i<nstates; i++)
        for (int j=0; j<nstates; j++)
            transprob[i][j] = matrix->get_log(tree, states, i, j);
}


// Calculate transition probability within a local block
// Stores transition probabilities in dense matrix
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
                           const States &states, const LineageCounts *lineages,
                           double **transprob)
{
    TransMatrix matrix(model, states.size());
    matrix.calc_transition_probs(tree, model, states, lineages);
    get_transition_probs(tree, model, states, lineages, &matrix, transprob);
}


//=============================================================================
// functions for switch matrix calculation


// Returns the deterministic transitions that occur when switching between
// blocks.  Transitions are stored in the array 'next_states' such that
//   next_states[i] = j
// for deterministic transition i-->j.
// If state i has no deterministic transition
//   next_states[i] = -1.
void get_deterministic_transitions(
    const LocalTree *last_tree, const LocalTree *tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    NodeStateLookup &state2_lookup,
    const ArgModel *model, int *next_states, bool internal)
{
    // recomb_node in tree and last_tree
    // coal_node in last_tree
#ifdef DEBUG
    static int count=0;
    count++;
#endif

    //const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int nstates1 = states1.size();

    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i].node;
        const int time1 = states1[i].time;
        const int path1 = states1[i].pop_path;

        if ((node1 == spr.coal_node && time1 == spr.coal_time) ||
            (node1 == spr.recomb_node && time1 == spr.recomb_time)) {
            // not a deterministic case
            next_states[i] = -1;
        } else if (node1 == spr.coal_node && spr.coal_node == spr.recomb_node) {
            // "bubble" recomb; tree topology does not change but pop path of
            // recomb branch may
            if (time1 < spr.recomb_time || time1 > spr.coal_time) {
                next_states[i] = state2_lookup.lookup(mapping[node1], time1, path1);
            } else {
                if (model->pop_tree != NULL) {
                    // this should be true because state paths are meant to be consistent with node that
                    // they coalesce to
                    assert(model->paths_equal(path1, last_tree->nodes[node1].pop_path, time1, spr.coal_time));
                }
                next_states[i] = state2_lookup.lookup(mapping[node1], spr.coal_time, path1);
            }
        } else if (node1 != spr.recomb_node && spr.recomb_node == spr.coal_node) {
            next_states[i] = state2_lookup.lookup(mapping[node1], time1, path1);
        } else if (node1 != spr.recomb_node) {
            // SPR only removes a subset of descendents, if any
            // trace up from remaining leaf to find correct new state

            int node2;
            const LocalNode *node = &last_nodes[node1];
            bool disrupt = false;


            if (node->child[0] == -1) {
                // SPR can't disrupt leaf branch
                node2 = node1;

            } else {
                const int child1 = node->child[0];
                const int child2 = node->child[1];

                if (spr.recomb_node == child1 && spr.coal_node != spr.recomb_node) {
                    // right child is not disrupted
                    node2 = mapping[child2];
                    disrupt = true;
                } else if (spr.recomb_node == child2 && spr.coal_node != spr.recomb_node) {
                    // left child is not disrupted
                    node2 = mapping[child1];
                    disrupt = true;
                } else {
                    // node is not disrupted
                    node2 = mapping[node1];
                }
            }

            // optionally walk up
            if ((spr.coal_node == node1 && spr.coal_time < time1) ||
                (mapping[spr.coal_node] == node2 && spr.coal_time < time1) ||
                (disrupt && mapping[spr.coal_node] == node2 &&
                 spr.coal_time <= time1))
            {
                // coal occurs under us
                node2 = nodes[node2].parent;
            }

            if (internal && nodes[node2].age > time1) {
                // this is a deadend
                next_states[i] = -1;
                continue;
            }

            assert(nodes[node2].age <= time1);
            const int p = nodes[node2].parent;
            if (p != -1) {
                if (internal && time1 > nodes[p].age) {
                    // this is a deadend
                    next_states[i] = -1;
                    continue;
                }

                assert(time1 <= nodes[p].age);
            }
            int path2 = path1;
            if (internal && model->pop_tree != NULL) {
                int last_subtree_root = last_tree->nodes[last_tree->root].child[0];
                if (spr.coal_node == last_subtree_root) {
                    /* just check that the populations meet at new node*/
                    if (model->get_pop(spr.pop_path, spr.coal_time) !=
                        model->get_pop(path1, spr.coal_time)) {
                        next_states[i] = -1;
                        continue;
                    }
                    int subtree_root = tree->nodes[tree->root].child[0];
                    if (mapping[last_subtree_root] != subtree_root && mapping[last_subtree_root] != -1) {
                        int tmp = mapping[last_subtree_root];
                        if (!model->pop_tree->paths_equal(tree->nodes[tmp].pop_path,
                                                          path1, tree->nodes[tmp].age,
                                                          tree->nodes[tmp].parent == -1 ? -1 :
                                                          tree->nodes[tree->nodes[tmp].parent].age)) {
                            next_states[i] = -1;
                            continue;
                        }
                    } else if (mapping[last_subtree_root] == -1) {
                        // this means recomb is below subtree root
                        // first assert that last_subtree_root and subtree_root have the same children
                        assert((mapping[last_nodes[last_subtree_root].child[0]] ==
                                nodes[subtree_root].child[0] &&
                                mapping[last_nodes[last_subtree_root].child[1]] ==
                                nodes[subtree_root].child[1]) ||
                               (mapping[last_nodes[last_subtree_root].child[0]] ==
                                nodes[subtree_root].child[1] &&
                                mapping[last_nodes[last_subtree_root].child[1]] ==
                                nodes[subtree_root].child[0]));
                        // now check that paths of unbroken child are consistent
                        int other;
                        if (spr.recomb_node == last_nodes[last_subtree_root].child[0])
                            other = last_nodes[last_subtree_root].child[1];
                        else {
                            assert(spr.recomb_node == last_nodes[last_subtree_root].child[1]);
                            other = last_nodes[last_subtree_root].child[0];
                        }
                        int other_path = model->consistent_path(last_nodes[other].pop_path,
                                                                path1,
                                                                last_nodes[other].age,
                                                                last_nodes[last_subtree_root].age,
                                                                spr.coal_time);
                        assert(spr.coal_time == nodes[subtree_root].age);
                        assert(mapping[other] != -1);
                        assert(last_nodes[other].age == nodes[mapping[other]].age);
                        if (!model->paths_equal(nodes[mapping[other]].pop_path,
                                                other_path, last_nodes[other].age,
                                                nodes[subtree_root].age)) {
                            next_states[i] = -1;
                            continue;
                        }
                    }

                } else {
                    int c[2];
                    for (int j=0; j <= 1; j++)
                        c[j] = last_tree->nodes[last_subtree_root].child[j];
                    for (int c1=0; c1 <= 1; c1++) {
                        int c2 = (c1 == 0 ? 1 : 0);
                        if (c[c1] == spr.recomb_node) {
                            if (c[c2] == spr.coal_node) {
                                path2 =
                                    model->consistent_path(
                                        last_tree->nodes[c[c2]].pop_path,
                                        path1,
                                        spr.coal_time,
                                        last_tree->nodes[last_subtree_root].age,
                                        time1);
                            } else {
                                path2 = model->consistent_path(
                                    last_tree->nodes[c[c2]].pop_path, path1,
                                    last_tree->nodes[c[c2]].age,
                                    last_tree->nodes[last_subtree_root].age,
                                    time1);
                            }
                            break;
                        }
                    }
                }
            }

            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1, path2);
        } else {   // node1 == recomb_node but time1 != recomb_time
            // SPR is on same branch as new chromosome
            if (spr.recomb_time > time1) {
                // we move with SPR subtree
                next_states[i] =
                    state2_lookup.lookup(mapping[spr.recomb_node], time1, path1);
                //                else next_states[i] = -1;
            } else {  // new coal is on broken branch


                // SPR subtree moves out from underneath us
                // therefore the new chromosome coalesces with
                // the branch above the subtree

                // search up for parent
                const int parent = last_nodes[spr.recomb_node].parent;
                const int time2 = last_nodes[parent].age;

                // find other child
                const int *c = last_nodes[parent].child;
                const int other = (c[1] == spr.recomb_node ? c[0] : c[1]);

                // find new state in tree
                const int node2 = (other == spr.coal_node ?
                    nodes[mapping[other]].parent : mapping[other]);

                //bubbles on coal node were dealt with previously
                assert(spr.coal_node != node1);
                if (!internal) {
                    next_states[i] = state2_lookup.lookup(node2, time2, path1);
                    continue;
                }

                const int subtree_root = tree->nodes[tree->root].child[0];
                const int last_subtree_root = last_tree->nodes[last_tree->root].child[0];
                const int minage = tree->nodes[subtree_root].age;
                const int last_minage = last_tree->nodes[last_subtree_root].age;
                                                          // state1 is on recomb node; spr is not bubble

                if (minage > time2) {
                    next_states[i] = -1;
                    continue;
                }
                int path2;
                if (spr.coal_node == last_subtree_root) {
                    int node = mapping[last_subtree_root];
                    assert(node != -1);
                    int nodeparent = nodes[node].parent;
                    assert(nodeparent != -1);
                    if (! model->paths_equal(path1, nodes[node].pop_path,
                                             last_nodes[last_subtree_root].age,
                                             min(time1, spr.coal_time))) {
                        next_states[i] = -1;
                        continue;
                    }
                    if (nodeparent == subtree_root) {
                        if (time1 >= spr.coal_time) {
                            path2 = model->path_to_root(last_nodes, node1, time1);
                            path2 = model->consistent_path(path1, path2, spr.coal_time, time1, -1);
                        } else if (time1 < spr.coal_time) {
                            assert(minage == spr.coal_time);
                            path2 = path1;
                        }
                    } else {
                        assert(mapping[last_subtree_root] == subtree_root);
                        assert(minage == last_minage);
                        path2 = model->path_to_root(last_nodes, node1, time1);
                        path2 = model->consistent_path(path1, path2, minage, time1, -1);
                    }
                    //                    else next_states[i] = state2_lookup.lookup(node2, time2, nodes[node].pop_path);
                } else {
                    assert(minage == last_minage);
                    assert(mapping[last_subtree_root] == subtree_root);
                    path2 = path1;
                }
                next_states[i] = state2_lookup.lookup(node2, time2, path2);
            }
        }
    }
}


// NOTE: fixed for pop_path, but:
// on return, next_states[0] contains state node/time if coal at recombining
// node and go along with node to new coalescence
// next_states[1] contains node/time if escape recombining node and coalesce
// with parent instead.
// HOWEVER, cannot know population path- not enough information here-
// must be determined later
// NOTE ALSO: states1 and states2 are not used??
// this is called to compute the possible next_states in the case that previous
// state is same node/time as the existing recombination
// on exit, next_states[0] holds the next state in the "stay" case (coal with
//   recomb branch)
// and next_states[1] holds the next state in the "escape" case
//   (coal with parent of recomb branch)
void get_recomb_transition_switch(const ArgModel *model,
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    NodeStateLookup &state2_lookup,
    int next_states[2], int last_path, int minage)
{
    //const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;

    // make state lookup
    //NodeStateLookup state2_lookup(states2, nnodes);

    // SPR subtree moves out from underneath us
    // therefore therefore the new chromosome coalesces with
    // the branch above the subtree

    // search up for parent
    int parent = last_nodes[spr.recomb_node].parent;
    int time2 = last_nodes[parent].age;
    int node2;

    // find other child
    int other = last_tree->get_sibling(spr.recomb_node);

    // find new state in tree
    if (spr.coal_node == spr.recomb_node) {
        node2 = mapping[spr.coal_node];
        time2 = spr.coal_time;
    } else node2 = (other == spr.coal_node ?
                    nodes[mapping[other]].parent : mapping[other]);

    // stay case
    next_states[0] = state2_lookup.lookup(mapping[spr.recomb_node],
                                          spr.recomb_time, last_path);

    // escape case
    next_states[1] = state2_lookup.lookup(node2, time2,
                                          model->consistent_path(last_path,
                                              last_nodes[spr.recomb_node].pop_path,
                                              minage, spr.recomb_time, time2));
}


// note that recomb_parent_age and last_treelen are dependent on
// the coalescence point at last_tree
double calc_recomb(
    const LocalTree *last_tree, const ArgModel *model,
    const LineageCounts *lineages,
    const Spr &spr, const State &state1,
    double last_treelen,
    const bool internal=false)
{
    // get times
    int a = state1.time;
    const int k = spr.recomb_time;
    double last_treelen_b;
    int minage=0;
    int maintree_root=-1, subtree_root=-1;

    int root_age;
    if (internal) {
        subtree_root = last_tree->nodes[last_tree->root].child[0];
        maintree_root = last_tree->nodes[last_tree->root].child[1];

        // THIS IS A HACK.. need to revisit best way to deal with recombs on root node
        if (spr.recomb_node == spr.coal_node &&
            spr.recomb_node == maintree_root) {
            return 1.0;
        }
        minage = last_tree->nodes[subtree_root].age;
        root_age = last_tree->nodes[maintree_root].age;

        assert(spr.recomb_node != subtree_root);

        // detect sprs onto subtree root branch
        if (spr.coal_node == subtree_root) {
            if (a < spr.coal_time)
                return 0.0;
            if (spr.recomb_node == maintree_root) {
                if (state1.node != maintree_root) {
                    return 0.0;
                }
            }
        }

        // detect sprs from anywhere in maintree to anywhere in subtree
        // state1 must be topologically above recomb point
        if (state1.time == spr.recomb_time) {
            int ptr = spr.coal_node;
            int ptr2 = -1, ptr3 = -1;
            while (ptr != last_tree->root) {
                ptr2 = ptr;
                ptr = last_tree->nodes[ptr].parent;
            }
            // now ptr2 is either subtree_root or maintree_root,
            // whichever is above coal_node
            ptr = spr.recomb_node;
            while (ptr != last_tree->root) {
                ptr3 = ptr;
                ptr = last_tree->nodes[ptr].parent;
            }
            // now ptr3 is either subtree_root or maintree_root, whichever is
            // above recomb_node
            if (ptr2 == subtree_root && ptr3 == maintree_root) {
                // check is needed, ensure state1 is not topologically below recomb
                ptr = last_tree->nodes[state1.node].parent;
                while (ptr != last_tree->root) {
                    if (ptr == spr.recomb_node)
                        return 0.0;
                    ptr = last_tree->nodes[ptr].parent;
                }
            }
        }

        last_treelen += model->times[a] -
            model->times[last_tree->nodes[subtree_root].age];

        if (a > root_age) {
            // add wrapped branch
            last_treelen += model->times[a] - model->times[root_age];

            // add basal branch
            last_treelen_b = last_treelen + model->time_steps[a];
        } else {
            // add basal branch
            last_treelen_b = last_treelen + model->time_steps[root_age];
        }
    } else {  //external
        root_age = last_tree->nodes[last_tree->root].age;
        last_treelen = get_treelen_branch(
            last_tree, model->times, model->ntimes,
            state1.node, state1.time, last_treelen, false);
        last_treelen_b = last_treelen + get_basal_branch(
            last_tree, model->times, model->ntimes,
            state1.node, state1.time);
    }

    double p;
    if (model->smc_prime) {
        int pop = model->get_pop(spr.pop_path, k);
        int nbranches_above = lineages->nbranches_pop[pop][2*k];
        int nbranches_below = (k == 0 ? 0 : lineages->nbranches_pop[pop][2*k-1]);
        int pop2 = model->get_pop(state1.pop_path, k);
        int new_root_age = max(a, root_age);
        int nrecomb = lineages->ncoals_pop[pop][k];
        assert(k <= new_root_age);
        if (k == new_root_age) nbranches_above = 0;
        if (pop == pop2) {
            if (k > minage && k <= a)
                nbranches_below++;
            if (k >= minage && k < a)
                nbranches_above++;
            if (k == a) {
                nrecomb += 2;
            } else if (k >= minage && k < a) {
                nrecomb++;
            }
        }
        double blen_above = model->coal_time_steps[2*k]
            * (double)nbranches_above;
        double blen_below = (k == 0 ? 0.0 :
                             model->coal_time_steps[2*k-1]
                             * (double)nbranches_below);
        if (last_treelen == 0) return 0.0;
        p = (blen_above + blen_below)
            * (1.0 - exp(-max(model->rho * last_treelen, model->rho)))
            / (last_treelen * (double)nrecomb);
        assert(!isnan(p));
    } else {
        // probability of recombination rate and location
        int nbranches_k = lineages->nbranches[k] + int(k < a && k >= minage);
        int nrecombs_k = lineages->nrecombs[k] + int(k <= a && k >= minage) +
            int(k == a && k >= minage)  - int(k >= max(root_age, a));
        p = nbranches_k * model->time_steps[k] /
            (nrecombs_k * last_treelen_b) *
            (1.0 - exp(-max(model->rho * last_treelen, model->rho)));
        if (nrecombs_k <= 0 || nbranches_k <= 0) {
            printError("counts %d %d %e\n",
                       nrecombs_k, nbranches_k, p);
            assert(false);
        }
    }
    //    printf("calc_recomb %e\n", p);fflush(stdout);
    return p;
}


//calc_recoal_sums sets *sums to the probability of not coalescing between
// spr.recomb_time and spr.coal_time, but does not consider the extra
// branch that may exist depending on where the current branch is being
// threaded.
// sums2[i] represents the probability of the broken branch not coalescing
//   at (half-interval) time i within the spr path so that it can be added
//   to the log probability of not coalescing with the rest of the tree in
// cases where this extra branch exists and is in the same population
// both sums only pertain to the population path taken by the Spr event
void calc_recoal_sums(const ArgModel *model, const LineageCounts *lineages,
                      const Spr &spr, const int recomb_parent_age,
                      const int recomb_node_path,
                      double *sums, double **sums2,
                      int minage)
{
    // get info
    const int k = spr.recomb_time;
    const int j = spr.coal_time;
    //    printf("calc_recoal_sums minage=%i\n", minage);

    // probability of not coalescing before time j
    double sum = 0.0;
    for (int m=2*k; m<2*j-1; m++) {
        int spr_pop = model->get_pop(spr.pop_path, (m+1)/2);
        int recomb_node_pop = model->get_pop(recomb_node_path, (m+1)/2);
        int nbranches_m = lineages->nbranches_pop[spr_pop][m] -
            int((!model->smc_prime) &&
                m/2<recomb_parent_age &&
                spr_pop == recomb_node_pop);
        sum += (model->coal_time_steps[m] * nbranches_m
		/ (2.0 * model->popsizes[spr_pop][m]));
    }
    *sums = sum;

    for (int path=0; path < model->num_pop_paths(); path++) {
        sum = 0.0;
        sums2[path][2*k]=sum;
        for (int m=2*k; m<2*j-1; m++) {
            if (m >= 2*minage) {
                int spr_pop = model->get_pop(spr.pop_path, (m+1)/2);
                int this_pop = model->get_pop(path, (m+1)/2);
                if (spr_pop == this_pop)
                    sum += model->coal_time_steps[m] / (2.0 * model->popsizes[spr_pop][m]);
            }
            sums2[path][m+1] = sum;
        }
    }
}


double calc_recoal(
    const LocalTree *last_tree, const ArgModel *model,
    const LineageCounts *lineages,
    const Spr &spr, const State &state1,
    const int recomb_parent_age,
    const int recomb_parent_path,
    const bool internal=false)
{

    // get times
    const int a = state1.time;
    const int k = spr.recomb_time;
    const int j = spr.coal_time;
    const int a_pop = model->get_pop(state1.pop_path, j);
    const int coal_pop = model->get_pop(spr.pop_path, j);

    const int *nbranches = lineages->nbranches_pop[coal_pop];
    const int *ncoals = lineages->ncoals_pop[coal_pop];
    const int recomb_parent_pop = model->get_pop(recomb_parent_path, j);
    int subtree_root=-1;
    int minage=0;
    //    static int function_count=0;
    if (internal) {
        subtree_root = last_tree->nodes[last_tree->root].child[0];
        minage = last_tree->nodes[subtree_root].age;
    }

    // probability of coalescing on chosen branch
    int nbranches_j = nbranches[2*j];
    if (!model->smc_prime) {
        nbranches_j -= int(j < recomb_parent_age && recomb_parent_pop == coal_pop);
        nbranches_j += int(j < a && a_pop == coal_pop && j >= minage);
    }
    nbranches_j += int(j <= a && a_pop == coal_pop && j >= minage);
    int ncoals_j = ncoals[j]
        - int((!model->smc_prime) &&
              j <= recomb_parent_age && recomb_parent_pop == coal_pop)
        - int((!model->smc_prime) &&
              j == recomb_parent_age && recomb_parent_pop == coal_pop)
        + int(j <= a && a_pop == coal_pop && j >= minage)
        + int(j == a && a_pop == coal_pop && j >= minage);
    bool over = false;

    if (internal) {
        int maintree_root = last_tree->nodes[last_tree->root].child[1];
        if (spr.recomb_node == maintree_root) {
            // special cases for properly calculating ncoals_j and nbranches_j
            if (spr.coal_time >= last_tree->nodes[subtree_root].age) {
                over = true;
                nbranches_j = 1;
                ncoals_j++;
            }
        }
    }
    if (nbranches_j == 0) return 0.0;
    double p = 1.0 / ncoals_j;

    // probability of coalescing in time interval j
    int b1=0;
    if (j < model->ntimes - 2) {
        double Z = 0.0;
        if (j>k) {
            b1 = nbranches[2*j-1]
                - int((!model->smc_prime) &&
                      j-1 < recomb_parent_age && recomb_parent_pop == coal_pop)
                + int(j-1 < a && a_pop == coal_pop && j-1 >= minage);
            if (over)
                b1 = 1;
            Z = model->coal_time_steps[2*j-1] *  b1 /
                (2.0 * model->popsizes[coal_pop][2*j-1]);
        }

        p *= 1.0 - exp(- model->coal_time_steps[2*j] * nbranches_j /
                       (2.0 * model->popsizes[coal_pop][2*j])
                       - Z);
    }

    // asserts
    if (ncoals_j <= 0) {
        assert(false);
        return 0.0;
        printError("counts %d %d %e\n",
                   ncoals_j, nbranches_j, p);
    }
    assert(!isnan(p) && p>0);
    /*    printf("%i calc_recoal %e %i %i %i %i %i %i %i\n", function_count++, p, a, k, j,
           nbranches_j, ncoals_j, b1, over); fflush(stdout);
    printf("j=%i minage=%i rpa=%i nb=%i int=%i\n",
    j, minage, recomb_parent_age, nbranches[2*j-1], internal);*/
    return p;
}


double calc_recomb_recoal(
    const LocalTree *last_tree, const ArgModel *model,
    const LineageCounts *lineages,
    const Spr &spr, const State state1,
    const int recomb_parent_age, const int recomb_parent_path,
    double last_treelen,
    const bool internal=false)
{
    // get times
    int a = state1.time;
    const int k = spr.recomb_time;
    const int j = spr.coal_time;

    //    printf("calc_recomb_recoal internal=%i\n", internal);
    int minage=0;
    if (internal) {
        int subtree_root = last_tree->nodes[last_tree->root].child[0];
        minage = last_tree->nodes[subtree_root].age;
    }

    // calculate probability of recombination
    double p = calc_recomb(last_tree, model, lineages, spr, state1,
                           last_treelen, internal);

    // probability of not coalescing before time j
    double sum = 0.0;
    int a_path = state1.pop_path;
    int spr_path = spr.pop_path;
    int recomb_node_path = last_tree->nodes[spr.recomb_node].pop_path;
    int a_pop = model->get_pop(a_path, k);
    int spr_pop = model->get_pop(spr_path, k);
    int recomb_node_pop = model->get_pop(recomb_node_path, k);
    for (int m=2*k; m<2*j-1; m++) {
        if (m%2 == 1 && model->pop_tree != NULL) {
            int t = (m+1)/2;
            a_pop = model->get_pop(a_path, t);
            spr_pop = model->get_pop(spr_path, t);
            recomb_node_pop = model->get_pop(recomb_node_path, t);
        }
        int nbranches_m = lineages->nbranches_pop[spr_pop][m]
            - int((!model->smc_prime) && m/2<recomb_parent_age && spr_pop == recomb_node_pop)
            + int(m/2 < a && spr_pop == a_pop && m/2 >= minage);
        /*        printf("%i %i %i\n", m, lineages->nbranches_pop[spr_pop][m],
                  nbranches_m);*/
        sum += model->coal_time_steps[m] * nbranches_m /
            (2.0 * model->popsizes[spr_pop][m]);
    }
    p *= exp(-sum);

    p *= calc_recoal(last_tree, model, lineages, spr, state1,
                     recomb_parent_age, recomb_parent_path, internal);
    p *= model->path_prob(spr.pop_path, spr.recomb_time,
                          spr.coal_time);
    /*    printf("calc_recomb_recoal %f %f %f\n", exp(-sum),
           calc_recoal(last_tree, model, lineages, spr, state1,
                       recomb_parent_age, recomb_parent_path, internal),
           model->path_prob(spr.pop_path, spr.recomb_time,
                            spr.coal_time)); fflush(stdout);*/
    return p;
}


// Calculate the switch transition matrix.
void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages,
    TransMatrixSwitch *transmat_switch, bool internal)
{
    const int nstates1 = states1.size();
    const int nstates2 = states2.size();
    int recomb_parent_age;
    int recomb_parent_path;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int root = tree->root;
    const int last_root = last_tree->root;
    const int maintree_root = internal ? nodes[root].child[1] : -1;
    const int subtree_root = internal ? nodes[root].child[0] : -1;
    //    const int last_maintree_root = internal ? last_nodes[last_root].child[1] : -1;
    const int last_subtree_root = internal ? last_nodes[last_root].child[0] : -1;
    const int minage1 = internal ? last_nodes[last_subtree_root].age : 0;
    const int minage2 = internal ?  nodes[subtree_root].age : 0;
    static int count=0;
    count++;

    //    printf("calc_transition_probs_switch internal=%i\n", internal);
    for (int i=0; i < max(1,nstates1); i++)
        transmat_switch->recombsrc[i] = transmat_switch->recoalsrc[i] = -1;
    transmat_switch->init_probs();

    double last_treelen = (internal ?
        get_treelen_internal(last_tree, model->times, model->ntimes) :
        get_treelen(last_tree, model->times, model->ntimes, false));

    // handle the case of threading an internal branch
    if (internal) {
        // remove from the top case
        if (nstates1 == 0) {
            // switching between two completely specified blocks
            if (nstates2 == 0) {
                transmat_switch->determ[0] = 0;
                transmat_switch->determprob[0] = 1.0;
                return;
            }

            assert(spr.coal_node == last_root);
            int target_path=-1;
            if (subtree_root == mapping[spr.recomb_node]) {
                target_path = model->consistent_path(last_nodes[spr.recomb_node].pop_path,
                                                     spr.pop_path,
                                                     last_nodes[spr.recomb_node].age,
                                                     spr.recomb_time,
                                                     spr.coal_time);
            } else if (subtree_root == mapping[spr.coal_node]) {
                target_path = last_nodes[spr.coal_node].pop_path;
            } else {
                assert(maintree_root == mapping[spr.recomb_node]);
                assert(mapping[last_root] == -1);
                int inv_mapping=-1;
                for (int i=0; i < tree->nnodes; i++)
                    if (mapping[i] == subtree_root)
                        inv_mapping = i;
                assert(inv_mapping != -1);
                target_path = model->consistent_path(last_tree->nodes[inv_mapping].pop_path,
                                                     last_tree->nodes[last_root].pop_path,
                                                     last_tree->nodes[inv_mapping].age,
                                                     last_tree->nodes[last_root].age,
                                                     spr.coal_time);
            }

            for (int j=0; j<nstates2; j++) {
                if (states2[j].node == maintree_root &&
                    states2[j].time == spr.coal_time &&
                    model->paths_equal(states2[j].pop_path, target_path,
                                       nodes[subtree_root].age,
                                       spr.coal_time)) {
                    transmat_switch->determ[0] = j;
                    transmat_switch->determprob[0] = 1.0;
                    return;
                }
            }


            /*            if (spr.recomb_time >= tree->nodes[subtree_root].age)
                target_path =
                    model->consistent_path(tree->nodes[subtree_root].pop_path,
                                           spr.pop_path,
                                           tree->nodes[subtree_root].age,
                                           spr.recomb_time,
                                           spr.coal_time);*/

            for (int j=0; j<nstates2; j++) {
                if (states2[j].node == maintree_root &&
                    states2[j].time == spr.coal_time &&
                    model->paths_equal(states2[j].pop_path, target_path,
                                       nodes[subtree_root].age,
                                       spr.coal_time)) {
                    transmat_switch->determ[0] = j;
                    transmat_switch->determprob[0] = 1.0;
                    return;
                }
            }

            assert(false);
        }

        // fall off the top case
        if (nstates2 == 0) {
            fill(transmat_switch->determ, transmat_switch->determ+nstates1, 0);
            for (int i=0; i<nstates1; i++) {
                if (states1[i].node == spr.recomb_node &&
                    states1[i].time > spr.recomb_time) {
                    recomb_parent_age = states1[i].time;
                    recomb_parent_path = last_nodes[spr.recomb_node].pop_path;
                } else {
                    int recomb_parent = last_nodes[spr.recomb_node].parent;
                    recomb_parent_age = last_nodes[recomb_parent].age;
                    recomb_parent_path = last_nodes[spr.recomb_node].pop_path;
                }
                assert(mapping[last_subtree_root] != -1);
                if (model->paths_equal(states1[i].pop_path,
                                       nodes[mapping[last_subtree_root]].pop_path,
                                       last_nodes[last_subtree_root].age,
                                       nodes[mapping[last_subtree_root]].parent == -1 ?
                                       -1 :
                                       nodes[nodes[mapping[last_subtree_root]].parent].age))
                    transmat_switch->determprob[i] =
                        calc_recomb_recoal(last_tree, model, lineages, spr,
                                           states1[i], recomb_parent_age, recomb_parent_path,
                                           last_treelen, internal);
                else transmat_switch->determprob[i] = 0.0;
            }
            return;
        }
    }  // if internal


    // get deterministic transitions
    NodeStateLookup state2_lookup(states2, minage2, model->pop_tree);
    get_deterministic_transitions(last_tree, tree, spr, mapping,
                                  states1, states2, state2_lookup,
                                  model, transmat_switch->determ,
                                  internal);

    // calculate terms for deterministic probabilities
    recomb_parent_age = last_nodes[last_nodes[spr.recomb_node].parent].age;
    int numpath = (int)model->num_pop_paths();
    double sums;
    double **sums2 = new double* [numpath];
    for (int i=0; i < numpath; i++)
        sums2[i] = new double[2*model->ntimes];
    calc_recoal_sums(model, lineages, spr, recomb_parent_age,
                     last_nodes[spr.recomb_node].pop_path,
                     &sums, sums2, minage1);
    int npop = model->num_pops();
    int recoals_size = model->ntimes * npop;
    double recoals[recoals_size];
    for (int a=0; a < recoals_size; a++)
        recoals[a] = -1.0;

    for (int i=0; i<nstates1; i++) {
        int j = transmat_switch->determ[i];
        if (j >= 0) {

            // deterministic
            if (states1[i].node == spr.recomb_node &&
                states1[i].time > spr.recomb_time) {
                // prev coal on recomb branch but at greater time- deterministic
                // coal at parent
                recomb_parent_age = states1[i].time;
                recomb_parent_path = last_nodes[spr.recomb_node].pop_path;
                transmat_switch->determprob[i] =
                calc_recomb_recoal(last_tree, model, lineages, spr,
                                   states1[i], recomb_parent_age,
                                   recomb_parent_path,
                                   last_treelen, internal);
            } else {
                int pop = model->get_pop(states1[i].pop_path, states1[i].time);
                int recoal_pos = pop * model->ntimes + states1[i].time;
                int recomb_parent = last_nodes[spr.recomb_node].parent;
                recomb_parent_age = last_nodes[recomb_parent].age;
                recomb_parent_path = last_nodes[spr.recomb_node].pop_path;
                if (recoals[recoal_pos] < 0)
                    recoals[recoal_pos] = calc_recoal(last_tree, model, lineages,
                                                      spr, states1[i],
                                                      recomb_parent_age,
                                                      recomb_parent_path, internal);
                transmat_switch->determprob[i] =
                    calc_recomb(last_tree, model, lineages, spr, states1[i],
                                last_treelen, internal) *
                    exp(-sums
                        -sums2[states1[i].pop_path][
                           max(min(2*spr.coal_time-1, 2*states1[i].time),
                                   2*spr.recomb_time)]) *
                    recoals[recoal_pos] *
                    model->path_prob(spr.pop_path, spr.recomb_time,
                                     spr.coal_time);
                assert(!isnan(transmat_switch->determprob[i]) &&
                       transmat_switch->determprob[i] > 0);
                /*                printf("d2: %i %i %e %e %i %i %i\n",
                       i, j, -sums, -sums2[max(min(2*spr.coal_time-1, 2*states1[i].time),
                                               2*spr.recomb_time)],
                       2*spr.coal_time-1, 2*states1[i].time,
                       2*spr.recomb_time);*/
            }
        }
    }
    for (int i=0; i < numpath; i++)
        delete [] sums2[i];
    delete [] sums2;


    // find probabilitistic transition source states
    for (int i=0; i < nstates1; i++) {
        if (states1[i].node == spr.recomb_node &&
            states1[i].time == spr.recomb_time) {
            transmat_switch->recombsrc[i] = states1[i].pop_path;
            int recomb_next_states[2];
            get_recomb_transition_switch(model, tree, last_tree, spr, mapping,
                                         state2_lookup, recomb_next_states,
                                         states1[i].pop_path, minage1);
            int need_path = -1;
            if (internal) {
                int last_subroot = last_tree->nodes[last_tree->root].child[0];
                int subroot = tree->nodes[tree->root].child[0];
                int inv_mapping = mapping[last_subroot];
                if (inv_mapping != subroot) {
                    need_path = tree->nodes[inv_mapping].pop_path;
                    if (!model->paths_equal(need_path, states1[i].pop_path,
                                            last_tree->nodes[last_subroot].age,
                                            min(states1[i].time,
                                                nodes[inv_mapping].parent == -1 ? model->ntimes :
                                                nodes[nodes[inv_mapping].parent].age)))
                        continue;
                }
            }

            // stay case (recomb above)
            int j = recomb_next_states[0];
            if (j != -1) {
                recomb_parent_age = last_nodes[
                    last_nodes[spr.recomb_node].parent].age;
                recomb_parent_path = last_nodes[spr.recomb_node].pop_path;
                transmat_switch->set(i, j,
                    calc_recomb_recoal(last_tree, model, lineages, spr,
                                       states1[i], recomb_parent_age,
                                       recomb_parent_path,
                                       last_treelen, internal));
                assert(!isnan(transmat_switch->get(i, j)));
            }

            // escape case (recomb below)
            j = recomb_next_states[1];
            if (j != -1) {
                recomb_parent_age = states1[i].time;
                recomb_parent_path = last_nodes[states1[i].node].pop_path;
                transmat_switch->set(i, j,
                    calc_recomb_recoal(last_tree, model, lineages, spr, states1[i],
                                       recomb_parent_age, recomb_parent_path,
                                       last_treelen, internal));
                assert(!isnan(transmat_switch->get(i, j)));
            }
        } else if (states1[i].node == spr.coal_node &&
                   states1[i].time == spr.coal_time) {
            transmat_switch->recoalsrc[i] = states1[i].pop_path;
            int node1 = states1[i].node;
            int time1 = states1[i].time;

            // determine if node1 is still here or not
            int node3;
            int last_parent = last_nodes[spr.recomb_node].parent;
            if (last_parent == node1) {
                // recomb breaks node1 branch, we need to use the other child
                const int *c = last_nodes[last_parent].child;
                node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
            } else{
                node3 = mapping[node1];
            }

            int parent = nodes[mapping[spr.recomb_node]].parent;
            assert(parent == nodes[node3].parent);

            for (int j=0; j<nstates2; j++) {
                const int node2 = states2[j].node;
                const int time2 = states2[j].time;
                //                const int path2 = states2[j].pop_path;

                bool possible=false;
                if (time2 == time1 &&
                    (node2 == node3 ||
                     (spr.recomb_node != spr.coal_node && node2 == parent)) &&
                    model->paths_equal(states2[j].pop_path,
                                       states1[i].pop_path,
                                       max(minage1, minage2), time1)) {
                    possible=true;
                } else if (node2 == mapping[spr.recomb_node] &&
                           time2 >= spr.recomb_time &&
                           time2 <= time1) {
                    if ((!internal) ||
                        mapping[last_subtree_root] == subtree_root) {
                        if (model->paths_equal(states2[j].pop_path,
                                               spr.pop_path,
                                               time2, time1) &&
                            model->paths_equal(states1[i].pop_path,
                                               spr.pop_path, time2, time1) &&
                            model->paths_equal(states2[j].pop_path,
                                               states1[i].pop_path,
                                               max(minage1, minage2), time2)) {
                            possible=true;
                        }
                    } else if (internal) {
                        int inv_mapping;
                        int c[2];
                        c[0] = last_nodes[last_subtree_root].child[0];
                        c[1] = last_nodes[last_subtree_root].child[1];
                        assert(c[0] == spr.recomb_node ||
                               c[1] == spr.recomb_node);
                        //find the node that corresponds to the subtree root in new tree
                        if (c[0] == spr.recomb_node) {
                            assert(mapping[c[1]] == subtree_root);
                            inv_mapping = c[1];
                        } else {
                            assert(mapping[c[0]] == subtree_root);
                            inv_mapping = c[0];
                        }
                        assert(nodes[subtree_root].age ==
                               last_nodes[inv_mapping].age);
                        if (states2[j].time > last_nodes[last_subtree_root].age) {
                            if (model->paths_equal(last_nodes[inv_mapping].pop_path,
                                                   states2[j].pop_path,
                                                   nodes[subtree_root].age,
                                                   last_nodes[last_subtree_root].age) &&
                                model->paths_equal(states1[i].pop_path,
                                                   states2[j].pop_path,
                                                   last_nodes[last_subtree_root].age,
                                                   states2[j].time))
                                possible=true;
                        } else {
                            if (model->paths_equal(last_nodes[inv_mapping].pop_path,
                                                   states2[j].pop_path,
                                                   nodes[subtree_root].age,
                                                   states2[j].time) &&
                                model->paths_equal(last_nodes[inv_mapping].pop_path,
                                                   nodes[states2[j].node].pop_path,
                                                   states2[j].time,
                                                   states1[i].time))
                                possible=true;
                        }
                    }
                }
                if (!possible) continue;

                // if internal, check if minage changed, implying that
                // broken branch is below
                // path must match along broken b
                if (model->pop_tree != NULL) {
                    assert(time2 <= time1);
                    assert(time2 >= minage2);
                    if ((!internal) || mapping[last_subtree_root] == subtree_root) {
                        assert(minage1 == minage2);
                        if (!model->paths_equal(states1[i].pop_path,
                                                states2[j].pop_path,
                                                minage1, time2)) {
                            continue;
                        }
                    } else {   // subtree_root has changed
                        assert(minage2 <= minage1);
                        assert(last_subtree_root == last_parent);
                        int other;
                        if (last_nodes[last_subtree_root].child[0] ==
                            spr.recomb_node) {
                            other = last_nodes[last_subtree_root].child[1];
                        } else {
                            assert(last_nodes[last_subtree_root].child[1] ==
                                   spr.recomb_node);
                            other = last_nodes[last_subtree_root].child[0];
                        }
                        assert(last_nodes[other].age == minage2);
                        int path1 =
                            model->consistent_path(last_nodes[other].pop_path,
                                                   states1[i].pop_path,
                                                   minage2, minage1,
                                                   time1, false);
                        if (path1 == -1) continue;

                        int path2 =
                            model->consistent_path(states2[j].pop_path,
                                                   spr.pop_path,
                                                   minage2, time2, time1, false);
                        if (path2 == -1) continue;
                        if (! model->paths_equal(path1, path2, minage2, time1))
                            continue;
                    }
                }

                // in all cases check that path is same between two states
                // for times where they overlap
                if (max(minage1, minage2) <= min(time1, time2)) {
                    assert(model->paths_equal(states1[i].pop_path,
                                              states2[j].pop_path,
                                              max(minage1, minage2),
                                              min(time1, time2)));
                }

                recomb_parent_age = last_nodes[last_nodes[spr.recomb_node].parent].age;
                recomb_parent_path = last_nodes[spr.recomb_node].pop_path;
                Spr spr2 = spr;
                spr2.coal_time = time2;
                transmat_switch->set(i, j,
                    calc_recomb_recoal(last_tree, model, lineages, spr2, states1[i],
                                       recomb_parent_age, recomb_parent_path,
                                       last_treelen, internal));
                assert(!isnan(transmat_switch->get(i, j)));
            }
        }
    }

    bool found=false;
    for (int i=0; i < nstates1; i++) {
        if (transmat_switch->determ[i] >= 0) found=true;
        if (transmat_switch->recombsrc[i] >= 0) found=true;
        if (transmat_switch->recoalsrc[i] >= 0) found=true;
        if (found) break;
    }
    assert(found);

}


void calc_transition_probs_switch_internal(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages,
    TransMatrixSwitch *transmat_switch)
{
    calc_transition_probs_switch(tree, last_tree, spr, mapping,
                                 states1, states2,
                                 model, lineages,
                                 transmat_switch, true);
}


void get_transition_probs_switch(const TransMatrixSwitch *matrix,
                                 double **transprob)
{
    for (int i=0; i<matrix->nstates1; i++) {
        for (int j=0; j<matrix->nstates2; j++) {
            transprob[i][j] = matrix->get_log(i, j);
        }
    }
}


void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, double **transprob)
{
    TransMatrixSwitch transmat_switch(states1.size(), states2.size(),
                                      model->num_pop_paths());
    calc_transition_probs_switch(tree, last_tree, spr, mapping,
                                 states1, states2, model, lineages,
                                 &transmat_switch);
    get_transition_probs_switch(&transmat_switch, transprob);
}


//=============================================================================
// prior for state space


double calc_state_priors(const ArgModel *model,
                         int time, int pop_path,
                         int **nbranches_pop,
                         int **ncoals_pop, int minage)
{
    const int b = time;

    if (b < minage)
        return 0.0;

    // probability of not coalescing before time b
    double sum = 0.0;
    int pop=model->get_pop(pop_path, minage);
    for (int m=2*minage; m<2*b-1; m++) {
        if (m%2==1) pop = model->get_pop(pop_path, (m+1)/2);
        sum += model->coal_time_steps[m] * nbranches_pop[pop][m]
            / (2.0 * model->popsizes[pop][m]);
    }
    pop = model->get_pop(pop_path, b);
    double p = exp(-sum) / ncoals_pop[pop][b];

    // probability of coalescing in time interval b
    if (b < model->ntimes - 2) {
        double Z = 0.0;
        if (b>minage)
            Z = model->coal_time_steps[2*b-1] *
                nbranches_pop[pop][2*b-1]/
                ( 2.0 * model->popsizes[pop][2*b-1] );
        p *= 1.0 - exp(- (model->coal_time_steps[2*b] *
                          nbranches_pop[pop][2*b])/
                       (2.0 * model->popsizes[pop][2*b]) - Z);
    } else {
        // b = ntimes -1, guaranteed coalescence
    }
    if (model->pop_tree != NULL)
        p *= model->path_prob(pop_path, minage, b);
    return p;
}


void calc_state_priors(const States &states, const LineageCounts *lineages,
                       const ArgModel *model, double *priors,
                       const int minage)
{
    const int nstates = states.size();

    // special case
    if (nstates == 0) {
        priors[0] = 1.0;
        return;
    }

    for (int i=0; i<nstates; i++)
        priors[i] = calc_state_priors(model,
            states[i].time, states[i].pop_path,
            lineages->nbranches_pop, lineages->ncoals_pop, minage);
}


void calc_state_priors_log(const States &states, const LineageCounts *lineages,
                       const ArgModel *model, double *priors,
                       const int minage)
{
    const int nstates = max((int) states.size(), 1);
    calc_state_priors(states, lineages, model, priors, minage);
    for (int i=0; i<nstates; i++)
        priors[i] = log(priors[i]);
}



//=============================================================================
// assert transition matrices


/*
double calc_recomb_prob(const LocalTree *tree, const ArgModel *model)
{
    double treelen = get_treelen(tree, model->times, model->ntimes, false);
    return 1.0 - exp(-max(model->rho * treelen, model->rho));
}


bool assert_transmat(const LocalTree *tree, const ArgModel *model,
                     const TransMatrix *matrix)
{
    const int nleaves = tree->get_num_leaves();
    const int nnodes = tree->nnodes;
    const int nnodes2 = nnodes + 2;


    LineageCounts lineages(model->ntimes, model->num_pops());
    States states;
    get_coal_states(tree, model->ntimes, states);
    int nstates = states.size();
    int displace[nnodes2];
    int newleaf = nleaves;

    // get local tree we can modify
    LocalTree tree2(tree->nnodes, tree->capacity + 2);
    tree2.copy(*tree);


    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            State state1 = states[i];
            State state2 = states[j];
            Spr spr;
            int sister_age = tree2.nodes[state1.node].age;

            add_tree_branch(&tree2, state1.node, state1.time, state1.pop_path);

            double p = -INFINITY;

            // recomb could be on new branch
            // recoal is state2
            spr.recomb_node = newleaf;
            spr.coal_node = state2.node;
            spr.coal_time = state2.time;

            // sum over possible recombination times
            for (int rtime=0; rtime<=min(state1.time, state2.time); rtime++) {
                spr.recomb_time = rtime;
                p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
            }


            if (state1.node == state2.node) {
                // recomb could be state1.node
                // recoal is on new branch or parent of state1.node
                spr.recomb_node = state1.node;

                if (state2.time < state1.time) {
                    spr.coal_node = newleaf;
                } else if (state2.time >= state1.time) {
                    spr.coal_node = tree2.nodes[newleaf].parent;
                }
                spr.coal_time = state2.time;

                // sum over possible recombination times
                for (int rtime=sister_age;
                     rtime<=min(state1.time, state2.time); rtime++) {
                    spr.recomb_time = rtime;
                    p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
                }
            }

            // calculate probability of recombination
            double recomb_prob = calc_recomb_prob(&tree2, model);
            p += log(recomb_prob);

            if (i == j)
                p = logadd(p, log(1.0 - recomb_prob));


            double p2 = matrix->get_log(tree, states, i, j);

            // compare probabilities
            printLog(LOG_MEDIUM, "> (%d,%d)->(%d,%d): %e = %e; %e\n",
                     state1.node, state1.time,
                     state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9)) {
                return false;
            }

            remove_tree_branch(&tree2, newleaf, model, displace);
        }
    }

    return true;
}



bool assert_transmat_switch(const LocalTree *tree, const Spr &_spr,
                            const ArgModel *model,
                            const TransMatrixSwitch *matrix)
{
    const int nleaves = tree->get_num_leaves();
    const int nnodes = tree->nnodes;
    const int ntimes = model->ntimes;

    // get local tree we can modify
    LocalTree tree1(tree->nnodes, tree->capacity + 2);
    tree1.copy(*tree);

    // build next tree from last tree
    LocalTree tree2(nnodes);
    tree2.copy(tree1);
    apply_spr(&tree2, _spr, model->pop_tree);

    // define mapping
    int mapping[nnodes];
    for (int i=0; i<nnodes; i++)
        mapping[i] = i;
    mapping[tree1.nodes[_spr.recomb_node].parent] = -1;

    // get states
    States states1, states2;
    get_coal_states(&tree1, ntimes, states1);
    get_coal_states(&tree2, ntimes, states2);
    int nstates1 = states1.size();
    int nstates2 = states2.size();

    // get lineages
    LineageCounts lineages(ntimes, model->num_pops());

    // add/remove branch data
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;

    printLog(LOG_MEDIUM, "> matrix recombsrc=%d(%d,%d), recoalsrc=%d(%d,%d)\n",
             matrix->recombsrc,
             states1[matrix->recombsrc].node,
             states1[matrix->recombsrc].time,
             matrix->recoalsrc,
             states1[matrix->recoalsrc].node,
             states1[matrix->recoalsrc].time);

    for (int i=0; i<nstates1; i++) {
        for (int j=0; j<nstates2; j++) {
            State state1 = states1[i];
            State state2 = states2[j];
            double p = -INFINITY;
            Spr spr = _spr;

            double p2 = matrix->get_log(i, j);
            if (p2 == -INFINITY)
                continue;

            // add new branch and modify spr
            add_tree_branch(&tree1, state1.node, state1.time, state1.pop_path);
            add_tree_branch(&tree2, state2.node, state2.time, state1.pop_path);
            add_spr_branch(&tree2, &tree1, state2, state1,
                           &spr, mapping, newleaf, displaced, newcoal);

            // calculate probability of recombination
            double recomb_prob = calc_recomb_prob(&tree1, model);
            p = log(recomb_prob);
            p += calc_spr_prob(model, &tree1, spr, lineages);

            remove_tree_branch(&tree1, newleaf, model);
            remove_tree_branch(&tree2, newleaf, model);


            // compare probabilities
            printLog(LOG_MEDIUM, "> %d,%d (%d,%d)->(%d,%d): %e = %e; %e\n",
                     i, j,
                     state1.node, state1.time,
                     state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9)) {
                draw_local_tree(&tree1, stdout);
                return false;
            }
        }
    }

    return true;
}


bool assert_transmat_internal(const LocalTree *tree, const ArgModel *model,
                              const TransMatrix *matrix)
{
    LineageCounts lineages(model->ntimes, model->num_pops());
    States states;
    get_coal_states_internal(tree, model->ntimes, states);
    int nstates = states.size();
    int subtree_root = tree->nodes[tree->root].child[0];

    // get local tree we can modify
    LocalTree tree2(tree->nnodes);
    tree2.copy(*tree);


    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            State state1 = states[i];
            State state2 = states[j];
            Spr spr;
            int sister_age = tree2.nodes[state1.node].age;

            Spr add_spr(subtree_root, tree->nodes[subtree_root].age,
                        state1.node, state1.time);
            apply_spr(&tree2, add_spr, model->pop_tree);

            // calculate probability of recombination
            double recomb_prob = calc_recomb_prob(&tree2, model);

            double p = -INFINITY;

            // recomb could be on new branch
            // recoal is state2
            spr.recomb_node = subtree_root;
            spr.coal_node = state2.node;
            spr.coal_time = state2.time;
            assert(tree2[subtree_root].age <= spr.coal_time);

            // sum over possible recombination times
            for (int rtime=tree2[subtree_root].age;
                 rtime<=min(state1.time, state2.time); rtime++) {
                spr.recomb_time = rtime;
                p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
            }


            if (state1.node == state2.node) {
                // recomb could be state1.node
                // recoal is on new branch or parent of state1.node
                spr.recomb_node = state1.node;

                if (state2.time < state1.time) {
                    spr.coal_node = subtree_root;
                } else if (state2.time >= state1.time) {
                    spr.coal_node = tree2.nodes[subtree_root].parent;
                }
                spr.coal_time = state2.time;

                // sum over possible recombination times
                for (int rtime=sister_age;
                     rtime<=min(state1.time, state2.time); rtime++) {
                    spr.recomb_time = rtime;
                    p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
                }
            }

            p += log(recomb_prob);

            if (i == j)
                p = logadd(p, log(1.0 - recomb_prob));


            double p2 = matrix->get_log(tree, states, i, j);

            printLog(LOG_MEDIUM, "> (%d,%d)->(%d,%d): %e = %e; %e\n",
                     state1.node, state1.time,
                     state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9))
                return false;

            Spr remove_spr(subtree_root, tree2[subtree_root].age,
                           tree2.root, model->ntimes+1);
            apply_spr(&tree2, remove_spr, model->pop_tree);
        }
    }

    return true;
}


bool assert_transmat_switch_internal(const LocalTree *last_tree,
                                     const LocalTree *tree,
                                     const Spr &_spr, const int *_mapping,
                                     ArgModel *model,
                                     TransMatrixSwitch *transmat_switch)
{
    LineageCounts lineages(model->ntimes, model->num_pops());
    States states, last_states;
    get_coal_states_internal(last_tree, model->ntimes, last_states);
    get_coal_states_internal(tree, model->ntimes, states);
    int nstates1 = last_states.size();
    int nstates2 = states.size();
    int last_subtree_root = last_tree->nodes[last_tree->root].child[0];
    int subtree_root = tree->nodes[tree->root].child[0];

    for (int i=0; i<nstates1; i++) {
        for (int j=0; j<nstates2; j++) {
            State state1 = last_states[i];
            State state2 = states[j];

            double p2 = transmat_switch->get_log(i, j);
            if (p2 == -INFINITY)
                continue;

            // get local tree we can modify
            LocalTree last_tree2(last_tree->nnodes);
            last_tree2.copy(*last_tree);
            LocalTree tree2(tree->nnodes);
            tree2.copy(*tree);
            Spr spr = _spr;
            int mapping[tree->nnodes];
            for (int k=0; k<tree->nnodes; k++)
                mapping[k] = _mapping[k];


            // add branches
            Spr add_spr(last_subtree_root, last_tree2[last_subtree_root].age,
                        state1.node, state1.time);
            apply_spr(&last_tree2, add_spr, model->pop_tree);
            Spr add_spr2(subtree_root, tree2[subtree_root].age,
                         state2.node, state2.time);
            apply_spr(&tree2, add_spr2, model->pop_tree);
            add_spr_branch(&tree2, &last_tree2, state2, state1,
                           &spr, mapping, subtree_root, last_subtree_root);

            double recomb_prob = calc_recomb_prob(&last_tree2, model);
            double p = log(recomb_prob);
            p += calc_spr_prob(model, &last_tree2, spr, lineages);

            printLog(LOG_MEDIUM, "> (%d,%d)->(%d,%d): %e = %e; %e\n",
                     state1.node, state1.time,
                     state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9))
                return false;
        }
    }

    return true;
}
*/


//=============================================================================
// C interface
extern "C" {
    /*
double **new_transition_probs(int nnodes, int *ptree,
                              int *ages, double treelen,
                              intstate *istates, int nstates,
                              int ntimes, double *times, double *time_steps,
                              int *nbranches, int *nrecombs, int *ncoals,
                              double *popsizes, double rho)
{

    // setup model, local tree, states
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);
    LineageCounts lineages(ntimes, model.num_pops());
    lineages.count(&tree, model.pop_tree);
    States states;
    make_states(istates, nstates, states);

    double **transprob = new_matrix<double>(nstates, nstates);
    calc_transition_probs(&tree, &model, states, &lineages, transprob);
    return transprob;
}


double **new_transition_probs_switch(
    int *ptree, int *last_ptree, int nnodes,
    int recomb_node, int recomb_time, int coal_node, int coal_time,
    int *ages_index, int *last_ages_index,
    double treelen, double last_treelen,
    intstate *istates1, int nstates1,
    intstate *istates2, int nstates2,

    int ntimes, double *times, double *time_steps,
    int *nbranches, int *nrecombs, int *ncoals,
    double *popsizes, double rho)
{
    // setup model
    ArgModel model(ntimes, times, popsizes, rho, 0.0);

    // setup local trees
    LocalTree tree(ptree, nnodes, ages_index);
    LocalTree last_tree(last_ptree, nnodes, last_ages_index);
    Spr spr(recomb_node, recomb_time, coal_node, coal_time);
    int mapping[nnodes];
    make_node_mapping(last_ptree, nnodes, recomb_node, mapping);
    LineageCounts lineages(ntimes, model.num_pops());
    lineages.count(&last_tree, model.pop_tree);

    // setup states
    States states1, states2;
    make_states(istates1, nstates1, states1);
    make_states(istates2, nstates2, states2);

    double **transprob = new_matrix<double>(nstates1, nstates2);
    calc_transition_probs_switch(&tree, &last_tree, spr, mapping,
        states1, states2, &model, &lineages, transprob);
    return transprob;
}


void delete_transition_probs(double **transmat, int nstates)
{
    delete_matrix<double>(transmat, nstates);
}


bool arghmm_assert_transmat(int nnodes, int *ptree, int *ages,
                            int ntimes, double *times,
                            double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);

    States states;
    get_coal_states(&tree, ntimes, states);

    LineageCounts lineages(ntimes);
    lineages.count(&tree, model.pop_tree);

    TransMatrix transmat(&model, states.size());
    calc_transition_probs(&tree, &model, states, &lineages, &transmat);

    return assert_transmat(&tree, &model, &transmat);
}



bool arghmm_assert_transmat_switch(
    int nnodes, int *ptree, int *ages,
    int recomb_node, int recomb_time, int coal_node, int coal_time,
    int ntimes, double *times,
    double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree last_tree(ptree, nnodes, ages);

    // build next tree from last tree
    LocalTree tree(nnodes);
    tree.copy(last_tree);
    Spr spr(recomb_node, recomb_time, coal_node, coal_time);
    apply_spr(&tree, spr);

    // define mapping
    int mapping[nnodes];
    for (int i=0; i<nnodes; i++)
        mapping[i] = i;
    mapping[last_tree.nodes[spr.recomb_node].parent] = -1;

    // get states
    States states, last_states;
    get_coal_states(&last_tree, ntimes, last_states);
    get_coal_states(&tree, ntimes, states);

    // get lineages
    LineageCounts lineages(ntimes);
    lineages.count(&last_tree, model.pop_tree);

    TransMatrixSwitch transmat_switch(last_states.size(), states.size(),
                                      model.num_pop_paths());
    calc_transition_probs_switch(&tree, &last_tree,
        spr, mapping, last_states, states, &model, &lineages,
        &transmat_switch);

    return assert_transmat_switch(&last_tree, spr, &model, &transmat_switch);
}


bool arghmm_assert_transmat_internal(int nnodes, int *ptree, int *ages,
                                     int ntimes, double *times,
                                     double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);

    // randomly choose branch to remove
    int node = irand(nnodes);
    Spr remove_spr(node, tree[node].age, tree.root, ntimes+1);
    apply_spr(&tree, remove_spr);
    int *c = tree[tree.root].child;
    if (c[0] != node) {
        int tmp = c[0];
        c[0] = c[1];
        c[1] = tmp;
    }


    States states;
    get_coal_states_internal(&tree, ntimes, states);

    LineageCounts lineages(ntimes);
    lineages.count(&tree, model.pop_tree, true);

    TransMatrix transmat(&model, states.size());
    transmat.calc_transition_probs(&tree, &model, states, &lineages, true);

    return assert_transmat_internal(&tree, &model, &transmat);
}


bool arghmm_assert_transmat_switch_internal(
    LocalTrees *trees, int ntimes, double *times, double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);

    const int maxtime = model.ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    int node = irand(trees->nnodes);
    sample_arg_removal_path(trees, node, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime);

    // build matrices
    ArgHmmMatrixIter matrix_iter(&model, NULL, trees);
    matrix_iter.set_internal(true);

    LocalTree const *last_tree = NULL;
    for (matrix_iter.begin(); matrix_iter.more(); matrix_iter.next()) {
        ArgHmmMatrices &matrices = matrix_iter.ref_matrices();
        const LocalTreeSpr *tree_spr = matrix_iter.get_tree_spr();
        const LocalTree *tree = tree_spr->tree;

        if (matrices.transmat_switch) {
            if (!assert_transmat_switch_internal(last_tree, tree, tree_spr->spr,
                                                 tree_spr->mapping, &model,
                                                 matrices.transmat_switch))
                return false;
        }

        last_tree = tree;
    }

    return true;
}
*/

} // extern C

} // namespace argweaver


