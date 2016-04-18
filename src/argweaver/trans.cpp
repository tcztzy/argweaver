//=============================================================================
// transitions

#include "matrices.h"
#include "total_prob.h"
#include "thread.h"
#include "trans.h"


namespace argweaver {


void calc_coal_rates_partial_tree(const ArgModel *model, const LocalTree *tree,
                                  const LineageCounts *lineages,
                                  double ***coal_rates, bool *do_path)
{
    for (int path=0; path < model->num_pop_paths(); path++) {
        if (do_path != NULL && do_path[path] == false) continue;
        for (int path2=0; path2 < model->num_pop_paths(); path2++) {
            if (do_path != NULL && do_path[path2] == false) continue;
            coal_rates[path][path2][-1] = 0.0;
            for (int i=0; i<2*model->ntimes-1; i++) {
                int pop = model->get_pop(path, (i+1)/2);
                int pop2 = model->get_pop(path2, (i+1)/2);
                int nbranch = lineages->nbranches_pop[pop][i/2];
                if (pop != pop2) nbranch--;
                coal_rates[path][path2][i] = model->coal_time_steps[i] *
                    nbranch / (2.0 * model->popsizes[pop][i]);
            }
            coal_rates[path][path2][2 * model->ntimes - 1] = 1.0;
        }
    }
}


// Calculate transition probability within a local block
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix,
    bool internal, int minage)
{
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
    if (num_paths > 1) {
        for (int i=0; i < num_paths; i++)
            have_pop_path[i]=false;
        for (unsigned int i=0; i < states.size(); i++)
            have_pop_path[states[i].pop_path] = true;
    } else have_pop_path[0] = true;

    // set internal branch resampling flag
    matrix->internal = internal;
    matrix->minage = minage;

    // get coalescent rates for each time sub-interval
    double coal_rates_alloc[num_paths * num_paths * (2*ntimes+1)];
    double **coal_rates[num_paths];
    int idx=1;
    for (int i=0; i < num_paths; i++) {
        coal_rates[i] = new double* [num_paths];
        for (int j=0; j < num_paths; j++) {
            coal_rates[i][j] = &coal_rates_alloc[idx];
            idx += 2 * ntimes + 1;
        }
    }
    calc_coal_rates_partial_tree(model, tree, lineages, coal_rates,
                                 have_pop_path);

    // compute cumulative coalescent rates
    double C_alloc[num_paths * num_paths * (2*ntimes+2)];
    double **C[num_paths];
    idx=2;
    for (int i=0; i < num_paths; i++) {
        C[i] = new double* [num_paths];
        for (int j=0; j < num_paths; j++) {
            if (have_pop_path[i] && have_pop_path[j]) {
                C[i][j] = &C_alloc[idx];
                C[i][j][-2] = 0.0;
                C[i][j][-1] = 0.0;
                for (int b=0; b < 2*ntimes - 1; b++)
                    C[i][j][b] = C[i][j][b-1] + coal_rates[i][j][b];
            }
            idx += 2*ntimes + 2;
        }
    }

    // determine tree information: root, root age, tree length
    int root_age_index;
    double root_age;
    double treelen;
    if (internal) {
        const int *c = tree->nodes[tree->root].child;
        const int subtree_root = c[0];
        const int maintree_root = c[1];
        const double subtree_age = times[tree->nodes[subtree_root].age];
        root_age_index = tree->nodes[maintree_root].age;
        root_age = times[root_age_index];
        // NOTE: subtree_age is discounted in advance to offset +times[b]
        treelen = get_treelen_internal(tree, times, ntimes) - subtree_age;
        matrix->minage = max(matrix->minage, tree->nodes[subtree_root].age);
    } else {
        root_age_index = tree->nodes[tree->root].age;
        root_age = times[root_age_index];
        treelen = get_treelen(tree, times, ntimes, false);
    }

    // calculate transition matrix terms
    for (int b=0; b<ntimes-1; b++) {
        // get tree length
        double treelen2 = treelen + times[b];
        double treelen2_b;
        if (b > root_age_index) {
            // add wrapped branch
            treelen2 += times[b] - root_age;

            // add basal branch
            treelen2_b = treelen2 + time_steps[b];
        } else {
            // add basal branch
            treelen2_b = treelen2 + time_steps[root_age_index];
        }
        matrix->norecombs[b] = exp(-max(rho * treelen2, rho));
        matrix->D[b] = (1.0 - exp( -rho * treelen2)) / treelen2_b;
    }

    for (int path=0; path < num_paths; path++) {
        if (! have_pop_path[path]) continue;
        for (int b=0; b < ntimes-1; b++) {
            double path_prob = model->path_prob(path, 0, b);
            matrix->path_prob[path][b] = path_prob;
            int pop = model->get_pop(path, b);
            for (int path2=0; path2 < num_paths; path2++) {
                if (! have_pop_path[path2]) continue;
                //                int pop2 = model->get_pop(path2, b);
                double term = C[path][path2][2*b-1] + log(
                   time_steps[b] * (nbranches[b] + 1.0) /
                   (nrecombs[b] + 1.0) / path_prob);
                if (b == 0)
                    matrix->lnB[path][path2][b] = term;
                else
                    matrix->lnB[path][path2][b] = logadd(matrix->lnB[path][path2][b-1], term);
                matrix->lnE2[path][path2][b] = -C[path][path2][2*b-2] +
                   (b < ntimes - 2 ?
                    log(1 - exp(-coal_rates[path][path2][2*b]-coal_rates[path][path2][2*b-1])) : 0.0);
                matrix->lnNegG1[path][path2][b] =
                    C[path][path2][2*b-1] + log( - time_steps[b] / path_prob * (
                     (nbranches[b] / (nrecombs[b] + 1.0 + int(b < root_age_index)))
                      - (nbranches[b] + 1.0) / (nrecombs[b] + 1.0)));
            }
            matrix->G2[path][b] = (b<ntimes-2 ? 1.0 - exp(-coal_rates[path][path][2*b]) : 1.0) *
                time_steps[b] *
                   (nbranches[b] + 1.0) / (nrecombs[b] + 1.0) / path_prob;
            matrix->G3[path][b] = (b<ntimes-2 ? 1.0 - exp(-coal_rates[path][path][2*b]) : 1.0) *
                time_steps[b] *
                (nbranches[b] / (nrecombs[b] + 1.0 + int(b < root_age_index)) /
                 path_prob);
            matrix->E[path][b] = 1.0 / ncoals_pop[pop][b];
        }
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
    calc_transition_probs(tree, model, states, lineages, &matrix);
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

                if (spr.recomb_node == child1) {
                    // right child is not disrupted
                    node2 = mapping[child2];
                    disrupt = true;
                } else if (spr.recomb_node == child2) {
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

            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1, path1);

        } else {   // node1 == recomb_node but time1 != recomb_time
            // SPR is on same branch as new chromosome
            if (spr.recomb_time > time1) {
                // we move with SPR subtree
                next_states[i] = state2_lookup.lookup(
                                 mapping[spr.recomb_node], time1, path1);

            } else {
                // SPR should not be able to coal back onto same branch
                // this would be a self cycle
                assert(spr.coal_node != node1);

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
                const int minage = (internal ?
                            tree->nodes[tree->nodes[tree->root].child[0]].age :
                                    0);
                const int path2 = model->consistent_path(path1,
                                      last_tree->nodes[spr.recomb_node].pop_path,
                                                         minage, spr.recomb_time, time2);
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
    node2 = (other == spr.coal_node ?
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

    int root_age;
    if (internal) {
        int subtree_root = last_tree->nodes[last_tree->root].child[0];
        int maintree_root = last_tree->nodes[last_tree->root].child[1];
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


    // probability of recombination rate and location
    int nbranches_k = lineages->nbranches[k] + int(k < a);
    int nrecombs_k = lineages->nrecombs[k] + int(k <= a) +
        int(k == a) - int(k >= max(root_age, a));
    double p = nbranches_k * model->time_steps[k] /
        (nrecombs_k * last_treelen_b) *
        (1.0 - exp(-max(model->rho * last_treelen, model->rho)));

    if (nrecombs_k <= 0 || nbranches_k <= 0) {
        printError("counts %d %d %e\n",
                   nrecombs_k, nbranches_k, p);
        assert(false);
    }

    return p;
}


//calc_recoal_sums sets *sums to the probability of not coalescing between
// spr.recomb_time and spr.coal_time, but does not consider the extra
// branch that may exist depending on where the current branch is being
// threaded.
// sums2[i] represents the probability of the broken branch not coalescing
//   at (half-interval) time i within the spr path so that it can be added
//   to the log probability of not coalescing with the rest of the tree in
// cases where this extra branch exists and is in the same path
// both sums only pertain to the population path taken by the Spr event
void calc_recoal_sums(const ArgModel *model, const LineageCounts *lineages,
                      const Spr &spr, const int recomb_parent_age,
                      const int recomb_node_path,
                      double *sums, double *sums2)
{
    // get info
    const int k = spr.recomb_time;
    const int j = spr.coal_time;

    // probability of not coalescing before time j
    double sum = 0.0;
    for (int m=2*k; m<2*j-1; m++) {
        int spr_pop = model->get_pop(spr.pop_path, (m+1)/2);
        int recomb_node_pop = model->get_pop(recomb_node_path, (m+1)/2);
        int nbranches_m = lineages->nbranches_pop[spr_pop][m/2] -
            int(m/2<recomb_parent_age &&
                spr_pop == recomb_node_pop);
        sum += (model->coal_time_steps[m] * nbranches_m
		/ (2.0 * model->popsizes[spr_pop][m]));
    }
    *sums = sum;

    sum = 0.0;
    sums2[2*k] = sum;
    for (int m=2*k; m<2*j-1; m++) {
        int spr_pop = model->get_pop(spr.pop_path, (m+1)/2);
        sum += model->coal_time_steps[m] / (2.0 * model->popsizes[spr_pop][m]);
        sums2[m+1] = sum;
    }
}


double calc_recoal(
    const LocalTree *last_tree, const ArgModel *model,
    const LineageCounts *lineages,
    const Spr &spr, const State &state1,
    const int recomb_parent_age,
    const int recomb_parent_pop,
    const bool internal=false)
{

    // get times
    const int a = state1.time;
    const int a_pop = model->get_pop(state1.pop_path, a);
    const int k = spr.recomb_time;
    const int j = spr.coal_time;
    const int coal_pop = model->get_pop(spr.pop_path, j);

    const int *nbranches = lineages->nbranches_pop[coal_pop];
    const int *ncoals = lineages->ncoals_pop[coal_pop];

    // probability of coalescing on chosen branch
    int nbranches_j = nbranches[j]
        - int(j < recomb_parent_age && recomb_parent_pop == coal_pop)
        + int(j < a && a_pop == coal_pop);
    int ncoals_j = ncoals[j]
        - int(j <= recomb_parent_age && recomb_parent_pop == coal_pop)
        - int(j == recomb_parent_age) + int(j <= a) + int(j == a);
    bool over = false;
    if (internal) {
        int subtree_root = last_tree->nodes[last_tree->root].child[0];
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
    double p = 1.0 / ncoals_j;

    // probability of coalescing in time interval j
    if (j < model->ntimes - 2) {
        double Z = 0.0;
        if (j>k) {
            int b1 = nbranches[j-1]
                - int(j-1 < recomb_parent_age && recomb_parent_pop == coal_pop)
                + int(j-1 < a && a_pop == coal_pop);
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
    if (ncoals_j <= 0 || nbranches_j <= 0) {
        printError("counts %d %d %e\n",
                   ncoals_j, nbranches_j, p);
        assert(false);
    }
    assert(!isnan(p) && p>0);

    return p;
}


double calc_recomb_recoal(
    const LocalTree *last_tree, const ArgModel *model,
    const LineageCounts *lineages,
    const Spr &spr, const State state1,
    const int recomb_parent_age, const int recomb_parent_pop,
    double last_treelen,
    const bool internal=false)
{
    // get times
    int a = state1.time;
    const int k = spr.recomb_time;
    const int j = spr.coal_time;

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
        int nbranches_m = lineages->nbranches[m/2]
            - int(m/2<recomb_parent_age && spr_pop == recomb_node_pop)
            + int(m/2 < a && spr_pop == a_pop);
        sum += model->coal_time_steps[m] * nbranches_m /
            (2.0 * model->popsizes[spr_pop][m]);
    }
    p *= exp(-sum);

    p *= calc_recoal(last_tree, model, lineages, spr, state1,
                     recomb_parent_age, recomb_parent_pop, internal);
    p *= model->path_prob(spr.pop_path, spr.recomb_time,
                          spr.coal_time);
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
    int recomb_parent_pop;
    int minage = internal ? tree->nodes[tree->nodes[tree->root].child[0]].age : 0;

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

            assert(spr.coal_node == last_tree->root);
            const int maintree_root = tree->nodes[tree->root].child[1];
            const int subtree_root = tree->nodes[tree->root].child[0];
            const int target_path =
                model->consistent_path(tree->nodes[subtree_root].pop_path,
                                       spr.pop_path,
                                       tree->nodes[subtree_root].age,
                                       spr.recomb_time,
                                       spr.coal_time);

            for (int j=0; j<nstates2; j++) {
                if (states2[j].node == maintree_root &&
                    states2[j].time == spr.coal_time &&
                    states2[j].pop_path == target_path) {
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
                    recomb_parent_pop = model->get_pop(states1[1].pop_path,
                                                       states1[i].time);
                } else {
                    int recomb_parent = last_tree->nodes[spr.recomb_node].parent;
                    recomb_parent_age = last_tree->nodes[recomb_parent].age;
                    recomb_parent_pop = model->get_pop(last_tree->nodes[recomb_parent].pop_path,
                                                       recomb_parent_age);
                }
                transmat_switch->determprob[i] = calc_recomb_recoal(
                    last_tree, model, lineages, spr,
                    states1[i], recomb_parent_age, recomb_parent_pop,
                    last_treelen, internal);
            }
            return;
        }
    }  // if internal


    // get deterministic transitions
    NodeStateLookup state2_lookup(states2, minage, model->pop_tree);
    get_deterministic_transitions(last_tree, tree, spr, mapping,
                                  states1, states2, state2_lookup,
                                  model, transmat_switch->determ,
                                  internal);

    // calculate terms for deterministic probabilities
    recomb_parent_age = last_tree->nodes[
        last_tree->nodes[spr.recomb_node].parent].age;
    double sums, sums2[model->ntimes*2+1];
    calc_recoal_sums(model, lineages, spr, recomb_parent_age,
                     last_tree->nodes[spr.recomb_node].pop_path,
                     &sums, sums2);
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
                recomb_parent_pop = model->get_pop(states1[i].pop_path,
                                                   recomb_parent_age);
                transmat_switch->determprob[i] =
                calc_recomb_recoal(last_tree, model, lineages, spr,
                                   states1[i], recomb_parent_age, recomb_parent_pop,
                                   last_treelen, internal);
            } else {
                int pop = model->get_pop(states1[i].pop_path, states1[i].time);
                int recoal_pos = pop * model->ntimes + states1[i].time;
                int recomb_parent = last_tree->nodes[spr.recomb_node].parent;
                recomb_parent_age = last_tree->nodes[recomb_parent].age;
                recomb_parent_pop = model->get_pop(last_tree->nodes[recomb_parent].pop_path,
                                                   recomb_parent_age);
                if (recoals[recoal_pos] < 0)
                    recoals[recoal_pos] = calc_recoal(last_tree, model, lineages,
                                                      spr, states1[i],
                                                      recomb_parent_age,
                                                      recomb_parent_pop);
                transmat_switch->determprob[i] =
                    calc_recomb(last_tree, model, lineages, spr, states1[i],
                                last_treelen, internal) *
                    exp(-sums
                        -sums2[max(min(2*spr.coal_time-1, 2*states1[i].time),
                                   2*spr.recomb_time)]) *
                    recoals[recoal_pos] *
                    model->path_prob(spr.pop_path, spr.recomb_time,
                                     spr.coal_time);
            }
        }
    }




    // find probabilitistic transition source states
    for (int i=0; i < nstates1; i++) {
        if (states1[i].node == spr.recomb_node &&
            states1[i].time == spr.recomb_time) {
            transmat_switch->recombsrc[i] = states1[i].pop_path;
            int recomb_next_states[2];
            get_recomb_transition_switch(model, tree, last_tree, spr, mapping,
                                         state2_lookup, recomb_next_states,
                                         states1[i].pop_path, minage);

            // stay case (recomb above)
            int j = recomb_next_states[0];
            if (j != -1) {
                recomb_parent_age = last_tree->nodes[
                    last_tree->nodes[spr.recomb_node].parent].age;
                recomb_parent_pop =
                    model->get_pop(last_tree->nodes[spr.recomb_node].pop_path,
                                   recomb_parent_age);
                transmat_switch->set(i, j,
                    calc_recomb_recoal(last_tree, model, lineages, spr,
                                       states1[i], recomb_parent_age,
                                       recomb_parent_pop,
                                       last_treelen, internal));
                assert(!isnan(transmat_switch->get(i, j)));
            }

            // escape case (recomb below)
            j = recomb_next_states[1];
            if (j != -1) {
                recomb_parent_age = states1[i].time;
                transmat_switch->set(i, j,
                    calc_recomb_recoal(last_tree, model, lineages, spr, states1[i],
                                       recomb_parent_age, recomb_parent_pop,
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
            int last_parent = last_tree->nodes[spr.recomb_node].parent;
            if (last_parent == node1) {
                // recomb breaks node1 branch, we need to use the other child
                const int *c = last_tree->nodes[last_parent].child;
                node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
            } else{
                node3 = mapping[node1];
            }

            int parent = tree->nodes[mapping[spr.recomb_node]].parent;
            assert(parent == tree->nodes[node3].parent);

            for (int j=0; j<nstates2; j++) {
                const int node2 = states2[j].node;
                const int time2 = states2[j].time;

                if (!((node2 == mapping[spr.recomb_node]
                       && time2 >= spr.recomb_time) ||
                      (node2 == node3 && time2 == time1) ||
                      (node2 == parent && time2 == time1)))
                    // not a probabilistic transition
                    continue;

                if (model->pop_tree != NULL) {
                    const int path2 =
                        model->consistent_path(states2[j].pop_path,
                                               tree->nodes[node2].pop_path,
                                               minage, states2[j].time,
                                               states1[i].time, false);
                    if (path2 == -1) continue;
                    if (! model->paths_equal(states1[i].pop_path,
                                             path2,
                                             minage, states1[i].time))
                        continue;
                }

                recomb_parent_age = last_tree->nodes[last_tree->nodes[spr.recomb_node].parent].age;
                recomb_parent_pop =
                    model->get_pop(last_tree->nodes[spr.recomb_node].pop_path,
                                   recomb_parent_age);
                Spr spr2 = spr;
                spr2.coal_time = time2;
                transmat_switch->set(i, j,
                    calc_recomb_recoal(last_tree, model, lineages, spr2, states1[i],
                                       recomb_parent_age, recomb_parent_pop,
                                       last_treelen, internal));
                assert(!isnan(transmat_switch->get(i, j)));
            }
        }
    }
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
        sum += model->coal_time_steps[m] * nbranches_pop[pop][m/2]
            / (2.0 * model->popsizes[pop][m]);
    }
    pop = model->get_pop(pop_path, b);
    double p = exp(-sum) / ncoals_pop[pop][b];

    // probability of coalescing in time interval b
    if (b < model->ntimes - 2) {
        double Z = 0.0;
        if (b>minage)
            Z = model->coal_time_steps[2*b-1] *
                nbranches_pop[pop][b-1]/
                ( 2.0 * model->popsizes[pop][2*b-1] );
        p *= 1.0 - exp(- (model->coal_time_steps[2*b] *
                          nbranches_pop[pop][b])/
                       (2.0 * model->popsizes[pop][2*b]) - Z);
    } else {
        // b = ntimes -1, guaranteed coalescence
    }

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
    calc_transition_probs(&tree, &model, states, &lineages, &transmat, true);

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


