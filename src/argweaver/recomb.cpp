
#include "local_tree.h"
#include "matrices.h"

namespace argweaver {

using namespace std;


//=============================================================================
// Sample recombinations


/* relative probability of a recombination on a fully realized local tree
   ( assumes threading already completed )
 */
double recomb_prob_smcPrime_unnormalized(const ArgModel *model,
                                         const LocalTree *tree,
                                         const LineageCounts &lineages,
                                         const Spr &recomb) {
    const int recomb_time = recomb.recomb_time;
    const int recomb_path = recomb.pop_path;
    const int coal_time = recomb.coal_time;
    if (model->get_pop(recomb_path, recomb_time) !=
        model->get_pop(tree->nodes[recomb.recomb_node].pop_path, recomb_time))
        return 0.0;
    if (model->get_pop(recomb_path, coal_time) !=
        model->get_pop(tree->nodes[recomb.coal_node].pop_path, coal_time))
        return 0.0;
    assert(recomb_time <= coal_time);

    double nocoal_rate = 0.0;
    double coal_rate = 0.0;
    int pop = -1;
    for (int i=recomb_time*2; i <= 2*coal_time; i++) {
        pop = model->get_pop(recomb_path, (i+1)/2);
        double rate = model->coal_time_steps[i] * lineages.nbranches_pop[pop][i]
            / (2.0 * model->popsizes[pop][i]);
        if (i < 2*coal_time - 1)
            nocoal_rate += rate;
        else coal_rate += rate;
    }
    return model->time_steps[recomb_time] * lineages.nbranches[recomb_time]
        / lineages.nrecombs[recomb_time]
        * model->path_prob(recomb_path, recomb_time, coal_time)
        * exp(-nocoal_rate) * (1.0 - exp(-coal_rate))
        / lineages.ncoals_pop[pop][coal_time];
}


// assumes consistency between tree, last_state, state, recomb
double recomb_prob_unnormalized(const ArgModel *model, const LocalTree *tree,
                                const LineageCounts &lineages,
                                const State &last_state,
                                const State &state,
                                const NodePointPath &recomb, bool internal)
{
    const int k = recomb.time;
    const int j = state.time;

    int root_time;
    int recomb_parent_age;
    int recomb_node_path;
    int minage=0;

    if (internal) {
        int subtree_root = tree->nodes[tree->root].child[0];
        int maintree_root = tree->nodes[tree->root].child[1];
        minage = tree->nodes[subtree_root].age;
        root_time = max(tree->nodes[maintree_root].age, last_state.time);
        recomb_parent_age = (recomb.node == subtree_root ||
                             tree->nodes[recomb.node].parent == -1 ||
                             recomb.node == last_state.node) ?
            last_state.time :
            tree->nodes[tree->nodes[recomb.node].parent].age;
        if (recomb.node == subtree_root) {
            recomb_node_path = last_state.pop_path;
        } else if (recomb.node == last_state.node) {
            recomb_node_path = tree->nodes[last_state.node].pop_path;
        } else assert(0);
    } else {
        root_time = max(tree->nodes[tree->root].age, last_state.time);
        recomb_parent_age = (recomb.node == -1 ||
                             tree->nodes[recomb.node].parent == -1 ||
                             recomb.node == last_state.node) ?
            last_state.time :
            tree->nodes[tree->nodes[recomb.node].parent].age;
        if (recomb.node == -1) {
            recomb_node_path = last_state.pop_path;
        } else if (recomb.node == last_state.node) {
            recomb_node_path = tree->nodes[last_state.node].pop_path;
        } else assert(0);
    }

    int nbranches_k = lineages.nbranches[k]
        + int(k < last_state.time && k >= minage);
    int nrecombs_k = lineages.nrecombs[k]
        + int(k <= last_state.time && k >= minage)
        + int(k == last_state.time && k >= minage)
        - int(k == root_time);

    double precomb = nbranches_k * model->time_steps[k] / nrecombs_k;

    // probability of not coalescing before time j-1
    double coal_sum=0.0;
    double nocoal_sum = 0.0;
    for (int m=2*k; m<=2*j; m++) {
        int pop_t = (m+1)/2;   // round up to get population assignment
        int pop = model->get_pop(recomb.path, pop_t);
        int coal_t = m/2; // round down to determine which time interval branch is in
        int nbranches_m = lineages.nbranches_pop[pop][m]
            + int( coal_t < last_state.time &&
                   coal_t >= minage &&
                   model->get_pop(last_state.pop_path, pop_t) == pop)
            - int( (! model->smc_prime) &&
                   coal_t < recomb_parent_age &&
                   model->get_pop(recomb_node_path, pop_t)==pop);
        double rate = (model->coal_time_steps[m] * nbranches_m
                       / (2.0 * model->popsizes[pop][m]));
        if (m >= 2*j-1)
            coal_sum += rate;
        else nocoal_sum += rate;
    }
    if (coal_sum == 0.0) return 0.0;
    double pcoal = ( j >= model->ntimes - 2 ? 1.0 : (1.0 - exp(-coal_sum)));

    // probability of recoalescing on a choosen branch
    int pop = model->get_pop(recomb.path, j);
    int recomb_parent_pop = model->get_pop(recomb_node_path, j);
    int last_state_pop = model->get_pop(last_state.pop_path, j);
    int ncoals_j = lineages.ncoals_pop[pop][j]
        - int( (!model->smc_prime) && j <= recomb_parent_age && recomb_parent_pop==pop)
        - int( (!model->smc_prime) && j == recomb_parent_age && recomb_parent_pop==pop)
        + int(j <= last_state.time && j >= minage && pop==last_state_pop)
        + int(j == last_state.time && j >= minage && pop==last_state_pop);

    // Also return 0 (or maybe assert error?) if recomb path does not end in correct pop
    /*    printf("%i recomb_prob_unnormalized %.1e %i %.1e %.1e %.1e\n",
          count++, pcoal, ncoals_j, precomb, exp(-nocoal_sum), model->path_prob(recomb.path, k, j));*/

    double p = pcoal / ncoals_j * precomb * exp(- nocoal_sum)
        * model->path_prob(recomb.path, k, j);
    return p;
}


// Returns the possible recombination events that are compatiable with
// the transition (last_state -> state).
void get_possible_recomb(const ArgModel *model, const LocalTree *tree,
                         const State last_state, const State state,
                         bool internal, vector<NodePointPath> &candidates)
{
    // represents the branch above the new node in the tree.
    const int new_node = -1;
    int end_time = min(state.time, last_state.time);

    // note: cannot normally compare pop_paths directly but OK in following
    // line because get_coal_states() assures that the paths are consistent
    // with branch that is being coalesced to, so when the node is the same
    // the paths should be comparable
    if (state.node == last_state.node && state.pop_path == last_state.pop_path) {
        // y = v, k in [0, min(timei, last_timei)]
        // y = node, k in Sr(node)
        for (int k=tree->nodes[state.node].age; k<=end_time; k++)
            candidates.push_back(NodePointPath(state.node, k,
                                               tree->nodes[state.node].pop_path));
    }

    if (internal) {
        const int subtree_root = tree->nodes[tree->root].child[0];
        const int subtree_root_age = tree->nodes[subtree_root].age;
        for (int k=subtree_root_age; k<=end_time; k++) {
            if (model->get_pop(state.pop_path, k) !=
                model->get_pop(last_state.pop_path, k)) break;
            candidates.push_back(NodePointPath(subtree_root, k, state.pop_path));
        }
    } else {
        for (int k=0; k<=end_time; k++) {
            if (model->get_pop(state.pop_path, k) !=
                model->get_pop(last_state.pop_path, k)) break;
            candidates.push_back(NodePointPath(new_node, k, state.pop_path));
        }
    }
}


// if using SMC' model, this will not sample invisible recombinations.
// those can be sampled later with sample_invisible_recombinations
void sample_recombinations(
    const LocalTrees *trees, const ArgModel *model,
    ArgHmmMatrixIter *matrix_iter,
    int *thread_path, vector<int> &recomb_pos, vector<NodePointPath> &recombs,
    bool internal)
{
    States states;
    LineageCounts lineages(model->ntimes, model->num_pops());
    vector <NodePointPath> candidates;
    vector <double> probs;

    // loop through local blocks
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {

        // get local block information
        ArgHmmMatrices &matrices = matrix_iter->ref_matrices();
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        lineages.count(tree, model->pop_tree, internal);
        matrices.states_model.get_coal_states(tree, states);
        int next_recomb = -1;

        // don't sample recombination if there is no state space
        if (internal && states.size() == 0)
            continue;

        int start = matrix_iter->get_block_start();
        int end = matrix_iter->get_block_end();
        if (matrices.transmat_switch || start == trees->start_coord) {
            // don't allow new recomb at start if we are switching blocks
            start++;
        }

        // loop through positions in block
        for (int i=start; i<end; i++) {

            if (thread_path[i] == thread_path[i-1]) {
                // no change in state, recombination is optional
                // under SMC prime invisible recombs are sampled later
                if (model->smc_prime) continue;

                if (i > next_recomb) {
                    // sample the next recomb pos
                    int last_state = thread_path[i-1];
                    TransMatrix *m = matrices.transmat;
                    int a = states[last_state].time;
                    double self_trans = m->get(
                        tree, states, last_state, last_state);
                    double rate = 1.0 - (m->norecombs[a] / self_trans);

                    // NOTE: the min prevents large floats from overflowing
                    // when cast to int
                    next_recomb = int(min(double(end), i + expovariate(rate)));
                }

                if (i < next_recomb)
                    continue;
            }


            // sample recombination
            next_recomb = -1;
            State state = states[thread_path[i]];
            State last_state = states[thread_path[i-1]];

            // there must be a recombination
            // either because state changed or we choose to recombine
            // find candidates
            candidates.clear();
            get_possible_recomb(model, tree, last_state, state, internal,
                                candidates);

            // compute probability of each candidate
            probs.clear();
            for (vector<NodePointPath>::iterator it=candidates.begin();
                 it != candidates.end(); ++it) {
                probs.push_back(recomb_prob_unnormalized(
                    model, tree, lineages, last_state, state, *it, internal));
            }

            // sample recombination
            recomb_pos.push_back(i);
            int r = sample(&probs[0], probs.size());
            recombs.push_back(candidates[r]);
            /*            printf("%i\t%i\n", recomb_pos[recomb_pos.size()-1],
                          recombs[recombs.size()-1].time);*/
            assert(recombs[recombs.size()-1].time <= min(state.time,
                                                         last_state.time));
        }
    }
}


// this is meant to be called on a full ARG (after threading is complete)
// assumes no invisible recombinations currently exist; resamples all of them
// there can be more than one per position
void sample_invisible_recombinations(const ArgModel *model, LocalTrees *trees,
                                     vector<int> &recomb_pos,
                                     vector<Spr> &recombs) {
    if (!model->smc_prime) {
        fprintf(stderr, "sample_invisible_recombinations should only be called for smc Prime model\n");
        exit(0);
    }
    LineageCounts lineages(model->ntimes, model->num_pops());
    recomb_pos.clear();
    recombs.clear();
    int end = trees->start_coord;
    int idx=0;
    vector<Spr> possible_recombs;
    vector<double> recomb_probs;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        int start = end;
        end += it->blocklen;
        lineages.count(tree, model->pop_tree, false);
        double treelen = get_treelen(tree, model->times, model->ntimes, true);
        double rho = model->get_local_rho(start, &idx);
        double d_term = (1.0 - exp(-rho * treelen)) / treelen;
        possible_recombs.clear();
        recomb_probs.clear();
        double total_prob = 0.0;
        for (int node=0; node < tree->nnodes; node++) {
            if (node == tree->root) continue;
            int pop_path = tree->nodes[node].pop_path;
            int parent = tree->nodes[node].parent;
            int minage = tree->nodes[node].age;
            int maxage = tree->nodes[parent].age;
            for (int k=minage; k <= maxage; k++) {
                for (int j=k; j <= maxage; j++) {
                    possible_recombs.push_back(Spr(node, k, node, j, pop_path));
                    double prob = recomb_prob_smcPrime_unnormalized(model,
                                        tree, lineages, possible_recombs.back())
                        * d_term;
                    recomb_probs.push_back(prob);
                    total_prob += prob;
                }
                // can also recombine onto parent or sister at coalescence time
                possible_recombs.push_back(Spr(node, k, parent,
                                               maxage, pop_path));
                double prob = recomb_prob_smcPrime_unnormalized(model,
                            tree, lineages, possible_recombs.back()) * d_term;
                recomb_probs.push_back(prob);
                total_prob += prob;
                int sib = tree->nodes[parent].child[0];
                if (sib == node)
                    sib = tree->nodes[parent].child[1];
                possible_recombs.push_back(Spr(node, k, sib,
                                               maxage, pop_path));
                prob = recomb_prob_smcPrime_unnormalized(model,
                            tree, lineages, possible_recombs.back()) * d_term;
                recomb_probs.push_back(prob);
                total_prob += prob;
            }
        }
        assert(total_prob >= 0.0 && total_prob <= 1.0);
        for (int i=start; i < end; i++) {
            int next_recomb = min( i + int(expovariate(total_prob)), end);
            if (next_recomb < end) {
                double r = frand(total_prob);
                total_prob = 0.0;
                bool found=false;
                for (unsigned int j=0; j < possible_recombs.size(); j++) {
                    total_prob += recomb_probs[j];
                    if (r <= total_prob) {
                        recombs.push_back(possible_recombs[j]);
                        found=true;
                        break;
                    }
                }
                assert(found);
                recomb_pos.push_back(next_recomb);
            }
            i = next_recomb;
        }
    }
}


} // namespace argweaver
