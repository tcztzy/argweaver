

#include "states.h"
#include "pop_model.h"

namespace argweaver {

NodeStateLookup::NodeStateLookup(const States &states, int minage,
                                 const PopulationTree *pop_tree=NULL) :
    states(states) {
    npath = ( pop_tree == NULL ? 1 : pop_tree->num_pop_paths() );
    nnode=0;
    mintime=-1;
    maxtime=-1;

    for (unsigned int i=0; i < states.size(); i++) {
        if (states[i].node >= nnode)
            nnode = states[i].node+1;
        if (i==0 || states[i].time < mintime) mintime = states[i].time;
        if (i==0 || states[i].time > maxtime) maxtime = states[i].time;
    }
    ntime = maxtime - mintime + 1;
    table_size = npath * nnode * ntime;
    lookup_table = new int[table_size];
    for (int i=0; i < table_size; i++) lookup_table[i] = -1;

    for (unsigned int i=0; i < states.size(); i++) {
        int t = states[i].time;
        int node = states[i].node;
        if (pop_tree == NULL) {
            lookup_table[node*ntime + t - mintime] = i;
        } else {
            const set<int> *paths =
                pop_tree->get_equivalent_paths(states[i].pop_path, minage, t);
            for (set<int>::iterator it=paths->begin(); it != paths->end();
                 it++) {
                int idx = (*it)*nnode*ntime + node*ntime + t - mintime;
                assert(lookup_table[idx] == -1);
                lookup_table[idx] = i;
            }
        }
    }
}

NodeStateLookup::~NodeStateLookup() {
    delete [] lookup_table;
}


// Converts integer-based states to State class
// Not updated for multiple populations but only used for python
void make_states(intstate *istates, int nstates, States &states) {
    states.clear();
    for (int i=0; i<nstates; i++)
        states.push_back(State(istates[i][0], istates[i][1]));
}


// Converts state class represent to integer-based
// Not updated for multiple populations but only used for python
void make_intstates(States states, intstate *istates)
{
    const int nstates = states.size();
    for (int i=0; i<nstates; i++) {
        istates[i][0] = states[i].node;
        istates[i][1] = states[i].time;
    }
}



void get_coal_states(const LocalTree *tree, int ntimes, States &states,
                     bool internal, const PopulationTree *pop_tree,
                     int start_pop)
{
    if (internal)
        get_coal_states_internal(tree, ntimes, states, 0, pop_tree);
    else
        get_coal_states_external(tree, ntimes, states, 0, pop_tree, start_pop);
}


// Only used for logging purposes; not necessarily accurate
int get_num_coal_states(const LocalTree *tree, int ntimes, bool internal)
{
    if (internal)
        return get_num_coal_states_internal(tree, ntimes);
    else
        return get_num_coal_states_external(tree, ntimes);
}


// Returns the possible coalescing states for a tree
//
// NOTE: Do not allow coalescing at top time
// states for the same branch are clustered together and ages are given
// in increasing order
void get_coal_states_external(const LocalTree *tree, int ntimes, States &states,
			      int minage, const PopulationTree *pop_tree,
                              int start_pop)
{
    states.clear();
    const LocalNode *nodes = tree->nodes;

    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = max(minage, nodes[i].age);
        const int parent = nodes[i].parent;

        int max_time;
        if (parent == -1)
            max_time = ntimes - 2;
        else max_time = nodes[parent].age;
        if (pop_tree == NULL) {
            for ( ; time <= max_time; time++)
                states.push_back(State(i, time));
        } else {
            for ( ; time <= max_time; time++) {
                int target_path = pop_tree->path_to_root(nodes, i, time);
                int end_pop = pop_tree->path_pop(target_path, time);
                // loop over all unique paths
                for (int p=0;
                     p < pop_tree->num_paths(minage, start_pop, time, end_pop);
                     p++) {
                    double path_prob =
                        pop_tree->subpath_prob(minage, start_pop, time, end_pop, p);
                    if (path_prob > 0.0 &&
                        (pop_tree->max_migrations < 0 ||
                         (pop_tree->subpath_num_mig(minage, start_pop, time, end_pop, p)
                          <= pop_tree->max_migrations))) {
                        int path1 =
                            pop_tree->unique_path(minage, start_pop, time,
                                                  end_pop, p);
                        // use a path that is consistent with the branch
                        // we are coalescing to
                        int path =
                            pop_tree->consistent_path(path1, target_path,
                                                      minage, time, -1);
                        states.push_back(State(i, time, path));
                    }
                }
            }
        }
    }
}

// Returns the number of possible coalescing states for a tree
// NOTE: is not accurate for multiple populations; but currently only used
// for logging purposes
int get_num_coal_states_external(const LocalTree *tree, int ntimes)
{
    int nstates = 0;
    const LocalNode *nodes = tree->nodes;

    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = nodes[i].age;
        const int parent = nodes[i].parent;

        if (parent == -1) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                nstates++;
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                nstates++;
        }
    }

    return nstates;
}


// Returns the possible coalescing states for a tree that is in
// subtree-maintree format for internal branch resampling
//
// NOTE: Do not allow coalescing at top time
// states for the same branch are clustered together and ages are given
// in increasing order
void get_coal_states_internal(const LocalTree *tree, int ntimes,
                              States &states, int minage,
                              const PopulationTree *pop_tree)
{
    states.clear();
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;

    if (nodes[tree->root].age < ntimes) {
        // we have a fully specified tree, there are no states
        return;
    }

    int subtree_root = nodes[tree->root].child[0];
    minage = max(minage, nodes[subtree_root].age);

    // ignore root and whole subtree
    bool ignore[nnodes];
    fill(ignore, ignore + nnodes, false);
    ignore[tree->root] = true;

    // recurse down subtree
    int stack[nnodes];
    int stacki = 0;
    stack[stacki++] = subtree_root;
    while (stacki > 0) {
        // pop stack
        int node = stack[--stacki];
        ignore[node] = true;

        // push children on stack
        if (!nodes[node].is_leaf()) {
            stack[stacki++] = nodes[node].child[0];
            stack[stacki++] = nodes[node].child[1];
        }
    }
    int start_pop = (pop_tree == NULL ? 0 :
                     pop_tree->get_pop(nodes[subtree_root].pop_path, minage) );

    // iterate over the branches of the tree
    for (int i=0; i<nnodes; i++) {
        int time = max(nodes[i].age, minage);
        const int parent = nodes[i].parent;

        // skip subtree and root node
        if (ignore[i])
            continue;

        int max_time = nodes[parent].age;
        if (parent == tree->root) {
            // no parent, allow coalescing up basal branch until ntimes-2
            max_time = ntimes - 2;
        }

        if (pop_tree == NULL) {
            for (; time<=max_time; time++) {
                states.push_back(State(i, time));
            }
        } else {
            assert(time <= nodes[tree->root].age);
            for (; time<=max_time; time++) {
                int target_path = pop_tree->path_to_root(nodes, i, time);
                int end_pop = pop_tree->path_pop(target_path, time);
                // loop over all unique paths
                for (int p=0;
                     p < pop_tree->num_paths(minage, start_pop, time, end_pop);
                     p++) {
                    double path_prob =
                        pop_tree->subpath_prob(minage, start_pop, time, end_pop, p);
                    if (path_prob > 0.0 &&
                        (pop_tree->max_migrations < 0 ||
                         (pop_tree->subpath_num_mig(minage, start_pop, time, end_pop, p)
                          <= pop_tree->max_migrations))) {
                        int path1 =
                            pop_tree->unique_path(minage, start_pop, time,
                                                  end_pop, p);
                       // use a path that is consistent with the branch
                       // we are coalescing to
                        int path =
                            pop_tree->consistent_path(path1, nodes[i].pop_path,
                                                      minage, time, -1);
                        states.push_back(State(i, time, path));
                    }
                }
            }
        }
    }
}


// Returns the number of possible coalescing states for a tree
// NOTE: is not accurate, does not exclude states below the subtree and
// has not been updated for multiple populations.
// Currently it is only used for logging purposes
// TODO: update this?
int get_num_coal_states_internal(const LocalTree *tree, int ntimes, int minage)
{
    int nstates = 0;
    const LocalNode *nodes = tree->nodes;
    int subtree_root = nodes[tree->root].child[0];
    minage = max(minage, nodes[subtree_root].age);

    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = max(nodes[i].age, minage);
        const int parent = nodes[i].parent;

        // skip subtree root branch and root node
        if (i == subtree_root || i == tree->root)
            continue;

        if (parent == tree->root) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                nstates++;
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                nstates++;
        }
    }

    return nstates;
}


int find_state(const States &states, const State target, const ArgModel *model,
               int minage) {
    for (unsigned int i=0; i < states.size(); i++) {
        if (states[i].node != target.node) continue;
        if (states[i].time != target.time) continue;
        if (model->pop_tree == NULL ||
            states[i].pop_path == target.pop_path) return i;
        if (model->paths_equal(states[i].pop_path, target.pop_path,
                               minage, target.time))
            return i;
    }
    return -1;
}



//=============================================================================
// C interface

extern "C" {


void arghmm_get_nstates(LocalTrees *trees, int ntimes, bool internal,
                        int *nstates)
{
    States states;

    // iterate over local trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); it++) {
        LocalTree *tree = it->tree;

        get_coal_states(tree, ntimes, states, internal);
        int n = states.size();
        for (int j=0; j<it->blocklen; j++)
            nstates[i+j] = n;
        i += it->blocklen;
    }
}


// Returns state-spaces, useful for calling from python
intstate **get_state_spaces(LocalTrees *trees, int ntimes, bool internal)
{
    States states;

    // allocate state space
    intstate **all_states = new intstate* [trees->get_num_trees()];

    // iterate over local trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); it++) {
        LocalTree *tree = it->tree;
        get_coal_states(tree, ntimes, states);
        int nstates = states.size();
        all_states[i] = new intstate [nstates];

        for (int j=0; j<nstates; j++) {
            all_states[i][j][0] = states[j].node;
            all_states[i][j][1] = states[j].time;
        }
        i++;
    }

    return all_states;
}


// Deallocate state space memory
void delete_state_spaces(intstate **all_states, int ntrees)
{
    for (int i=0; i<ntrees; i++)
        delete [] all_states[i];
    delete [] all_states;
}


} // extern "C"

} // namespace argweaver
