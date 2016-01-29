//=============================================================================
// ArgHmm states


#ifndef ARGWEAVER_STATES_H
#define ARGWEAVER_STATES_H

#include "common.h"
#include "local_tree.h"

namespace argweaver {

class PopulationTree;

// A state in the ArgHmm
//
// Each state represents a node and time where coalescing is allowed
class State
{
public:
    State(int node=0, int time=0, int pop_path=0, double state_prior=-1) :
      node(node), time(time), pop_path(pop_path), state_prior(state_prior) {}

      // NOTE do not compare probabilities
    inline bool operator==(const State &other) const {
        return (node == other.node) && (time == other.time) &&
        (pop_path == other.pop_path);
    }

    void set(const int &_node, const int &_time, const int &_pop_path,
             const double &_state_prior) {
        node = _node;
        time = _time;
        pop_path = _pop_path;
        state_prior = _state_prior;
    }

    void set_null()
    {
        node = -1;
        time = -1;
        pop_path = -1;
        state_prior = -1;
    }

    bool is_null() const
    {
        return node == -1 && time == -1 && pop_path == -1;
    }

    int node;
    int time;
    int pop_path;
    double state_prior;
};

// A state space for a local block
typedef vector<State> States;


// This data structure provides a mapping from (node, time, path) tuples to
// the corresponding state index.
class NodeStateLookup
{
public:
    NodeStateLookup(const States &states, int nnodes, int npaths) :
        states(states),
        nstates(states.size()),
        nnodes(nnodes)
    {
        const int nstates = states.size();
        const int MAXTIME = 1000000;

        // allocate lookup arrays
        node_path_offset = new int[nnodes * npaths];
        state_lookup = new int[nstates];
        nstates_per_node_path = new int[nnodes * npaths];

        // count number of states per node/path and mintime per node/path
        int node_path_mintimes[nnodes * npaths];

        // initialize arrays
        for (int i=0; i<nnodes; i++) {
            nstates_per_node_path[i] = 0;
            node_path_mintimes[i] = MAXTIME;
        }

        for (int i=0; i<nstates; i++) {
            state_lookup[i] = -1;
            int idx = states[i].pop_path * nnodes + states[i].node;
            nstates_per_node_path[idx]++;
            node_path_mintimes[idx] = min(node_path_mintimes[idx],
                                          states[i].time);
        }

        // setup node_offsets
        int offset = 0;
        int idx=0;
        for (int i=0; i < npaths; i++) {
            for (int j=0; j<nnodes; j++) {
                node_path_offset[idx] = offset - node_path_mintimes[idx];
                offset += nstates_per_node_path[idx++];
            }
        }

        // set states
        for (int i=0; i<nstates; i++) {
            int j = node_path_offset[states[i].node + states[i].pop_path * nnodes]
                + states[i].time;
            assert(j >=0 && j<nstates);
            assert(state_lookup[j] == -1);
            state_lookup[j] = i;
        }
    }

    ~NodeStateLookup()
    {
        // clean up lookup arrays
        delete [] node_path_offset;
        delete [] state_lookup;
        delete [] nstates_per_node_path;
    }

    // Returns the state index for state (node, time)
    inline int lookup(int node, int time, int path) const {
        int idx = path * nnodes + node;
        if (nstates_per_node_path[idx] == 0)
            return -1;
        const int i = node_path_offset[idx] + time;
        if (i < 0 || i >= nstates)
            return -1;
        const int statei = state_lookup[i];
        if (statei == -1)
            return -1;
        if (states[statei].node != node ||
            states[statei].time != time ||
            states[statei].pop_path != path)
            return -1;
        return statei;
    }

protected:
    const States &states;
    int nstates;
    int nnodes;
    int *node_path_offset;
    int *state_lookup;
    int *nstates_per_node_path;
};



// A simple representation of a state, useful for passing from python
// Not updated for multiple populations
typedef int intstate[2];


// Converts integer-based states to State class
void make_states(intstate *istates, int nstates, States &states);

// Converts state class represent to integer-based
void make_intstates(States states, intstate *istates);

void get_coal_states(const LocalTree *tree, int ntimes, States &states,
                     bool internal=false, PopulationTree *pop_tree=NULL,
                     int start_pop=-1);
void get_coal_states_external(const LocalTree *tree, int ntimes, States &states,
                              int minage=0, PopulationTree *pop_tree=NULL,
                              int start_pop=-1);
void get_coal_states_internal(const LocalTree *tree, int ntimes,
                              States &states, int minage=0,
                              PopulationTree *pop_tree=NULL,
                              int start_pop=-1);


// NOTE: the three get_num_coal_states functions below are not quite accurate
// and may count some states which are not possible and not included in the
// actual states used. However these functions are currently only used for
// logging purposes.
int get_num_coal_states(const LocalTree *tree, int ntimes, bool internal=false,
                        PopulationTree *pop_tree=NULL, int start_pop=-1);
int get_num_coal_states_external(const LocalTree *tree, int ntimes,
                                 int minage=0, PopulationTree *pop_tree=NULL,
                                 int start_pop=-1);
int get_num_coal_states_internal(const LocalTree *tree, int ntimes,
                                 int minage=0, PopulationTree *pop_tree=NULL,
                                 int start_pop=-1);


class StatesModel
{
public:
    StatesModel(int ntimes=0, bool internal=false, int minage=0,
             PopulationTree *pop_tree = NULL) :
        ntimes(ntimes),
        internal(internal),
        minage(minage),
        pop_tree(pop_tree)

    {}

        void set(int _ntimes, bool _internal, int _minage,
                 PopulationTree *_pop_tree = NULL) {
        ntimes = _ntimes;
        internal = _internal;
        minage = _minage;
        pop_tree = _pop_tree;
    }

    void get_coal_states(const LocalTree *tree, States &states,
                         int start_pop=-1) const {
        if (!internal)
            get_coal_states_external(tree, ntimes, states, minage, pop_tree,
                                     start_pop);
        else
            get_coal_states_internal(tree, ntimes, states, minage, pop_tree,
                                     start_pop);
    }

    int ntimes;
    bool internal;
    int minage;
    PopulationTree *pop_tree;
};


} // namespace argweaver

#endif // ARGWEAVER_STATES_H
