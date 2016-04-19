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
    State(int node=0, int time=0, int pop_path=0) :
      node(node), time(time), pop_path(pop_path) {}

    inline bool operator==(const State &other) const {
        return (node == other.node) && (time == other.time) &&
        (pop_path == other.pop_path);
    }

    void set(const int &_node, const int &_time, const int &_pop_path) {
        node = _node;
        time = _time;
        pop_path = _pop_path;
    }

    void set_null()
    {
        node = -1;
        time = -1;
        pop_path = -1;
    }

    bool is_null() const
    {
        return node == -1 && time == -1;
    }

    int node;
    int time;
    int pop_path;
};

// A state space for a local block
typedef vector<State> States;


// This data structure provides a mapping from (node, time, path) tuples to
// the corresponding state index.
class NodeStateLookup
{
public:
  NodeStateLookup(const States &states, int minage,
                  const PopulationTree *pop_tree);
  ~NodeStateLookup();

    // Returns the state index for state (node, time)
  inline int lookup(int node, int time, int path=0) {
      if (node < 0 || node >= nnode ||
          time < mintime || time > maxtime ||
          path < 0 || path >= npath) return -1;
      return lookup_table[path*nnode*ntime + node*ntime + time - mintime];
  }

protected:
  const States &states;
  int nnode;
  int mintime;
  int maxtime;
  int ntime;
  int npath;
  int *lookup_table;
};


// A simple representation of a state, useful for passing from python
// Not updated for multiple populations
typedef int intstate[2];


// Converts integer-based states to State class
void make_states(intstate *istates, int nstates, States &states);

// Converts state class represent to integer-based
void make_intstates(States states, intstate *istates);

void get_coal_states(const LocalTree *tree, int ntimes, States &states,
                     bool internal=false, const PopulationTree *pop_tree=NULL,
                     int start_pop=0);
void get_coal_states_external(const LocalTree *tree, int ntimes, States &states,
                              int minage=0, const PopulationTree *pop_tree=NULL,
                              int start_pop=0);
void get_coal_states_internal(const LocalTree *tree, int ntimes,
                              States &states, int minage=0,
                              const PopulationTree *pop_tree=NULL);


// NOTE: the three get_num_coal_states functions below are not quite accurate
// and may count some states which are not possible and not included in the
// actual states used. However these functions are currently only used for
// logging purposes.
int get_num_coal_states(const LocalTree *tree, int ntimes, bool internal=false);
int get_num_coal_states_external(const LocalTree *tree, int ntimes);
int get_num_coal_states_internal(const LocalTree *tree, int ntimes,
                                 int minage=0);


class StatesModel
{
public:
    StatesModel(int ntimes=0, bool internal=false, int minage=0,
                PopulationTree *pop_tree = NULL, int start_pop=0) :
        ntimes(ntimes),
        internal(internal),
        minage(minage),
        start_pop(start_pop),
        pop_tree(pop_tree)

    {}

    void set(int _ntimes, bool _internal, int _minage,
                 PopulationTree *_pop_tree = NULL, int _start_pop = 0) {
        ntimes = _ntimes;
        internal = _internal;
        minage = _minage;
        pop_tree = _pop_tree;
        start_pop = _start_pop;
    }

    void set_start_pop(int _start_pop, const PopulationTree *_pop_tree) {
       start_pop = _start_pop;
       pop_tree = _pop_tree;
    }


    void get_coal_states(const LocalTree *tree, States &states) const {
        if (!internal)
            get_coal_states_external(tree, ntimes, states, minage, pop_tree,
                                     start_pop);
        else
            get_coal_states_internal(tree, ntimes, states, minage, pop_tree);

    }

    int ntimes;
    bool internal;
    int minage;
    int start_pop;
    const PopulationTree *pop_tree;
};


} // namespace argweaver

#endif // ARGWEAVER_STATES_H
