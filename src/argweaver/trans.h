//=============================================================================
// transitions

#ifndef ARGWEAVER_TRANS_H
#define ARGWEAVER_TRANS_H

#include "common.h"
#include "local_tree.h"
#include "model.h"
#include "states.h"
#include "MultiArray.h"

namespace argweaver {


class PopulationTree;

// A compressed representation of the transition matrix.
//
// This transition matrix is used in the chromosome threading HMM within
// one non-recombining alignment block.
class TransMatrix
{
public:
    TransMatrix(const ArgModel *model, int nstates) :
        nstates(nstates),
        internal(false),
        smc_prime(false)
    {
        initialize(model, nstates);
    }

    ~TransMatrix()
    {
        delete C1_prime;
        delete Q1_prime;
        if (!smc_prime) {
            for (int i=0; i < npaths; i++) {
                delete [] lnB[i];
                delete [] lnE2[i];
                delete [] lnNegG1[i];
            }
            delete [] E;
            delete [] lnB;
            delete [] lnE2;
            delete [] lnNegG1;
            delete [] G2;
            delete [] G3;
        } else {
            delete B0_prime;
            delete B1_prime;
            delete B2_prime;
            delete C0_prime;
            delete Q0_prime;
            delete E0_prime;
            delete E1_prime;
            delete E2_prime;
            delete F0_prime;
            delete F1_prime;
            delete F2_prime;
            delete G0_prime;
            delete G1_prime;
            delete G2_prime;
            delete L0_prime;
            delete L1_prime;
            delete L2_prime;
            delete K0_prime;
            delete K1_prime;
            delete K2_prime;
            delete RK0_prime;
            delete RK2_prime;
            /*            delete RK0diff;
                          delete RK2diff;*/
            delete [] self_recomb;
        }
        delete [] path_prob;
        delete [] data_alloc;
    }

    // allocate space for transition matrix
    // and initialize paths_equal matrix
    void initialize(const ArgModel *model, int nstates);

    // Probability of transition from state i to state j.
    inline double get(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        int minage = 0;
        if (internal) {
            if (nstates == 0)
                return 1.0;
            const int subtree_root = tree->nodes[tree->root].child[0];
            const int subtree_age = tree->nodes[subtree_root].age;
            minage = subtree_age;
        }

        const int node1 = states[i].node;
        const int a = states[i].time;
        const int a_path = states[i].pop_path;
        const int node2 = states[j].node;
        const int b = states[j].time;
        const int b_path = states[j].pop_path;
        const int c = tree->nodes[node2].age;
        const int c_path = tree->nodes[node2].pop_path;

        return get_time(a, b, c, a_path, b_path, c_path,
                        minage, node1 == node2, i);
    }

    // Returns the probability of transition from state1 with time 'a'
    // to state2 with time 'b'.  The probability also depends on whether
    // the node changes between states ('same_node') or whether there is
    // a minimum age ('minage') allowed for the state.
    // last argument "state_a" is only used in same_node case when smc_prime
    // is turned on- it is necessary for calculating probability of an
    // invisible recomb on the tree
    double get_time(int a, int b, int c,
                    int path_a, int path_b, int path_c,
                    int minage, bool same_node, int state_a=-1) const;

    // Log probability of transition from state i to state j.
    inline double get_log(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        return log(get(tree, states, i, j));
    }
    void assert_transmat(const LocalTree *tree,
                         const ArgModel *model,
                         const States &states,
                         bool internal, int minage) const;


    void calc_transition_probs_smcPrime(const LocalTree *tree, const ArgModel *model,
                                        const States &states,
                                        const LineageCounts *lineages,
                                        bool internal0=false,
                                        int minage0=0);
    void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
                               const States &states,
                               const LineageCounts *lineages,
                               bool internal0=false,
                               int minage0=0);


    int ntimes;
    int nstates;
    int npaths;
    bool internal;  // If true, this matrix is for threading an internal branch.
    int minage;     // Minimum age of a state we can consider (due to threading
                    // an internal branch).

    bool smc_prime;

    double *data_alloc;

    // Intermediate terms in calculating entries in the full transition matrix

    // these are used only for SMC' model
    MultiArray *C0_prime;
    MultiArray *Q0_prime;
    MultiArray *C1_prime;
    MultiArray *Q1_prime;

    MultiArray *B0_prime;
    MultiArray *B1_prime;
    MultiArray *B2_prime;
    MultiArray *G0_prime;
    MultiArray *G1_prime;
    MultiArray *G2_prime;
    MultiArray *E0_prime;
    MultiArray *E1_prime;
    MultiArray *E2_prime;
    MultiArray *F0_prime;
    MultiArray *F1_prime;
    MultiArray *F2_prime;
    MultiArray *L0_prime;
    MultiArray *L1_prime;
    MultiArray *L2_prime;
    MultiArray *K0_prime;
    MultiArray *K1_prime;
    MultiArray *K2_prime;
    MultiArray *RK0_prime;
    MultiArray *RK2_prime;
    /*    MultiArray *RK0diff;
          MultiArray *RK2diff;*/
    double *self_recomb;

    // these are used for both SMC and SMC'
    double *D;
    double **path_prob;
    double *norecombs;
    PopulationTree *pop_tree;

    // these used for SMC only
    double **E;
    double ***lnB;
    double ***lnE2;
    double ***lnNegG1;
    double **G2;
    double **G3;

 private:
    double get_l_term(int d, int path_d, int a, int path_a) const;
    double get_k_term(int d, int path_d, int a, int path_a) const;
    double get_b_term(int d, int path_d, int a, int path_a) const;
    double get_b_term_old(int d, int path_d, int a, int path_a) const;
    double get_rk_term(int d, int path_d, int a, int path_a) const;
    double selfRecombDiffTimeProb(int min_k, int max_k, int max_d,
                                  int path_d, int a, int path_a) const;
    double self_recomb_prob_slow(int a, int path_a, int minage,
                                 int d, int path_d) const;
    double self_recomb_prob_slow_sum(int a, int path_a, int min_d,
                                     int max_d, int path_d) const;
    double self_recomb_prob(int a, int path_a, int min_d, int max_d,
                            int path_d) const;
    void calc_self_recomb_probs_smcPrime(const LocalTree *tree,
                                         const States &states);

};


// A compressed representation of the switch transition matrix.
//
// This transition matrix is used in the chromosome threading HMM to go between
// one non-recombining alignment block to the next (i.e. switching blocks).
class TransMatrixSwitch
{
public:
    TransMatrixSwitch(int nstates1, int nstates2, int npaths, bool alloc=true) :
        nstates1(nstates1),
        nstates2(nstates2),
        npaths(npaths),
        own_data(false)
    {
        if (alloc)
            allocate(nstates1, nstates2, npaths);
    }

    ~TransMatrixSwitch()
    {
        if (own_data) {
            delete [] determ;
            delete [] determprob;
            delete [] recoalrow;
            delete [] recombrow;
            delete [] recombsrc;
            delete [] recoalsrc;
        }
    }

    // Allocate matrix with dimensions (nstates1, nstates2).
    void allocate(int _nstates1, int _nstates2, int _npaths)
    {
        nstates1 = _nstates1;
        nstates2 = _nstates2;
        npaths = _npaths;

        // NOTE: nstates1 and nstates2 might be zero
        // we still calculate transitions for a state space of size zero

        own_data = true;
        determ = new int [max(nstates1, 1)];
        determprob = new double [max(nstates1, 1)];
        recoalrow = new double [max(nstates2, 1) * npaths];
        recombrow = new double [max(nstates2, 1) * npaths];
        recombsrc = new int [max(nstates1, 1)];
        recoalsrc = new int[max(nstates1, 1)];
    }

    // Log probability of transition from state i to state j.
    inline double get_log(int i, int j) const
    {
        return log(get(i, j));
    }

    // Probability of transition from state i to state j.
    inline double get(int i, int j) const
    {
        if (recoalsrc[i] >= 0) {
            return recoalrow[recoalsrc[i] * nstates2 + j];
        } else if (recombsrc[i] >= 0) {
            return recombrow[recombsrc[i] * nstates2 + j];
        } else if (determ[i] == j) {
            return determprob[i];
        }
        return 0.0;
    }
    void init_probs() {
        for (int i=0; i < max(nstates2,1) * npaths; i++)
            recombrow[i] = recoalrow[i] = 0.0;
    }

    inline void set(int i, int j, double val) {
        if (recombsrc[i] >= 0) {
            recombrow[recombsrc[i] * nstates2 + j] = val;
        } else if (recoalsrc[i] >= 0) {
            recoalrow[recoalsrc[i] * nstates2 + j] = val;
        } else {
            assert(0);
        }
    }


    int nstates1;   // Number of states in beginning block
    int nstates2;   // Number of states in the ending block
    int npaths;     // Maximum number of population paths (there can
                    //  be one recombsrc/recoalsrc per path)

    //recoalsrc[i] = -1 if state i is not recoal src
    //recoalsrc[i] = j means state i is recoal src with path j
    int *recoalsrc;
    int *recombsrc;
    bool own_data;  // If true, delete matrix data when object is deleted
    int *determ;    // Compressed representation of deterministic transitions
                    // determ[i] = j indicates (i --> j) is a deterministic
                    // transition.  determ[i] = -1, means row i has many
                    // transitions.
    double *determprob;  // Unnormalized probability determprob[i] of
                         // transition (i -> determ[i])
    double *recoalrow;   // Transition probabilities for row recoalsrc
    double *recombrow;   // Transition probabilities for row recombsrc
};



//=============================================================================

/*void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix,
    bool internal=false, int minage=0);*/
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
                          const States &states, const LineageCounts *lineages,
                          double **transprob);
void get_transition_probs(const LocalTree *tree, const ArgModel *model,
                           const States &states, const LineageCounts *lineages,
                           const TransMatrix *matrix, double **transprob);


void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages,
    TransMatrixSwitch *transmat_switch, bool internal=false);
void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, double **transprob);
void get_transition_probs_switch(const TransMatrixSwitch *matrix,
                                 double **transprob);
void calc_transition_probs_switch_internal(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages,
    TransMatrixSwitch *transmat_switch);
void calc_transition_probs_switch_internal2(
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages,
    TransMatrixSwitch *transmat_switch);


double calc_state_priors(const ArgModel *model,
    int time, int pop_path,
    int **nbranches_pop, int **ncoals_pop,
    int minage);
void calc_state_priors(const States &states, const LineageCounts *lineages,
    const ArgModel *model, double *priors,
    const int minage=0);

void get_deterministic_transitions(
    const LocalTree *last_tree, const LocalTree *tree,
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    NodeStateLookup &state2_lookup,
    const ArgModel *model, int *next_states, bool internal=false);


bool assert_transmat(const LocalTree *tree, const ArgModel *model,
                     const TransMatrix *matrix);


void get_recomb_transition_switch(const ArgModel *model,
    const LocalTree *tree, const LocalTree *last_tree,
    const Spr &spr, const int *mapping,
    NodeStateLookup &state2_lookup,
    int next_states[2], int last_path, int minage2);

} // namespace argweaver

#endif // ARGWEAVER_TRANS_H
