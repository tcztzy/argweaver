//=============================================================================
// transitions

#ifndef ARGWEAVER_TRANS_H
#define ARGWEAVER_TRANS_H

#include "common.h"
#include "local_tree.h"
#include "model.h"
#include "states.h"

namespace argweaver {


// A compressed representation of the transition matrix.
//
// This transition matrix is used in the chromosome threading HMM within
// one non-recombining alignment block.
class TransMatrix
{
public:
    TransMatrix(const ArgModel *model, int nstates) :
        nstates(nstates),
        internal(false)
    {
        initialize(model);
    }

    ~TransMatrix()
    {
        if (npaths > 1) {
            for (int i=0; i < npaths; i++) {
                for (int j=0; j < npaths; j++) {
                    for (int k=0; k < ntimes; k++) {
                        delete [] paths_equal[i][j][k];
                    }
                    delete [] paths_equal[i][j];
                }
                delete [] paths_equal[i];
            }
            delete [] paths_equal;

            for (int i=0; i < npaths; i++) {
                for (int j=0; j < npaths; j++) {
                    delete [] paths_equal_max[i][j];
                }
                delete [] paths_equal_max[i];
            }
            delete [] paths_equal_max;
        }
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
        delete [] path_prob;
        delete [] data_alloc;
    }

    // allocate space for transition matrix
    // and initialize paths_equal matrix
    void initialize(const ArgModel *model)
    {
        ntimes = model->ntimes;
        npaths = model->num_pop_paths();
        int data_len = ntimes*2 + npaths*ntimes*4 + npaths*npaths*ntimes*3;
        data_alloc = new double [data_len];
        E = new double* [npaths];
        lnB = new double** [npaths];
        lnE2 = new double** [npaths];
        lnNegG1 = new double** [npaths];
        G2 = new double* [npaths];
        G3 = new double* [npaths];
        path_prob = new double* [npaths];
        D = &data_alloc[0];
        norecombs = &data_alloc[ntimes];
        int idx=ntimes*2;
        for (int i=0; i < npaths; i++) {
            E[i] = &data_alloc[idx]; idx += ntimes;
            G2[i] = &data_alloc[idx]; idx += ntimes;
            G3[i] = &data_alloc[idx]; idx += ntimes;
            path_prob[i] = &data_alloc[idx]; idx += ntimes;
            lnB[i] = new double* [npaths];
            lnE2[i] = new double* [npaths];
            lnNegG1[i] = new double* [npaths];
            for (int j=0; j < npaths; j++) {
                lnB[i][j] = &data_alloc[idx]; idx += ntimes;
                lnE2[i][j] = &data_alloc[idx]; idx += ntimes;
                lnNegG1[i][j] = &data_alloc[idx]; idx += ntimes;
            }
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
                        minage, node1 == node2);
    }

    // Returns the probability of transition from state1 with time 'a'
    // to state2 with time 'b'.  The probability also depends on whether
    // the node changes between states ('same_node') or whether there is
    // a minimum age ('minage') allowed for the state.
    inline double get_time(int a, int b, int c,
                           int path_a, int path_b, int path_c,
                           int minage, bool same_node) const
    {
        if (a < minage || b < minage)
            return 0.0;

        const int p = ( npaths == 1 ? ntimes :
                        paths_equal_max[path_a][path_b][minage] );
        if (p == -1) return 0.0;

        double term1 = D[a] * E[path_b][b] * path_prob[path_b][b];
        double minage_term = 0.0;
        double prob;
        if (minage > 0) {
            minage_term = exp(lnE2[path_b][path_b][b] +
                              lnB[path_b][path_b][minage-1]);
        }
        if (p < a && p < b) {
            prob = term1 * (exp(lnE2[path_b][path_b][b] + lnB[path_b][path_b][p])
                            - minage_term);
        } else if (a <= p && a < b) {
            prob = term1 * (exp(lnE2[path_b][path_b][b] + lnB[path_b][path_b][a]) -
                            exp(lnE2[path_b][path_b][b] + lnNegG1[path_b][path_b][a])
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
        if (a == b && (npaths == 1 || paths_equal[path_a][path_b][minage][a]))
            prob += norecombs[a];

        if (npaths > 1 && a < b && !(paths_equal[path_b][path_c][a][b] &&
                                     paths_equal[path_a][path_b][minage][a])) {
            if (isnan(prob))
                assert(false);
            return prob;
        }
        if (npaths > 1 && b <= a && !(paths_equal[path_a][path_c][b][a] &&
                                      paths_equal[path_b][path_a][minage][b]))
            return prob;

        term1 = D[a] * E[path_c][b] * path_prob[path_c][b];
        minage_term = 0.0;
        if (c > 0)
            minage_term = exp(lnE2[path_c][path_b][b] + lnB[path_c][path_b][c-1]);

        if (a < b) {
            prob += term1 * (exp(lnE2[path_c][path_b][b] + lnB[path_c][path_b][a]) -
                            exp(lnE2[path_c][path_b][b] + lnNegG1[path_c][path_b][a])
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
    }

    // Log probability of transition from state i to state j.
    inline double get_log(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        return log(get(tree, states, i, j));
    }


    int ntimes;     // Number of time steps in model.
    int nstates;    // Number of states in HMM.
    int npaths;

    double *data_alloc;

    double *D;       // Intermediate terms in calculating entries in the full
    double **E;      // transition matrix.
    double ***lnB;
    double ***lnE2;
    double ***lnNegG1;
    double **G2;
    double **G3;
    double **path_prob;
    double *norecombs;
    bool ****paths_equal;
    int ***paths_equal_max;

    bool internal;  // If true, this matrix is for threading an internal branch.
    int minage;     // Minimum age of a state we can consider (due to threading
                    // an internal branch).
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

void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix,
    bool internal=false, int minage=0);
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
