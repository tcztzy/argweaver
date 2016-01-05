//=============================================================================
// ArgHmm multiple population model


#ifndef ARGWEAVER_POP_MODEL_H
#define ARGWEAVER_POP_MODEL_H

// c/c++ includes
#include <math.h>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <vector>
#include <string>

// arghmm includes
#include "common.h"
#include "model.h"

namespace argweaver {

class ArgModel;

class PopulationPath {
 public:
    PopulationPath(int ntime) {
        path = vector<int>(ntime, -1);
        pathprob = 1.0;
    }
    PopulationPath(const PopulationPath &other) {
        pathprob = other.pathprob;
        path = vector<int>(other.path);
    }
    bool operator< (const PopulationPath &other) const {
        if (path.size() != other.path.size())
            return ( path.size() < other.path.size());
        if (pathprob != other.pathprob)
            return (pathprob < other.pathprob);
        for (unsigned int i=0; i < path.size(); i++)
            if (path[i] != other.path[i]) return (path[i] < other.path[i]);
        return false;
    }
    void set(int time, int pop, double prob=1.0) {
        path[time] = pop;
        pathprob *= prob;
    }
    void print() {
        int start_pop;
        for (unsigned int start_pop=0; start_pop < path.size(); start_pop++)
            if (path[start_pop] != -1) break;
        printf("path: %i (t=%i)", path[start_pop], start_pop);
        int prev = path[start_pop];
        for (unsigned int i=start_pop + 1; i < path.size(); i++) {
            if (path[i] == -1) break;
            if (path[i] != prev) {
                printf(", %i (t=%i)", path[i], i);
                prev = path[i];
            }
        }
        printf("\tprob=%f\n", pathprob);
    }

    double pathprob;
    vector<int> path;
};  /* class PopulationPath */


/* Represents the migration matrix for a single time period.
   Is simply a two-dimensional matrix giving the probability of
   migrating from each pop to each other pop.
 */
class MigMatrix {
 public:
   MigMatrix() {
       npop=0;
   }
   MigMatrix(int npop) : 
    npop(npop) {
        init();
   }
    MigMatrix(const MigMatrix &other) {
        npop = other.npop;
        init();
        if (npop > 0)
            std::copy(other.mat, other.mat + npop*npop, mat);
    }
    ~MigMatrix() {
        if (npop > 0) delete [] mat;
    }
    /*    void copy(MigMatrix &other) {
        std::copy(other.mat, other.mat + pop*pop, mat);
        }*/
    void set(int from_pop, int to_pop, double val) {
        mat[from_pop * npop + to_pop] = val;
    }
    void update(int from_pop, int to_pop, double val) {
        set(from_pop, to_pop, val);
        set(from_pop, from_pop, 1.0-val);
    }
    void resize(int new_npop) {
        if (new_npop != npop) {
            if (npop > 0) delete [] mat;
            npop = new_npop;
            init();
        }
    }
    double get(int from_pop, int to_pop) const {
        return mat[from_pop * npop + to_pop];
    }

    int npop;
    double *mat;

 protected:
    void init() {
        if (npop > 0) {
            mat = new double[npop*npop]();
            for (int i=0; i < npop; i++)
                set(i, i, 1.0);
        }
    }    
};



/* PopulationTree
   is implemented as a vector of migration matrices. One element for each
   half time interval. 
 */
class PopulationTree {

 public:
  PopulationTree(int npop, const ArgModel *model);
  PopulationTree(const PopulationTree &other);
  ~PopulationTree();
  void update_npop(int new_npop);
  void add_migration(int t, int from_pop, int to_pop, double prob);
  int set_up_population_paths();

  //  npop is the total number of non-ancestral populations
  int npop;
  const ArgModel *model;

  vector<MigMatrix> mig_matrix;


 // want a set of population paths for each starting and ending time (whole times only)
 // to/from every pair of populations
 // and a probability associated with each
 set<PopulationPath> ****paths;

 private:
    void getAllPopulationPathsRec(PopulationPath &curpath,
                                  int start_time, int cur_time, int end_time, 
                                  int start_pop, int cur_pop);

};  /* class PopulationTree */    

bool read_population_tree(FILE *infile, PopulationTree *pop_tree);

int get_closest_half_time(double treal, const double *time_steps, int ntime);
    

} // namespace argweaver


#endif // ARGWEAVER_POP_MODEL_H
