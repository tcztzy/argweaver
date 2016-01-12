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
        pop = vector<int>(ntime, -1);
        prob = 1.0;
    }
    PopulationPath(const PopulationPath &other) {
        prob = other.prob;
        pop = vector<int>(other.pop);
    }
    bool operator< (const PopulationPath &other) const {
        if (pop.size() != other.pop.size())
            return ( pop.size() < other.pop.size());
        if (prob != other.prob)
            return (prob < other.prob);
        for (unsigned int i=0; i < pop.size(); i++)
            if (pop[i] != other.pop[i]) return (pop[i] < other.pop[i]);
        return false;
    }
    void set(int time, int t_pop, double currprob=1.0) {
        pop[time] = t_pop;
        prob *= currprob;
    }
    int get(int time) const {
        return pop[time];
    }
    void print() const;

    double prob;
    vector<int> pop;
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


class PathProb {
 public:
   PathProb(int path_idx, double prob) : prob(prob) {
        path.clear();
        path.insert(path_idx);
    };
   void update_prob(const vector <PopulationPath> &all_paths);
   void print() const;
   set<int> path;
   double prob;
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
  void set_up_population_paths();
  void update_population_probs();
  bool paths_equal(int path1, int path2, int t1, int t2) const;
  void print_all_paths() const;
  void print_sub_path(vector<PathProb> &subpath) const;
  void print_sub_paths() const;
  void print() const;

  //  npop is the total number of non-ancestral populations
  int npop;
  const ArgModel *model;

  vector<MigMatrix> mig_matrix;


  // this is the set of distinct paths from t=0 to t=ntimes-1
  vector<PopulationPath> all_paths;

  //sub_paths[t1][t2][p1][p2] is a vector of all paths that lead from
  // time t1, pop p1 to time t2, pop p2. If multiple paths are identical
  // during that interval, then only the first is included in the vector,
  // and the probabilities are the identical paths are summed
  vector<PathProb> ****sub_paths;



 int get_pop(int path, int time) const {
     return all_paths[path].get(time);
 }

 int final_pop() const {
     return all_paths[0].get(model->ntimes - 1);
 }

 // Returns a path consistent with path1 from t1_start to t1_end,
 // and path2 from t2_start to t2_end.
 // Assumes t1_start <= t1_end and t2_start <= t2_end
 // Throws an error if no such path exists
 // Note that there is not always a unique answer;
 // the first one found is returned
 int consistent_path(int path1, int t1_start, int t1_end,
                     int path2, int t2_start, int t2_end) const;

 private:
    void getAllPopulationPathsRec(PopulationPath &curpath,
                                  int cur_time, int end_time, int cur_pop);
    int find_sub_path(int path, const vector<PathProb> &path_vec,
                      int t1, int t2);

};  /* class PopulationTree */

void read_population_tree(FILE *infile, PopulationTree *pop_tree);

int get_closest_half_time(double treal, const double *time_steps, int ntime);


} // namespace argweaver


#endif // ARGWEAVER_POP_MODEL_H
