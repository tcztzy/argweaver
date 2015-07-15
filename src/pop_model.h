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

// arghmm includes
#include "common.h"

namespace argweaver {

class MigrationMatrix
{
 public:
 MigrationMatrix() :
    npop(0), ntime(0), matrix(NULL) {}

 MigrationMatrix(int npop, int ntime) :
    npop(npop), ntime(ntime)
    {
	matrix = new double**[ntime];
	for (int i=0; i < ntime; i++) {
	    matrix[i] = new double*[npop];
	    for (int j=0; j < ntime; j++) {
		matrix[i][j] = new double[npop]();
	    }
	}
    }

 ~MigrationMatrix() {
     if (matrix != NULL) {
	 for (int i=0; i < ntime; i++) {
	     for (int j=0; j < npop; j++)
		 delete matrix[i][j];
	     delete matrix[i];
	 }
	 delete matrix;
     }
 }
     
 int npop;
 int ntime;
 double ***matrix;
};


class PopulationPath {
 public:
   PopulationPath(int ntime, int start_pop, double start_prob=1.0) : 
      ntime(ntime) {
        pop = new int[ntime]();
        pathprob = new double[ntime]();
        pop[0] = start_pop;
        pathprob[0]=start_prob;
    }
    PopulationPath(const PopulationPath &other) {
        pop = new int [ntime];
        pathprob = new double[ntime];
        for (int i=0; i < ntime; i++) {
            pop[i] = other.pop[i];
            pathprob[i] = other.pathprob[i];
        }
    }
    ~PopulationPath() {
        delete pop;
        delete pathprob;
    }
    int ntime;
    int *pop;
    double *pathprob;
};


class MigrationEvent {
 public:
   MigrationEvent(int from_pop, int to_pop, double prob, bool optimize=false) :
    from_pop(from_pop), to_pop(to_pop), prob(prob), optimize(optimize) {
    }

    int from_pop; 
    int to_pop;
    double prob;
    bool optimize;
};


class PopulationTree {
 public:
 PopulationTree(int npop, int ntime) :
    npop(npop), ntime(ntime) {
        mig_events = new vector<MigrationEvent>*[2*ntime-1];
        for (int i=0; i < ntime-1; i++) {
            mig_events[2*i+1] = new vector<MigrationEvent>[npop];
            for (int j=0; j < npop; j++)
                mig_events[2*i+1][j].clear();
        }
        pop_paths.clear();
        pathPtr=NULL;
  }

 ~PopulationTree() {
     for (int i=0; i < ntime-1; i++)
         delete mig_events[2*i+1];
     delete mig_events;
     if (pathPtr != NULL) {
         for (int i=0; i < ntime; i++) {
             for (int j=0; j < npop; j++) {
                 for (int k=0; k < ntime; k++)
                     delete pathPtr[i][j][k];
                 delete pathPtr[i][j];
             }
             delete pathPtr[i];
         }
         delete pathPtr;
     }
  }

 void addMigration(int t, int from_pop, int to_pop, double prob, bool optimize=false) {
     if (t%2 != 1) {
         fprintf(stderr, "Error: addMigration expects only odd times\n");
         assert(false);
     }
     if (t < 0 || t > 2*ntime-3) {
         fprintf(stderr, "Error: migratation index %i out of range\n", t);
         assert(false);
     }
     MigrationEvent mig = MigrationEvent(from_pop, to_pop, prob, optimize);
     mig_events[t][from_pop].push_back(mig);
 }
 int setUpPopulationPaths();

 int npop;
 int ntime;
 // mig_events[i][j] is list of all migration events at time i from pop j (i in half-time intervals)
 vector<MigrationEvent> **mig_events;
 // pop_paths is list of all possible paths through populations over time
 vector<PopulationPath> pop_paths;
 // pathPtr[i][j][k][l] is list of pointers to all pop_paths that start at time i, pop j and go
 // to time k, pop l
 vector<PopulationPath*> ****pathPtr;

 private:
 void getAllPopulationPathsRec(int idx, int start_pop, int time=0, double tol=1.0e-6);
};


/*
class PopulationTree {
 public:
    // 
 PopulationTree(int npop_leaf, int ntime, double *coal_time_steps) :
    npop_leaf(npop_leaf), ntime(ntime) {
    npop = 2*npop_present - 1;
    mig_matrix = new MigrationMatrix(npop, ntime);
    half_times = new double[ntime-1];
    int j=0;
    for (int i=0; i < ntime - 1; i++) {
	half_times[j++] = coal_time_steps[2*i+1];
    }
    active_times = new set<int>[ntime-1];
    for (int i=0; i < ntime - 1; i++) {
	for (int j=0; j < npop_leaf; j++)
	    active_times[i].insert(j);
    }
 }

 ~PopulationTree() {
     delete mig_matrix;
     delete half_times;
     delete active_times;
     delete active_pops;
 }

 void set_pop_divergence(int pop1, int pop2, int new_pop, double tdiv_exact) {
     int tdiv = get_closest_time(tdiv_exact, half_times, ntime-1);
     if (fabs(half_times[tdiv] - tdiv_exact) > 1) {
	 printLog(LOG_LOW, "Warning: setting tdiv of pops %i and %i to %.0f (closest interval to desired time %.0f)\n",
		  pop1, pop2, half_times[tdiv], tdiv_exact);
     }
     for (int i=0; i < npop; i++) {
	 mig_matrix.matrix[tdiv][pop1][i] = (i == new_pop ? 1.0 : 0.0);
	 mig_matrix.matrix[tdiv][pop2][i] = (i == new_pop ? 1.0 : 0.0);
	 if (i != pop1 && i != pop2 && i != new_pop) {
	     mig_matrix.matrix[tdiv][i][new_pop] += mig_matrix.matrix[tdiv][i][pop1];
	     mig_matrix.matrix[tdiv][i][new_pop] += mig_matrix.matrix[tdiv][i][pop2];
	     mig_matrix.matrix[tdiv][i][pop1] = 0.0;
	     mig_matrix.matrix[tdiv][i][pop2] = 0.0;
	 }
     }
     for (int i=tdiv+1; i < ntime - 1; i++) {
	 active_times[i].erase(pop1);
	 active_times[i].erase(pop2);
	 active_times[i].insert(new_pop);
     }
 }


 int npop;
 int ntime;
 set<int> *active_times;
 set<int> *active_pops;
 double *half_times;
 MigrationMatrix mig_matrix;

 };*/

} // namespace argweaver


#endif // ARGWEAVER_POP_MODEL_H
