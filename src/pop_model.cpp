#ifdef ARGWEAVER_MPI
#include "mpi.h"
#include "mcmcmc.h"
#endif

#include "logging.h"
#include "pop_model.h"



namespace argweaver {

void PopulationTree::getAllPopulationPathsRec(int idx, int start_pop, int time,
                                              double tol) 
{
    if (time == ntime-1) return;
    bool called=false;
    double prob = 0.0;
    double no_opt_prob=0.0;
    int mig_time = 2*time+1;
    vector<MigrationEvent> *mig_list = &mig_events[mig_time][start_pop];
    for (unsigned int i=0; i < mig_list->size(); i++) {
        prob += (*mig_list)[i].prob;
        if (! (*mig_list)[i].optimize)
            no_opt_prob += (*mig_list)[i].prob;
    }

    // migration events don't sum to 1; so can stay in current pop
    if (no_opt_prob < 1-tol) {
        pop_paths[idx].pathprob[time+1] = 
            pop_paths[idx].pathprob[time] * (1.0 - prob);
        pop_paths[idx].pop[time+1] = start_pop;
        getAllPopulationPathsRec(idx, start_pop, time+1, tol);
        called=true;
    }

    for (unsigned int i = 0; i < mig_list->size(); i++) {
        int new_idx;
        if (called) {
            new_idx = pop_paths.size();
            pop_paths.push_back(PopulationPath(pop_paths[idx]));
        } else {
            new_idx = idx;
            called = true;
        }
        pop_paths[new_idx].pop[time+1] = 
            (*mig_list)[i].to_pop;
        pop_paths[new_idx].pathprob[time+1] =
            pop_paths[new_idx].pathprob[time] * (*mig_list)[i].prob;
        getAllPopulationPathsRec(new_idx, (*mig_list)[i].to_pop, 
                                 time+1, tol);
    }
}


int PopulationTree::setUpPopulationPaths() {
    pop_paths.clear();
    for (int i=0; i < npop; i++) {
        pop_paths.push_back(PopulationPath(ntime, i, 0));
        getAllPopulationPathsRec(pop_paths.size()-1, 
                                 i, 0);
    }

    // ensure that all paths end in same population
    for (unsigned int i=1; i < pop_paths.size(); i++) {
        if (pop_paths[i].pop[ntime-1] != pop_paths[0].pop[ntime - 1]) {
            fprintf(stderr, "Error: populations do not converge\n");
            assert(false);
        }
    }

    printLog(LOG_LOW, "Number of population paths: %i\n", pop_paths.size());
    
    // now want 4-d array which contains pointers to all paths which
    // go from time i, pop j, to time k, pop l
    pathPtr = new vector<PopulationPath*>***[ntime];
    for (int i=0; i < ntime; i++) {
        pathPtr[i] = new vector<PopulationPath*>**[npop];
        for (int j=0; j < npop; j++) {
            pathPtr[i][j] = new vector<PopulationPath*>*[ntime];
            for (int k=0; k < ntime; k++) {
                pathPtr[i][j][k] = new vector<PopulationPath*>[npop];
                for (int l=0; l < npop; l++) {
                    pathPtr[i][j][k][l] = vector<PopulationPath*>();
                }
            }
        }
    }
    for (unsigned int idx=0; idx < pop_paths.size(); idx++) {
        PopulationPath *path = &pop_paths[idx];
        for (int i=0; i < ntime; i++) {
            for (int j=i; j < ntime; j++) {
                pathPtr[i][path->pop[i]][j][path->pop[j]].push_back(&pop_paths[idx]);
            }
        }
    }
    return pop_paths.size();
}

}  //namespace argweaver
