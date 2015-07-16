#ifdef ARGWEAVER_MPI
#include "mpi.h"
#include "mcmcmc.h"
#endif

#include <string.h>
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
        pop_paths.push_back(PopulationPath(ntime, i, 1.0));
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
    for (unsigned int i=0; i < pop_paths.size(); i++)
        pop_paths[i].print();

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


    // not implemented efficiently, as only called a few times
int PopulationTree::get_closest_half_time(double tgen, const double *time_steps, int ntime) {
    int closest=1;
    double dist = fabs(tgen - time_steps[closest]);
    for (int i=3; i < ntime; i+=2) {
        double tempdist = fabs(tgen - time_steps[i]);
        if (tempdist < dist) {
            dist = tempdist;
            closest = i;
        }
    }
    return closest;
}

PopulationTree::PopulationTree(string filename, const double *coal_time_steps, int ntime) : ntime(ntime){
    FILE *infile = fopen(filename.c_str(), "r");
    char str[10000], c;
    double time_steps[2*ntime-1];
    if (infile == NULL)
        exitError("Error opening population file %s\n", filename.c_str());
    npop = -1;
    time_steps[0] = 0;
    for (int i=1; i < 2*ntime-1; i++)
        time_steps[i] = coal_time_steps[i-1] + time_steps[i-1];
    while (1==fscanf(infile, "%s", str)) {
        if (str[0] == '#') {
            while ('\n' != (c=fgetc(infile)) && c != EOF);
            if (c == EOF) break;
            continue;
        }
        if (npop == -1) {
            if (strcmp(str, "npop") != 0)
                exitError("Error: Expected population file to start with 'npop'\n");
            if (1 != fscanf(infile, "%i", &npop))
                exitError("Premature end of population file\n");
            printLog(LOG_LOW, "Number of populations: %i\n", npop);
            mig_events = new vector<MigrationEvent>*[2*ntime-1];
            for (int i=0; i < ntime-1; i++) {
                mig_events[2*i+1] = new vector<MigrationEvent>[npop];
                for (int j=0; j < npop; j++)
                    mig_events[2*i+1][j].clear();
            }
        } else if (strcmp(str, "div")==0 || strcmp(str, "mig")==0) {
            double tgen;
            int from_pop, to_pop;
            double prob = 1.0;
            bool is_mig=(strcmp(str, "mig")==0);
            if (3 != fscanf(infile, "%lf %i %i", &tgen, &from_pop, &to_pop))
                exitError("Premature end of population file\n");
            if (is_mig) {
                if (1 != fscanf(infile, "%lf", &prob))
                    exitError("Premature end of population file\n");
            }
            // todo: need to get closest time
            int tidx = get_closest_half_time(tgen, time_steps, ntime);
            if (fabs(tgen - time_steps[tidx]) > 1)
                printLog(LOG_LOW, "Using time %f instead of %f for %s event in popfile\n",
                         time_steps[tidx], tgen, str);
            // for now, optimize all migrations but not divergences
            // except optimization not yet implemented anyway
            addMigration(tidx, from_pop, to_pop, prob, is_mig);
        }
        bool incomment=false;
        while ('\n' != (c=fgetc(infile)) && c != EOF) {
            if (c== '#') incomment=true;
            if (!isspace(c) && !incomment) {
                printError("Format error in population file\n");
            }
        }
        if (c==EOF) break;
    }
    fclose(infile);

    setUpPopulationPaths();
}

}  //namespace argweaver
