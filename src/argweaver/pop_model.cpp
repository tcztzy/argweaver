#ifdef ARGWEAVER_MPI
#include "mpi.h"
#include "mcmcmc.h"
#endif

#include <string.h>
#include "logging.h"
#include "pop_model.h"
#include "local_tree.h"


namespace argweaver {


void PopulationPath::print() const {
    printf("prob=%f, path=%i (t=0)", prob, pop[0]);
    for (unsigned int t=1; t < pop.size(); t++) {
        if (pop[t] != pop[t-1]) {
            printf(", %i (t=%i)", pop[t], t);
        }
    }
}

void PopulationTree::print_all_paths() const {
    for (unsigned int i=0; i < all_paths.size(); i++) {
        printf("path %i:", i);
        all_paths[i].print();
        printf("\n");
    }
}

void PathProb::print() const {
    set<int>::iterator it = path.begin();
    printf("prob=%f paths=(%i", prob, *it);
    for (++it; it != path.end(); ++it)
        printf(", %i", *it);
    printf(")");
}


void PopulationTree::print_sub_paths() const {
    for (int t1=0; t1 < model->ntimes; t1++) {
        for (int t2=t1; t2 < model->ntimes; t2++) {
            for (int p1=0; p1 < npop; p1++) {
                for (int p2=0; p2 < npop; p2++) {
                    printf("sub_path(%i, %i, %i, %i): ", t1, t2, p1, p2);
                    sub_paths[t1][t2][p1][p2].print();
                    printf("\n");
                }
            }
        }
    }
}

void PopulationTree::print() const {
    printf("All complete paths:\n");
    print_all_paths();
    printf("\nAll sub-paths:\n");
    print_sub_paths();
}


PopulationTree::PopulationTree(int npop, const ArgModel *model) :
    npop(npop), model(model) {
    int ntime2 = 2 * model->ntimes - 1;  // number of half-time intervals
    mig_matrix.clear();
    mig_matrix.resize(ntime2);
    for (int i=0; i < ntime2; i++)
        mig_matrix[i].resize(npop);
    sub_paths = NULL;
    max_matching_path = NULL;
}


PopulationTree::PopulationTree(const PopulationTree &other) {
    npop = other.npop;
    model = other.model;
    //    mig_matrix.init(npop);
    //    mig_matrix.copy(other.mig_matrix);
    mig_matrix = other.mig_matrix;
    sub_paths = NULL;
    max_matching_path = NULL;
    if (npop > 0) set_up_population_paths();
}


PopulationTree::~PopulationTree() {
    if (sub_paths != NULL) {
        for (int i=0; i < model->ntimes; i++) {
            for (int j=i; j < model->ntimes; j++) {
                for (int k=0; k < npop; k++) {
                    delete [] sub_paths[i][j][k];
                }
                delete [] sub_paths[i][j];
            }
            delete [] sub_paths[i];
        }
        delete sub_paths;
    }
    if (max_matching_path != NULL) {
        for (unsigned int i=0; i < all_paths.size(); i++) {
            for (unsigned int j=0; j < all_paths.size(); j++) {
                delete [] max_matching_path[i][j];
            }
            delete max_matching_path[i];
        }
        delete max_matching_path;
    }
}

void PopulationTree::update_npop(int new_npop) {
    if (npop != new_npop) {
        for (int i=0; i < 2*model->ntimes - 1; i++)
            mig_matrix[i].resize(new_npop);
        npop = new_npop;
     }
}


void PopulationTree::add_migration(int t, int from_pop, int to_pop, double prob) {
     if (t%2 != 1)
         exitError("Error: addMigration expects only odd times\n");
     if (t < 0 || t > 2 * model->ntimes - 3)
         exitError("Error: migratation index %i out of range\n", t);
     if (from_pop >= npop || from_pop < 0)
         exitError("Error: from_pop (%i) out of range\n", from_pop);
     if (to_pop >= npop || to_pop < 0)
         exitError("Error: to_pop (%i) out of range\n", to_pop);
     mig_matrix[t].update(from_pop, to_pop, prob);
 }




bool PopulationTree::paths_equal(int path1, int path2, int t1, int t2) const {
    if (path1 == path2) return true;
    for (int t=t1; t <= t2; t++)
        if (all_paths[path1].get(t) != all_paths[path2].get(t))
            return false;
    return true;
}


void PathProb::update_prob(const vector<PopulationPath> &all_paths) {
    prob = 0;
    for (set<int>::iterator it = path.begin(); it != path.end(); ++it)
        prob += all_paths[*it].prob;
}

void PopulationTree::update_population_probs() {
    int ntime = model->ntimes;
    for (unsigned int i=0; i < all_paths.size(); i++) {
        all_paths[i].prob = 1.0;
        for (int t=1; t < ntime; t++) {
            all_paths[i].prob *= mig_matrix[2*t-1].get(all_paths[i].pop[t-1],
                                                       all_paths[i].pop[t]);
        }
    }

    // update sub_paths probs
    for (int t1=0; t1 < ntime; t1++) {
        for (int t2=t1; t2 < ntime; t2++) {
            for (int p1=0; p1 < npop; p1++) {
                for (int p2=0; p1 < npop; p1++) {
                    sub_paths[t1][t2][p1][p2].update_probs(all_paths);
                }
            }
        }
    }
}


int PopulationTree::find_sub_path(int path,
                                  const SubPath &subpath,
                                  int t1, int t2) {
    for (unsigned int i=0; i < subpath.size(); i++) {
        int cur_path = subpath.first_path(i);
        if (paths_equal(path, cur_path, t1, t2))
            return i;
    }
    return -1;
}

void PopulationTree::getAllPopulationPathsRec(PopulationPath &curpath,
                                              int cur_time, int end_time,
                                              int cur_pop) {
    if (cur_time == end_time) {
        all_paths.push_back(curpath);
        return;
    }
    int t = cur_time * 2 + 1;
    bool called=false;
    for (int next_pop=0; next_pop < npop; next_pop++) {
        double prob = mig_matrix[t].get(cur_pop, next_pop);
        if (prob > 0.0) {
            if (called) {
                PopulationPath newpath(curpath);
                newpath.set(cur_time + 1, next_pop, prob);
                getAllPopulationPathsRec(newpath, cur_time + 1, end_time,
                                         next_pop);
            } else {
                called = true;
                curpath.set(cur_time + 1, next_pop, prob);
                getAllPopulationPathsRec(curpath, cur_time + 1, end_time,
                                         next_pop);
            }
        }
    }
}



void PopulationTree::set_up_population_paths() {
    int ntime = model->ntimes;

    all_paths.clear();
    for (int p1 = 0; p1 < npop; p1++) {
        PopulationPath p(ntime);
        p.set(0, p1);
        getAllPopulationPathsRec(p, 0, ntime-1, p1);
    }
    printLog(LOG_LOW, "Total number of possible paths: %i\n", all_paths.size());

    // check that all paths end in same population
    int lastpop=all_paths[0].pop[ntime - 1];
    for (unsigned int i=1; i < all_paths.size(); i++) {
        if (all_paths[i].get(ntime - 1) != lastpop)
            exitError("Error: populations do not converge by final time\n");
    }

    // now get all possible paths for each possible start/end time start/end pop
    sub_paths = new SubPath ***[ntime];
    for (int t1=0; t1 < ntime; t1++) {
        sub_paths[t1] = new SubPath **[ntime];
        for (int t2=t1; t2 < ntime; t2++) {
            sub_paths[t1][t2] = new SubPath *[npop];
            for (int p1=0; p1 < npop; p1++) {
                sub_paths[t1][t2][p1] = new SubPath[npop];
                for (int p2=0; p2 < npop; p2++) {
                    sub_paths[t1][t2][p1][p2].set_size(all_paths.size());

                    // find all paths which have start t1,p1 and end t2,p2
                    for (unsigned int i=0; i < all_paths.size(); i++) {
                        if (all_paths[i].get(t1) == p1 &&
                            all_paths[i].get(t2) == p2) {
                            double prob = all_paths[i].prob;

                            // path i has this start/end.
                           // search previous paths for identical sub-path
                            int other = find_sub_path(i,
                                                      sub_paths[t1][t2][p1][p2],
                                                      t1, t2);
                            if (other >= 0) {
                                sub_paths[t1][t2][p1][p2].add_path_to_subpath(
                                     i, other, prob);
                            } else {
                                sub_paths[t1][t2][p1][p2].new_subpath(i, prob);
                            }
                        }
                    }
                }
            }
        }
    }

    max_matching_path = new int **[all_paths.size()];
    for (unsigned int i=0; i < all_paths.size(); i++) {
        max_matching_path[i] = new int*[all_paths.size()];
        for (unsigned int j=0; j < all_paths.size(); j++) {
            max_matching_path[i][j] = new int[ntime];
            for (int t=0; t < ntime; t++) {
                if (all_paths[i].get(t) != all_paths[j].get(t))
                    max_matching_path[i][j][t] = -1;
                else {
                    max_matching_path[i][j][t] = t;
                    for (int t1=t+1; t1 < model->ntimes; t1++) {
                        if (all_paths[i].get(t1) == all_paths[j].get(t1))
                            max_matching_path[i][j][t] = t1;
                        else break;
                    }
                }
            }
        }
    }

}

 int PopulationTree::final_pop() const {
     return all_paths[0].get(model->ntimes - 1);
 }


// This could be made more efficient if necessary
int PopulationTree::consistent_path(int path1, int path2,
                                    int t1, int t2, int t3,
                                    bool require_exists) const {
    if (path1 == path2) return path1;
    if (path1 == -1 || path2 == -1) return -1;
    if (t1 < 0 || t1 >= model->ntimes) t1 = model->ntimes - 1;
    if (t2 < 0 || t2 >= model->ntimes) t2 = model->ntimes - 1;
    if (t3 < 0 || t3 >= model->ntimes) t3 = model->ntimes - 1;
    assert(t1 <= t2);
    assert(t2 <= t3);
    int p1 = all_paths[path1].get(t1);
    int p2 = all_paths[path1].get(t2);
    if (p2 != all_paths[path2].get(t2)) {
        if (require_exists)
            exitError("No consistent path found\n");
        return -1;
    }
    int p3 = all_paths[path2].get(t3);
    SubPath possible_paths;
    if (sub_paths[t1][t2][p1][p2].size() <=
        sub_paths[t2][t3][p2][p3].size())
        possible_paths = sub_paths[t1][t2][p1][p2];
    else possible_paths = sub_paths[t2][t3][p2][p3];
    for (unsigned int i=0; i < possible_paths.size(); i++) {
        int path = possible_paths.first_path(i);
        if (paths_equal(path, path1, t1, t2) &&
            paths_equal(path, path2, t2, t3))
            return path;
    }
    if (require_exists)
        exitError("No consistent path found\n");
    return -1;
}


int PopulationTree::path_to_root(const LocalNode *nodes, int node) const {
    assert(node != -1);
    int path = nodes[node].pop_path;
    int parent = nodes[node].parent;
    while (true) {
        path = consistent_path(path,
                               nodes[parent].pop_path,
                               nodes[node].age,
                               nodes[parent].age,
                               ( nodes[parent].parent == -1 ?
                                 model->ntimes - 1 :
                                 nodes[nodes[parent].parent].age));
        node = parent;
        parent = nodes[node].parent;
        if (parent == -1) break;
    }
    return path;
}

// not implemented efficiently, as only called a few times
int get_closest_half_time(double tgen, const double *time_steps, int ntime) {
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


/**
Read in a "population tree file". Format should be as in this example:

# Comments start with hash-tag
npop 2   # 2 populations
mig 50 0 1 0.01  # migration at time 50 from pop 0 to pop 1 w/prob 0.01
div 100 1 0  # divergence at time 100
### end example

In this example, there are two populations, so they are numbered 0 and 1. At 50
generations ago, there is a migration from population 0 to population 1 with
probability 0.01. This is looking backwards in time, so this actually represents
a forward-time migration from population 1 to 0.

The next line indicates that the populations diverged at generation 100. This is
implemented as a migration from 1 to 0 at time 100 with 100% probability. So
after time 100 there is only a single population, pop 0.

**/
void read_population_tree(FILE *infile, PopulationTree *pop_tree) {
    char str[10000], c;
    int ntimes = pop_tree->model->ntimes;
    double time_steps[2*ntimes-1];
    int npop=-1;
    time_steps[0] = 0.0;
    for (int i=1; i < 2*ntimes-1; i++)
        time_steps[i] = time_steps[i-1] + pop_tree->model->coal_time_steps[i-1];
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
            pop_tree->update_npop(npop);
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
            int tidx = get_closest_half_time(tgen, time_steps, 2*ntimes-1);
            if (fabs(tgen - time_steps[tidx]) > 1)
                printLog(LOG_LOW, "Using time %f instead of %f for %s event in "
                         "popfile\n",
                         time_steps[tidx], tgen, str);
            pop_tree->add_migration(tidx, from_pop, to_pop, prob);
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
    pop_tree->set_up_population_paths();
}

}  //namespace argweaver
