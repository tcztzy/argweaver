#ifdef ARGWEAVER_MPI
#include "mpi.h"
#include "mcmcmc.h"
#endif

#include <string.h>
#include "logging.h"
#include "pop_model.h"
#include "local_tree.h"
#include "Tree.h"

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

void UniquePath::print() const {
    set<int>::iterator it = path.begin();
    printf("start=%i end=%i prob=%f paths=(%i",
           start_time, end_time, prob, *it);
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
    max_matching_path_arr = NULL;
    min_matching_path_arr = NULL;
    max_migrations = -1;
}


PopulationTree::PopulationTree(const PopulationTree &other) {
    npop = other.npop;
    model = other.model;
    //    mig_matrix.init(npop);
    //    mig_matrix.copy(other.mig_matrix);
    mig_matrix = other.mig_matrix;
    sub_paths = NULL;
    max_matching_path_arr = NULL;
    min_matching_path_arr = NULL;
    if (npop > 0) set_up_population_paths();
    update_population_probs();
    max_migrations = -1;
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
        delete [] sub_paths;
    }
    if (max_matching_path_arr != NULL) {
        for (unsigned int i=0; i < all_paths.size(); i++) {
            for (unsigned int j=0; j < all_paths.size(); j++) {
                delete [] max_matching_path_arr[i][j];
            }
            delete [] max_matching_path_arr[i];
        }
        delete [] max_matching_path_arr;
    }
    if (min_matching_path_arr != NULL) {
        for (unsigned int i=0; i < all_paths.size(); i++) {
            for (unsigned int j=0; j < all_paths.size(); j++) {
                delete [] min_matching_path_arr[i][j];
            }
            delete [] min_matching_path_arr[i];
        }
        delete [] min_matching_path_arr;
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

int PopulationTree::max_matching_path(int path1, int path2, int t) const {
    /*    static int mult1 = all_paths.size() * model->ntimes;
    static int mult2 = model->ntimes;
    assert(path1 < (int)all_paths.size() && path1 >=0);
    assert(path2 < (int)all_paths.size() && path2 >=0);
    assert(t >= 0 && t < model->ntimes);
    return max_matching_path_arr[path1*mult1 + path2*mult2 + t];*/
    return max_matching_path_arr[path1][path2][t];
}

int PopulationTree::min_matching_path(int path1, int path2, int t) const {
    return min_matching_path_arr[path1][path2][t];
}


bool PopulationTree::paths_equal(int path1, int path2, int t1, int t2) const {
    if (path1 == path2) return true;
    if (t1 > model->ntimes - 1) t1 = model->ntimes - 1;
    if (t2 == -1 || t2 > model->ntimes - 1) t2 = model->ntimes - 1;
    assert(t1 <= t2);
    return ( max_matching_path(path1, path2, t1) >= t2 );
    //    return max_matching_path[path1][path2][t1] >= t2;
}


void UniquePath::update_prob(const vector<PopulationPath> &all_paths,
                             const vector<MigMatrix> &mig_matrix) {

    int p=first_path();
    assert(p >= 0 && p < (int)all_paths.size());
    prob = 1.0;
    num_mig=0;
    for (int t=start_time+1; t <= end_time; t++) {
        double thisprob = mig_matrix[2*t-1].get(all_paths[p].pop[t-1],
                                                all_paths[p].pop[t]);
        if (thisprob > 0 && thisprob < 0.5)
            num_mig++;
        prob *= thisprob;
    }
    assert(prob >= 0 && prob <= 1.0);
}

void PopulationTree::update_population_probs() {
    int ntime = model->ntimes;
    //    printf("update_population_probs\n");

    // update sub_paths probs
    for (int t1=0; t1 < ntime; t1++) {
        for (int t2=t1; t2 < ntime; t2++) {
            for (int p1=0; p1 < npop; p1++) {
                for (int p2=0; p2 < npop; p2++) {
                    sub_paths[t1][t2][p1][p2].update_probs(all_paths,
                                                           mig_matrix);
                }
            }
        }
    }

    for (unsigned int i=0; i < all_paths.size(); i++) {
        all_paths[i].prob = 1.0;
        for (int t=1; t < ntime; t++) {
            all_paths[i].prob *= mig_matrix[2*t-1].get(all_paths[i].pop[t-1],
                                                       all_paths[i].pop[t]);
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

    // check that all paths end in same population
    int lastpop=all_paths[0].pop[ntime - 1];
    for (unsigned int i=1; i < all_paths.size(); i++) {
        if (all_paths[i].get(ntime - 1) != lastpop)
            exitError("Error: populations do not converge by final time\n");
    }

    max_matching_path_arr = new int **[all_paths.size()];
    for (unsigned int i=0; i < all_paths.size(); i++) {
        max_matching_path_arr[i] = new int*[all_paths.size()];
        for (unsigned int j=0; j < all_paths.size(); j++) {
            max_matching_path_arr[i][j] = new int[ntime];
            for (int t=0; t < ntime; t++) {
                if (all_paths[i].get(t) != all_paths[j].get(t))
                    max_matching_path_arr[i][j][t] = -1;
                else {
                    max_matching_path_arr[i][j][t] = t;
                    for (int t1=t+1; t1 < model->ntimes; t1++) {
                        if (all_paths[i].get(t1) == all_paths[j].get(t1))
                            max_matching_path_arr[i][j][t] = t1;
                        else break;
                    }
                }
            }
        }
    }

    min_matching_path_arr = new int **[all_paths.size()];
    for (unsigned int i=0; i < all_paths.size(); i++) {
        min_matching_path_arr[i] = new int*[all_paths.size()];
        for (unsigned int j=0; j < all_paths.size(); j++) {
            min_matching_path_arr[i][j] = new int[ntime];
            for (int t=0; t < ntime; t++) {
                if (all_paths[i].get(t) != all_paths[j].get(t))
                    min_matching_path_arr[i][j][t] = -1;
                else {
                    min_matching_path_arr[i][j][t] = t;
                    for (int t1=t-1; t1 >= 0; t1--) {
                        if (all_paths[i].get(t1) == all_paths[j].get(t1))
                            min_matching_path_arr[i][j][t] = t1;
                        else break;
                    }
                }
            }
        }
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

                            // path i has this start/end.
                           // search previous paths for identical sub-path
                            int other = find_sub_path(i,
                                                      sub_paths[t1][t2][p1][p2],
                                                      t1, t2);
                            if (other >= 0) {
                                sub_paths[t1][t2][p1][p2].add_path_to_subpath(
                                     i, other);
                            } else {
                                sub_paths[t1][t2][p1][p2].new_subpath(t1, t2, i);
                            }
                        }
                    }
                }
            }
        }
    }


    update_population_probs();
}

 int PopulationTree::get_pop(int path, int time) const {
     if (time >= model->ntimes) return final_pop();
     return all_paths[path].get(time);
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
            assert(0);
            exitError("No consistent path found\n");
        return -1;
    }
    SubPath *possible_paths = &sub_paths[t1][t2][p1][p2];
    for (unsigned int i=0; i < possible_paths->size(); i++) {
        int path = possible_paths->first_path(i);
        if (paths_equal(path, path1, t1, t2)) {
            UniquePath *u = &possible_paths->unique_subs[i];
            for (set<int>::iterator it=u->path.begin(); it != u->path.end(); it++) {
                path = *it;
                assert(paths_equal(path, path1, t1, t2));
                if (paths_equal(path, path2, t2, t3))
                    return path;
            }
            break;
        }
    }
    if (require_exists) {
        assert(0);
        exitError("No consistent path found\n");
    }
    return -1;
}


 int PopulationTree::path_to_root(const LocalNode *nodes, int node, int time) const {
    assert(node != -1);
    int path = nodes[node].pop_path;
    int parent = nodes[node].parent;
    if (parent == -1) return path;
    int orig_age = (time < 0 ? nodes[node].age : time);
    assert(orig_age >= nodes[node].age);
    if (parent >= 0) assert(orig_age <= nodes[parent].age);
    while (true) {
        path = consistent_path(path,
                               nodes[parent].pop_path,
                               orig_age,
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



int PopulationTree::path_to_root(const spidir::Node *node, double time) const {
    assert(node != NULL);
    int path = node->pop_path;
    spidir::Node *parent = node->parent;
    if (parent == NULL) return path;
    int orig_age = (time < 0 ? node->age : time);
    assert(orig_age >= node->age);
    if (parent >= 0) assert(orig_age <= parent->age);
    while (true) {
        path = consistent_path(path,
                               parent->pop_path,
                               model->discretize_time(orig_age),
                               model->discretize_time(parent->age),
                               ( parent->parent == NULL ?
                                 model->ntimes - 1 :
                                 model->discretize_time(parent->parent->age)));
        node = parent;
        parent = node->parent;
        if (parent == NULL) break;
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
    pop_tree->update_population_probs();
}

}  //namespace argweaver
