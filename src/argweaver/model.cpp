#ifdef ARGWEAVER_MPI
#include "mpi.h"
#include "mcmcmc.h"
#endif

#include "model.h"
#include "pop_model.h"


namespace argweaver {


double get_delta_diff(double log_delta, const double *times, int ntimes, double maxtime) {
    double delta = exp(log_delta);
    int i=1;
    return get_time_point(i, ntimes-1, maxtime, delta) - times[i];
}

double get_delta(const double *times, int ntimes, double maxtime) {
    double min_log_delta=-10, max_log_delta=10.0, tol=1e-10, mid_log_delta=0;
    double min_diff, max_diff, mid_diff;
    min_diff = get_delta_diff(min_log_delta, times, ntimes, maxtime);
    max_diff = get_delta_diff(max_log_delta, times, ntimes, maxtime);
    assert(min_diff * max_diff < 0.0);
    while (max_log_delta - min_log_delta > tol) {
	mid_diff = get_delta_diff(mid_log_delta, times, ntimes, maxtime);
	if (min_diff * mid_diff > 0) {
	    min_diff = mid_diff;
	    min_log_delta = mid_log_delta;
	} else {
	    assert(max_diff * mid_diff > 0);
	    max_diff = mid_diff;
	    max_log_delta = mid_log_delta;
	}
	mid_log_delta = (min_log_delta + max_log_delta)/2.0;
    }
    double delta = exp(mid_log_delta);
    printLog(LOG_LOW, "using delta=%e\n", delta);
    return delta;
}

void get_coal_time_steps(const double *times, int ntimes,
			 double *coal_time_steps, bool linear,
			 double delta)
{
    // get midpoints
    double times2[2*ntimes+1];
    for (int i=0; i < ntimes; i++)
	times2[2*i] = times[i];
    if (linear) {
	for (int i=0; i < ntimes-1; i++)
	    times2[2*i+1] = 0.5*(times[i+1] + times[i]);
    } else {
	for (int i=0; i < ntimes-1; i++)
	    times2[2*i+1] = get_time_point(2*i+1, 2*ntimes-2, times[ntimes-1],
					   delta);
    }
    for (int i=0; i < 2*ntimes-2; i++) {
	coal_time_steps[i] = times2[min(i+1, 2*ntimes)] - times2[i];
	if (coal_time_steps[i] < 0) {
	    assert(0);
	}
    }
    coal_time_steps[2*ntimes-2] = INFINITY;
}



// returns true if regions in track are flush with one another
template <class T>
bool check_map(const Track<T> &track, int start, int end)
{
    // check that track is not empty
    typename Track<T>::const_iterator it = track.begin();
    if (it == track.end()) {
        printError("map is empty");
        return false;
    }

    // check that start and end cover desired range
    if (it->start > start || track.back().end < end) {
        printError("map does not cover entire region");
        return false;
    }

    int last = it->end;
    ++it;
    for (; it != track.end(); ++it) {
        if (it->start != last) {
            printError("map is not complete at %s:%d",
                       it->chrom.c_str(), it->start);
            return false;
        }
        last = it->end;
    }

    return true;
}


// returns true if regions in track are flush with one another
template <class T>
bool complete_map(Track<T> &track, string chrom, int start, int end, const T &default_value)
{
    // check for empty track
    if (track.size() == 0) {
        track.append(chrom, start, end, default_value);
        return true;
    }

    // check that start and end cover desired range
    typename Track<T>::iterator it = track.begin();
    if (it->start > start)
        track.insert(it, RegionValue<T>(chrom, start, it->start, default_value));
    if (track.back().end < end)
        track.append(chrom, track.back().end, end, default_value);

    it = track.begin();
    int last = it->end;
    ++it;
    for (; it != track.end(); ++it) {
        if (it->start > last) {
            it = track.insert(it, RegionValue<T>(chrom, last, it->start, default_value));
        } else if (it->start < last) {
            printError("map contains over laps %s:%d-%d",
                       chrom.c_str(), it->start, last);
            return false;
        }
        last = it->end;
    }

    return true;
}

int ArgModel::get_pop(int path, int time) const {
    if (pop_tree == NULL) return 0;
    if (time >= ntimes)
        time = ntimes-1;
    return pop_tree->get_pop(path, time);
}


int ArgModel::consistent_path(int path1, int path2,
                              int t1, int t2, int t3,
                              bool require_exists) const {
    if (pop_tree == NULL) return 0;
    return pop_tree->consistent_path(path1, path2, t1, t2, t3,
                                     require_exists);
}


void ArgModel::copy(const ArgModel &other) {
    owned = true;
    rho = other.rho;
    mu = other.mu;
    infsites_penalty = other.infsites_penalty;
    unphased = other.unphased;
    sample_phase = other.sample_phase;
    unphased_file = other.unphased_file;
    popsize_config = other.popsize_config;
    mc3 = other.mc3;

    if (other.pop_tree)
        pop_tree = new PopulationTree(*other.pop_tree);
    else pop_tree = NULL;

    // copy popsizes and times
    set_times(other.times, other.coal_time_steps, ntimes);
    if (other.popsizes)
        set_popsizes(other.popsizes);

    // copy maps
    if (other.mutmap.size() > 0)
        mutmap.insert(mutmap.begin(),
                      other.mutmap.begin(), other.mutmap.end());
    if (other.recombmap.size() > 0)
        recombmap.insert(recombmap.begin(),
                         other.recombmap.begin(), other.recombmap.end());
}

void ArgModel::clear() {
    if (owned) {
        delete [] times;
        delete [] time_steps;
        delete [] coal_time_steps;
        if (popsizes) {
            int npop = this->num_pops();
            for (int i=0; i < npop; i++)
                delete [] popsizes[i];
            delete [] popsizes;
        }
        if (pop_tree)
            delete pop_tree;
    }
}


// Initializes mutation and recombination maps for use
bool ArgModel::setup_maps(string chrom, int start, int end) {

    // check maps
    if (!complete_map(mutmap, chrom, start, end, mu)) {
        printError("mutation map has errors");
        return false;
    }
    if (!complete_map(recombmap, chrom, start, end, rho)) {
        printError("recombination map has errors");
        return false;
    }

    /*
    // setup default maps
    if (mutmap.size() == 0)
        mutmap.append(chrom, start, end, mu);
    if (recombmap.size() == 0)
        recombmap.append(chrom, start, end, rho);

    // check maps
    if (!check_map(mutmap, start, end)) {
        printError("mutation map has errors");
        return false;
    }
    if (!check_map(recombmap, start, end)) {
        printError("recombination map has errors");
        return false;
    }
    */

    // create new mut and recomb maps that share common boundaries
    int pos = start, pos2;
    unsigned int i = 0;
    unsigned int j = 0;
    Track<double> mutmap2;
    Track<double> recombmap2;
    while (i < mutmap.size() || j < recombmap.size()) {
        //printf("recomb[%d] = (%d, %d, %e), mut[%d] = (%d, %d, %e)\n",
        //       j, recombmap[j].start, recombmap[j].end, recombmap[j].value,
        //       i, mutmap[i].start, mutmap[i].end, mutmap[i].value);

        if (mutmap[i].end < recombmap[j].end) {
            // mutation region ends first
            pos2 = mutmap[i].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            i++;
        } else if (mutmap[i].end > recombmap[j].end) {
            // recombination region ends first
            pos2 = recombmap[j].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            j++;
        } else {
            // mutation and region region end together
            pos2 = recombmap[j].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            i++;
            j++;
        }
    }

    // copy over new maps
    mutmap.clear();
    recombmap.clear();
    mutmap.insert(mutmap.begin(), mutmap2.begin(), mutmap2.end());
    recombmap.insert(recombmap.begin(), recombmap2.begin(), recombmap2.end());

    return true;
}


void ArgModel::set_popsizes_random(double popsize_min,
                                   double popsize_max) {
    if (!popsizes)
        alloc_popsizes();
    int npop = this->num_pops();
#ifdef ARGWEAVER_MPI
    if (mc3.group_comm->Get_rank() == 0) {
#endif
        if (popsize_config.size() == 0) {
            for (int pop=0; pop < npop; pop++) {
                for (int i=0; i < 2*ntimes-1; i++) {
                    popsizes[pop][i] = frand(popsize_min, popsize_max);
                }
            }
            return;
        }
        list<PopsizeConfigParam> l = popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin();
             it != l.end(); ++it) {
            double popsize=frand(popsize_min, popsize_max);
            for (set<PopTime>::iterator it2 = it->intervals.begin();
                 it2 != it->intervals.end(); ++it2)
                popsizes[it2->pop][it2->time] = popsize;
        }
#ifdef ARGWEAVER_MPI
    }
    for (int pop=0; pop < npop; pop++) {
        mc3.group_comm->Bcast(popsizes[pop], 2*ntimes-1, MPI::DOUBLE, 0);
    }
#endif
}


void PopsizeConfig::addInterval(const char *name, int pop,
                                int time, int sample) {
    for (list<PopsizeConfigParam>::iterator it=params.begin();
	 it != params.end(); ++it) {
	if (strcmp(it->name.c_str(), name)==0) {
	    it->add_interval(pop, time);
	    if (it->sample != sample) {
		printError("Error in PopsizeConfig.add: got conflicting info"
                           " on whether to sample pop %s\n", name);
		exit(1);
	    }
	    return;
	}
    }
    // need to add new param
    PopsizeConfigParam p(string(name), sample, pop, time);
    params.push_back(p);
}

PopsizeConfig::PopsizeConfig(string filename, int ntimes, int npop,
                             double **popsizes) :
    sample(true),
    popsize_prior_alpha(1.0),
    popsize_prior_beta(1.0e-4),
    config_buildup(0),
    epsilon(0.01),
    pseudocount(0)
  {
    if (filename=="") {
        for (int pop=0; pop < npop; pop++) {
            for (int i=0; i < ntimes; i++) {
                char str[1000];
                sprintf(str, "N%d.%d", pop, i);
                if (i > 0)
                    addInterval(str, pop, 2*i-1, true);
                addInterval(str, pop, 2*i, true);
            }
	}
    } else {
	FILE *infile = fopen(filename.c_str(), "r");
	if (infile == NULL) {
	    fprintf(stderr, "Error opening popsize config file %s\n",
		    filename.c_str());
	    exit(0);
	}
	char *line;
        while (NULL != (line = fgetline(infile))) {
	    chomp(line);
	    vector<string>tokens;
	    split(line, '\t', tokens);
	    string intervalname = tokens[0];
            int pop=0;
            bool sample=true;
            if (tokens.size() < 2)
                exitError("Expect at least two cols on each line in popsize config file: param_name, time_idx [, pop, sample, init_val]");
            int time = atoi(tokens[1].c_str());
            if (time < 0 || time > 2*ntimes-1)
                exitError("time index out of range [0, 2*ntime-1] in popsize config file\n");
            if (tokens.size() >= 3) {
                pop = atoi(tokens[2].c_str());
                if (pop < 0 || pop >= npop)
                    exitError("pop out of range [0, npop-1] in popsize config file\n");
            }
            if (tokens.size() >= 4)
                sample = (atoi(tokens[3].c_str()) != 0);
            if (tokens.size() == 5)
                popsizes[pop][time] = atof(tokens[4].c_str());
            if (tokens.size() > 5)
                exitError("Too many columns in popsize config file; maximum is 5 (param_name, time_idx, pop, sample, init_val\n");
	    addInterval(tokens[0].c_str(), pop, time, sample);
	    delete [] line;
	}
	while (NULL != (line = fgetline(infile))) {
	    char *tline = trim(line);
	    if (strlen(tline) != 0) {
		fprintf(stderr, "Error: too many lines in popsize config file"
                        " %s; expected ntimes=%i entries.\n",
                        filename.c_str(), 2*ntimes-1);
		exit(0);
	    }
	}
    }
    printf("done set_popsize_config num_n_params=%i\n", (int)params.size());
}


void ArgModel::set_popsize_config_by_pop_tree() {
    if (pop_tree == NULL)
        exitError("Error: no population tree defined in set_popsize_config_by_pop_tree\n");
    int npop = this->num_pops();
    int nextpop = npop;
    char tmp[1000];
    popsize_config = PopsizeConfig();
    popsize_config.sample = true;
    for (int pop=0; pop < npop; pop++) {
        sprintf(tmp, "N%i", pop);

        for (int t=0; t < ntimes-1; t++) {

            popsize_config.addInterval(tmp, pop, 2*t, true);

            MigMatrix *migmat = &(pop_tree->mig_matrix[t]);
            if (migmat->get(pop, pop) == 0) break;  // pop becomes extinct

            for (int pop2=0; pop2 < npop; pop2++) {
                if (pop != pop2 && migmat->get(pop2, pop) == 1) {
                    sprintf(tmp, "N%i", nextpop++);
                    break;
                }
            }
            popsize_config.addInterval(tmp, pop, 2*t+1, true);
        }
    }
    printf("Done set_popsize_config_by_pop_tree numParam = %i\n", nextpop);
}

void ArgModel::read_population_tree(string pop_file) {
    if (pop_tree != NULL)
        delete pop_tree;
    pop_tree = new PopulationTree(1, this);
    FILE *infile = fopen(pop_file.c_str(), "r");
    if (infile == NULL)
        exitError("error opening population file %s\n", pop_file.c_str());
    argweaver::read_population_tree(infile, pop_tree);
    fclose(infile);
}


int ArgModel::num_pops() const {
    if (pop_tree == NULL) return 1;
    return pop_tree->npop;
}

int ArgModel::num_pop_paths() const {
    if (pop_tree == NULL) return 1;
    return pop_tree->num_pop_paths();
}

double ArgModel::path_prob(int path, int t1, int t2) const {
    if (pop_tree == NULL) return 1.0;
    return pop_tree->path_prob(path, t1, t2);
}

bool ArgModel::paths_equal(int path1, int path2, int t1, int t2) const {
    if (pop_tree == NULL || path1==path2) return true;
    return pop_tree->paths_equal(path1, path2, t1, t2);
}

int ArgModel::max_matching_path(int path1, int path2, int t) const {
    if (pop_tree == NULL) return ntimes-1;
    return pop_tree->max_matching_path(path1, path2, t);
}


} // namespace argweaver
