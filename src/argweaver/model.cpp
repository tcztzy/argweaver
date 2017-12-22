#ifdef ARGWEAVER_MPI
#include "mpi.h"
#include "mcmcmc.h"
#endif

#include "model.h"
#include "pop_model.h"
#include "local_tree.h"
#include "MultiArray.h"
#include "logging.h"
#include "sequences.h"
#include "track.h"

namespace argweaver {


double get_delta_diff(double log_delta, const double *times, int ntimes, double maxtime) {
    double delta = exp(log_delta);
    int i=1;
    return get_time_point(i, ntimes-1, maxtime, delta) - times[i];
}

// estimate delta from time vector
// if delta cannot be found, return -1 (likely implies linear times)
double get_delta(const double *times, int ntimes, double maxtime) {
    double min_log_delta=-10, max_log_delta=10.0, tol=1e-10, mid_log_delta=0;
    double min_diff, max_diff, mid_diff;
    min_diff = get_delta_diff(min_log_delta, times, ntimes, maxtime);
    max_diff = get_delta_diff(max_log_delta, times, ntimes, maxtime);
    if (min_diff * max_diff >= 0) return -1.0;
    while (max_log_delta - min_log_delta > tol) {
	mid_diff = get_delta_diff(mid_log_delta, times, ntimes, maxtime);
	if (min_diff * mid_diff > 0) {
	    min_diff = mid_diff;
	    min_log_delta = mid_log_delta;
	} else {
            if (max_diff * mid_diff > 0) return -1.0;
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


bool PopTime::operator< ( const PopTime &other ) const {
    if (pop != other.pop) return ( pop < other.pop );
    return ( time < other.time );
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


int ArgModel::consistent_path(int path1, int path2,
                              double t1, double t2, double t3,
                              bool require_exists) const {
    int t1d = discretize_time(t1);
    int t2d = t2 < -0.1 ? -1 : discretize_time(t2, t1d);
    int t3d = t3 < -0.1 ? -1 : discretize_time(t3, t2d);
    return consistent_path(path1, path2, t1d, t2d, t3d, require_exists);
}

void ArgModel::copy(const ArgModel &other) {
    owned = true;
    rho = other.rho;
    mu = other.mu;
    infsites_penalty = other.infsites_penalty;
    unphased = other.unphased;
    unphased_file = other.unphased_file;
    popsize_config = other.popsize_config;
    mc3 = other.mc3;
    smc_prime = other.smc_prime;

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

// file should have format;
// pop time0 popsize0
// pop time1 popsize1
// pop time2 popsize2
// means that pop has popsize0 from time 0 to time0,
// popsize1 from time0 to time1, etc
// last entry should have time >= maxtime.
// only two columns required if only single population
void ArgModel::read_population_sizes(string popsize_file) {
    if (!popsizes)
        alloc_popsizes();
    FILE *infile = fopen(popsize_file.c_str(), "r");
    char *line = NULL;
    if (infile == NULL)
        exitError("error opening popsize file %s\n", popsize_file.c_str());
    int pop_idx[num_pops()];
    double next_time[num_pops()];
    for (int i=0; i < num_pops(); i++) {
        pop_idx[i]=0;
        next_time[i] = 0.0;
    }
    while (NULL != (line = fgetline(infile))) {
        chomp(line);
        if (line[0] == '#') {
            delete [] line;
            continue;
        }
        vector<string> splitStr;
        split(line, '\t', splitStr);
        int pop=-1;
        double curr_time=-1, curr_size=-1;
        if (num_pops() == 2 && splitStr.size() == 2) {
            pop = 0;
            curr_time = atof(splitStr[0].c_str());
            curr_size = atof(splitStr[1].c_str());
        } else if (splitStr.size() == 3) {
            pop = atoi(splitStr[0].c_str());
            if (pop < 0 || pop >= num_pops())
                exitError("Error parsing population in popsize file\n");
            curr_time = atof(splitStr[1].c_str());
            curr_size = atof(splitStr[2].c_str());
        } else exitError("Error reading popsize file; format should be pop, time, size");
        if (curr_time < next_time[pop])
            exitError("Error reading popsize file; times should be increasing");
        assert(curr_size > 0);
        assert(pop >= 0 && pop < num_pops());
        assert(pop_idx[pop] < 2*ntimes - 2);
        while (curr_time > next_time[pop]) {
            popsizes[pop][pop_idx[pop]] = curr_size;
            next_time[pop] += coal_time_steps[pop_idx[pop]++];
            if (pop_idx[pop] >= 2*ntimes - 1) break;
        }
        delete [] line;
    }
    fclose(infile);
    for (int pop=0; pop < num_pops(); pop++)
        for (int i=0; i < 2*ntimes - 2; i++) {
            if (popsizes[pop][i] == 0.0)
                exitError("Error in read_population_sizes: some population sizes are zero or not set");
        }
}


void ArgModel::read_population_tree(FILE *infile) {
    if (pop_tree != NULL)
        delete pop_tree;
    pop_tree = new PopulationTree(1, this);
    argweaver::read_population_tree(infile, pop_tree);
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

int ArgModel::path_to_root(const LocalNode *nodes, int node, int time) const {
    if (pop_tree == NULL) return 0;
    return pop_tree->path_to_root(nodes, node, time);
}

int ArgModel::path_to_root(const spidir::Node *node, double time) const {
    if (pop_tree == NULL) return 0;
    return pop_tree->path_to_root(node, time);
}

ArgModel::ArgModel(const char *logfilename) {
    FILE *logfile = fopen(logfilename, "r");
    char *line = NULL;
    bool found_model = false;
    int npop=1;
    int numpath;
    char *pop_file=NULL;
    bool read_pop_file = false;
    pop_tree = NULL;
    smc_prime=true;
    if (logfile == NULL) {
        printError("Could not open log file %s\n", logfilename);
        abort();
    }
    times=NULL;
    while (NULL != (line = fgetline(logfile))) {
        chomp(line);
        if (str_starts_with(line, "command:")) {
            vector<string> splitStr;
            split(line, ' ', splitStr);
            for (int i=1; i < (int)splitStr.size()-1; i++) {
                if (strcmp(splitStr[i].c_str(), "--pop-tree-file")==0) {
                    pop_file = new char[strlen(splitStr[i+1].c_str())+1];
                    strcpy(pop_file, splitStr[i+1].c_str());
                    break;
                }
            }
        }
        if (str_starts_with(line, "model:")) found_model = true;
        if (found_model) {
            if (str_starts_with(line, "----------")) break;
            if (str_starts_with(line, "  mu = "))
                assert(1 == sscanf(line, "  mu = %le", &mu));
            if (str_starts_with(line, "  rho = "))
                assert(1 == sscanf(line, "  rho = %le", &rho));
            if (str_starts_with(line, "  ntimes = ")) {
                assert(1 == sscanf(line, "  ntimes = %i", &ntimes));
                assert(ntimes > 0);
                times = new double [ntimes];
            }
            if (str_starts_with(line, "  times = [")) {
                assert(times != NULL);
                assert(ntimes > 0);
                vector<string> splitStr;
                assert(line[strlen(line)-1] == ']');
                line[strlen(line)-1] = '\0';
                split(&line[11], ',', splitStr);
                assert((int)splitStr.size() == ntimes);
                for (int i=0; i < ntimes; i++)
                    times[i] = atof(splitStr[i].c_str());
                double delta = get_delta(times, ntimes, times[ntimes-1]);
                // if delta < 0 then using linear steps
                coal_time_steps = new double[2*ntimes];
                get_coal_time_steps(times, ntimes, coal_time_steps, delta < 0, delta);
            }
            if (str_starts_with(line, "  npop = ")) {
                if (pop_file != NULL) {
                    FILE *infile = fopen(pop_file, "r");
                    if (infile == NULL) {
                        fprintf(stderr,
                                "Warning: Could not open pop_file %s. "
                                "Migration probabilities will not be set correctyl!\n",
                                pop_file);
                    } else {
                        read_population_tree(infile);
                        fclose(infile);
                        read_pop_file = true;
                        npop = this->num_pops();
                        int temp_npop;
                        assert(1 == sscanf(line, " npop = %i", &temp_npop));
                        assert(temp_npop == npop);
                    }

                }
                if (!read_pop_file) {
                    assert(1 == sscanf(line, "  npop = %i", &npop));
                    assert(npop >= 1);
                    if (npop > 1) pop_tree = new PopulationTree(npop, this);
                }
            }
            if (str_starts_with(line, "  popsizes = [")) {
                alloc_popsizes();
                char *tmpstr = &line[14];
                vector<double> vals;
                bool done=false;
                while (true) {
                    assert(line[strlen(line)-1] == ',' ||
                           line[strlen(line)-1] == ']');
                    if (line[strlen(line)-1] == ']') done=true;
                    vector<string> splitStr;
                    tmpstr = trim(tmpstr);
                    split(tmpstr, ',', splitStr);
                    for (int i=0; i < (int)splitStr.size(); i++) {
                        if (splitStr[i].size() > 0)
                            vals.push_back(atof(splitStr[i].c_str()));
                    }
                    if (done) break;
                    assert(NULL != (line=fgetline(logfile)));
                    chomp(line);
                }
                assert((int)vals.size() == npop * ntimes ||
                       (int)vals.size() == npop * (2*ntimes - 1));
                if ((int)vals.size() == npop * ntimes) {
                    int idx=0;
                    for (int t=0; t < ntimes; t++) {
                        for (int pop=0; pop < npop; pop++) {
                            popsizes[pop][2*t] = vals[idx];
                            if (t != 0)
                                popsizes[pop][2*t-1] = vals[idx];
                            idx++;
                        }
                    }
                } else if ((int)vals.size() == npop * (2*ntimes - 1)) {
                    int idx=0;
                    for (int t=0; t < 2*ntimes - 1; t++)
                        for (int pop=0; pop < npop; pop++)
                            popsizes[pop][t] = vals[idx++];
                } else {
                    fprintf(stderr, "Error reading model from log file; wrong length of popsize vector");
                    assert(0);
                }
            }
            if (str_starts_with(line, "    numpath = ") && !read_pop_file) {
                assert(1 == sscanf(line, "    numpath = %i", &numpath));
                assert(pop_tree != NULL);
                MultiArray path(2, numpath, ntimes);
                for (int i=0; i < numpath; i++) {
                    char tmpstr[100];
                    vector<string> splitStr;
                    sprintf(tmpstr, "    path%i = [", i);
                    assert(NULL != (line=fgetline(logfile)));
                    chomp(line);
                    assert(str_starts_with(line, tmpstr));
                    assert(line[strlen(line)-1] == ']');
                    line[strlen(line)-1] = '\0';
                    split(&line[strlen(tmpstr)], ',', splitStr);
                    assert((int)splitStr.size() == ntimes);
                    for (int j=0; j < ntimes; j++) {
                        int p = atoi(splitStr[j].c_str());
                        assert(p >= 0 && p < npop);
                        path.set((double)p, i, j);
                    }
                }
                for (int t1=0; t1 < (ntimes-1); t1++) {
                    int t2=t1+1;
                    // this is not efficient but doesn't matter
                    for (int from_pop=0; from_pop < npop; from_pop++) {
                        set<int> to_pop;
                        for (int p=0; p < numpath; p++) {
                            if ((int)path.get(p, t1) == from_pop) {
                                to_pop.insert((int)path.get(p, t2));
                            }
                        }
                        double migprob = 0;
                        bool has_self = to_pop.find(from_pop) != to_pop.end();
                        // we don't have the actual migration/divergence
                        // probabilities in the log file, we are really just
                        // after the population tree structure. So set
                        // migration probability to 0.1 and divergence
                        // probability evenly spread across possible populations
                        if (has_self && to_pop.size() > 1)
                            migprob = 0.1/(to_pop.size()-1); //migration
                        else if (!has_self)
                            migprob = 1.0/(to_pop.size());  //divergence
                        if (migprob > 0.0) {
                            for (set<int>::iterator it=to_pop.begin();
                                 it!=to_pop.end(); it++) {
                                int p2 = *it;
                                if (p2 != from_pop)
                                    pop_tree->add_migration(2*t1+1, from_pop, p2, migprob);
                            }
                        }
                    }
                }
                pop_tree->set_up_population_paths();
                pop_tree->update_population_probs();
            }
        }
    }
    if (pop_file != NULL) delete[] pop_file;
}


void ArgModel::log_model() const {
    printLog(LOG_LOW, "\n");
    printLog(LOG_LOW, "model: \n");
    printLog(LOG_LOW, "  mu = %e\n", mu);
    printLog(LOG_LOW, "  rho = %e\n", rho);
    printLog(LOG_LOW, "  ntimes = %d\n", ntimes);
    printLog(LOG_LOW, "  times = [");
    for (int i=0; i<ntimes-1; i++)
        printLog(LOG_LOW, "%f,", times[i]);
    printLog(LOG_LOW, "%f]\n", times[ntimes-1]);
    printLog(LOG_LOW, "  npop = %d\n", num_pops());
    printLog(LOG_LOW, "  popsizes = [");
    for (int i=0; i < 2*ntimes-1; i++) {
        if (i != 0) printLog(LOG_LOW, "              ");
        for (int pop = 0; pop < num_pops(); pop++) {
            if (pop != 0) printLog(LOG_LOW, ",\t");
            printLog(LOG_LOW, "%.1f", popsizes[pop][i]);
        }
        printLog(LOG_LOW, "%c", i == 2*ntimes-2 ? ']' : ',');
        printLog(LOG_LOW, "\n");
    }
    if (pop_tree != NULL) {
        printLog(LOG_LOW, "    numpath = %d\n", num_pop_paths());
        for (int i=0; i < num_pop_paths(); i++) {
            printLog(LOG_LOW, "    path%d = [%d", i,
                   get_pop(i, 0));
            for (int j=1; j < ntimes; j++)
                printLog(LOG_LOW, ", %d", get_pop(i, j));
            printLog(LOG_LOW, "]\n");
        }
        printLog(LOG_LOW, " max_migrations = %i\n", pop_tree->max_migrations);
    }
    if (isLogLevel(LOG_HIGH)) {
        printLog(LOG_HIGH, "mutmap = [\n");
        for (unsigned int i=0; i<mutmap.size(); i++) {
            printLog(LOG_HIGH, "%d\t%d\t%e\n",
                     mutmap[i].start, mutmap[i].end,
                     mutmap[i].value);
        }
        printLog(LOG_HIGH, "]\n");

        printLog(LOG_HIGH, "recombmap = [\n");
        for (unsigned int i=0; i<recombmap.size(); i++) {
            printLog(LOG_HIGH, "%d\t%d\t%e\n",
                     recombmap[i].start, recombmap[i].end,
                     recombmap[i].value);
        }
        printLog(LOG_HIGH, "]\n");
    }

    printLog(LOG_LOW, "\n");
}

int time_index(double t, const double *times, int ntimes, int min_idx,
               double tol) {
    int min_time= (min_idx < 0 ? 0 : min_idx);
    int max_time = ntimes - 1;
    int mid_time = (max_time+min_time)/2;
    // check min_time first as it may often be the one
    if (fabs(t - times[min_time]) < tol) return min_time;
    assert(t > times[min_time]);
    // otherwise do binary search
    while (1) {
        if (fabs(t - times[mid_time]) < tol) return mid_time;
        if (times[mid_time] > t) {
            max_time = mid_time-1;
        } else {
            min_time = mid_time+1;
        }
        mid_time = (max_time + min_time)/2;
        if (max_time <= min_time) {
            if (fabs(t - times[max_time]) < tol) return max_time;
            if (fabs(t - times[min_time]) < tol) return min_time;
            assert(0);
        }
    }
    assert(0);
}


int ArgModel::discretize_time(double t, int min_idx, double tol) const {
    return time_index(t, times, ntimes, min_idx, tol);
}


// return iteration
int ArgModel::init_params_from_statfile(vector<string> header, const char *line,
                                         int popsize_em) {
    // parse stage and last iter
    vector<string> tokens;
    split(line, "\t", tokens);
    if (tokens.size() < 2) {
        printError("incomplete line in status file");
        return false;
    }
    assert(tokens.size() == header.size());

    string stage2 = tokens[0];
    int iter;
    if (sscanf(tokens[1].c_str(), "%d", &iter) != 1) {
        printError("iter column is not an integer");
        return false;
    }

    // NOTE: only resume resample stage for now
    if (stage2 != "resample")
        return -1;

    int npop = num_pops();
    if (popsize_config.sample) {
        list<PopsizeConfigParam> l=popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin(); it != l.end();
             ++it) {
            string popname=it->name;
            set<PopTime> popset = it->intervals;
            int found=false;
            for (unsigned int i=0; i < header.size(); i++) {
                if (header[i] == popname) {
                    found=true;
                    double tempN;
                    sscanf(tokens[i].c_str(), "%lf", &tempN);
                    for (set<PopTime>::iterator it2=popset.begin();
                         it2 != popset.end(); ++it2) {
                        int pop = it2->pop;
                        int time = it2->time;
                        if (pop < 0 || pop >= npop ||
                            time < 0 || time >= 2*ntimes-1) {
                            printError("Error in resume: popsize config does"
                                       " not match previous run\n");
                            abort();
                        }
                        popsizes[pop][time] = tempN;
                    }
                    break;
                }
            }
            if (!found) {
                printError("Error in resume: did not find pop %s in previous"
                           " run stats file\n",
                           popname.c_str());
                abort();
            }
        }
    } else if (popsize_em) {
        for (int pop=0; pop < npop; pop++)
            for (int i=0; i < ntimes - 1; i++)
                sscanf(tokens[i].c_str(), "%lf", &popsizes[pop][i]);
    }

    if (pop_tree != NULL) {
        for (unsigned int i=0; i < pop_tree->mig_params.size(); i++) {
            MigParam mp = pop_tree->mig_params[i];
            double migrate;
            bool found=false;
            for (unsigned int j=0; j < header.size(); j++) {
                if (mp.name == header[j]) {
                    found=true;
                    sscanf(tokens[j].c_str(), "%lf", &migrate);
                    pop_tree->mig_matrix[mp.time_idx].set(mp.from_pop,
                                                          mp.to_pop,
                                                          migrate);
                    double self_mig=1.0;
                    for (int k=0; k < npop; k++) {
                        if (k != mp.from_pop)
                            self_mig -= pop_tree->mig_matrix[mp.time_idx].get(mp.from_pop, k);
                    }
                    assert(self_mig > 0 && self_mig <= 1.0);
                    pop_tree->mig_matrix[mp.time_idx].set(mp.from_pop, mp.from_pop, self_mig);
                    pop_tree->update_population_probs();
                    break;
                }
            }
            if (!found) {
                printError("Error resuming run; could not find mig rate estimate %s in stats file\n",
                           mp.name.c_str());
                abort();
            }
        }
    }
    return iter;
}


void uncompress_model(ArgModel *model, const SitesMapping *sites_mapping,
		      double compress_seq)
{
    model->rho /= compress_seq;
    model->mu /= compress_seq;

    uncompress_track(model->mutmap, sites_mapping, compress_seq, true);
    uncompress_track(model->recombmap, sites_mapping, compress_seq, true);
}


void compress_model(ArgModel *model, const SitesMapping *sites_mapping,
                    double compress_seq)
{
    model->rho *= compress_seq;
    model->mu *= compress_seq;

    compress_track(model->mutmap, sites_mapping, compress_seq, true);
    compress_track(model->recombmap, sites_mapping, compress_seq, true);
}

} // namespace argweaver
