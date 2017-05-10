#ifndef ARGWEAVER_EST_POPSIZE_H
#define ARGWEAVER_EST_POPSIZE_H

#include "local_tree.h"
#include "model.h"

namespace argweaver {


 struct popsize_data {
     double *arr_alloc;
     int arr_size;
     double ***coal_counts;
     double ***nocoal_counts;
     double *coal_totals;
     double *nocoal_totals;
     ArgModel *model;
     int popsize_idx;
     int numleaf;
     double t1, t2;
     int min_t, max_t;
 };

void est_popsize_local_trees(const ArgModel *model, const LocalTrees *trees,
                             double *popsizes);
void mle_popsize(ArgModel *model, const struct popsize_data *data, double min_total=0);
void mle_popsize(ArgModel *model, const LocalTrees *trees, double min_total=0);
void one_popsize_like_and_dlike(int t, double log_popsize, struct popsize_data *data,
				double *likelihood, double *dlikelihood, double *dlikelihood2=NULL);
double one_popsize_likelihood(int t, double log_popsize, struct popsize_data *data);

//note log_popsize should be length ntimes-1 (one for each whole time interval)
double popsize_likelihood(double *log_popsize, struct popsize_data *data);

double one_popsize_likelihood(int t, double log_popsize, struct popsize_data *data);
double one_popsize_dlikelihood(int t, double log_popsize, struct popsize_data *data);

void popsize_sufficient_stats(struct popsize_data *data, ArgModel *model, const LocalTrees *trees, bool add=false);
void delete_popsize_data(struct popsize_data *data);

void update_popsize_hmc(ArgModel *model, const LocalTrees *trees);
void set_data_time(struct popsize_data *data, int t);
void no_update_popsize(ArgModel *model, const LocalTrees *trees);


// functions for adaptive HMC; may want to move to separate file
double dotProduct(double *x, int len);
double kineticEnergy(double *x, int len);
//void leapFrogL(double *theta, double *r, double epsilon, int L,

// Metropolis-Hastings resampling; no longer used in favor of HMC
void resample_popsizes_mh(ArgModel *model, const LocalTrees *trees,
                          bool sample_popsize_recomb, double heat=1.0);

}; // namespace argweaver

#endif // ARGWEAVER_EST_POPSIZE_H
