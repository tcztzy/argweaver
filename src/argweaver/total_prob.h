#ifndef ARGWEAVER_TOTAL_PROB_H
#define ARGWEAVER_TOTAL_PROB_H

#include "local_tree.h"
#include "model.h"

namespace argweaver {

void calc_coal_rates_spr(const ArgModel *model, const LocalTree *tree,
                         const Spr &spr, const LineageCounts &lineages,
                         double *coal_rates);

double calc_spr_prob(const ArgModel *model, const LocalTree *tree,
                     const Spr &spr, LineageCounts &lineages,
                     double treelen=-1.0,
		     double **num_coal=NULL, double **num_nocoal=NULL,
                     double coal_weight=1.0, int lineages_counted=false);

double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees);

// NOTE: trees should be uncompressed and sequences compressed
double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees,
                           const SitesMapping* sites_mapping,
                           const TrackNullValue *maskmap_uncompressed);

double calc_arg_prior_recomb_integrate(const ArgModel *model,
                                       const LocalTrees *trees,
                                       double **num_coal=NULL,
                                       double **num_nocoal=NULL);

double calc_arg_prior(const ArgModel *model, const LocalTrees *trees,
                      double **num_coal=NULL, double **num_ncoal=NULL);
double calc_arg_joint_prob(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *trees);



} // namespace argweaver

#endif // ARGWEAVER_TOTAL_PROB_H
