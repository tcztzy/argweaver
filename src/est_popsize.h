#ifndef ARGWEAVER_EST_POPSIZE_H
#define ARGWEAVER_EST_POPSIZE_H

#include "local_tree.h"
#include "model.h"

namespace argweaver {

void resample_popsizes(ArgModel *model, const LocalTrees *trees,
                       bool sample_popsize_recomb, double heat=1.0);

void est_popsize_local_trees(const ArgModel *model, const LocalTrees *trees,
                             double *popsizes);

} // namespace argweaver

#endif // ARGWEAVER_EST_POPSIZE_H
