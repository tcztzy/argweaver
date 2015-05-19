#ifdef ARGWEAVER_MPI
#include "mpi.h"
#endif

#include <algorithm>
#include <vector>
#include <map>

#include "est_popsize.h"
#include "total_prob.h"
#include "logging.h"
#include "model.h"


namespace argweaver {

using namespace std;


double resample_popsize_scale(ArgModel *model, const LocalTrees *trees,
                              bool sample_popsize_recomb, double heat,
                              double curr_like) {
    /* Prior for scale is gamma distribution w theta = 200000, k=1.
       This has mean 200,000 but basically flat from 0 to a million
     */
    list<PopsizeConfigParam> &l = model->popsize_config.params;
    double scale;

#ifdef ARGWEAVER_MPI
   MPI::Intracomm *comm = model->mc3.group_comm;
   int rank = comm->Get_rank();
   if (rank == 0) {
#endif
       double trans_k = 1000.0;
       double trans_theta = 0.001;
       scale = rand_gamma(trans_k, trans_theta);

       double trans_ratio = -2.0 * (trans_k - 1.0) * log(scale)
           + scale / trans_theta - 1.0/(scale * trans_theta);

       double popsize_sum=0.0;
       double prior_ratio = 0.0;
       double prior_theta = 200000.0;
       //       double prior_k = 1.0;
       for (list<PopsizeConfigParam>::iterator it=l.begin(); it != l.end();
            it++) {
           set<int>::iterator it2 = it->pops.begin();
           prior_ratio += model->popsizes[*it2] / prior_theta * (1.0 - scale);
           popsize_sum += model->popsizes[*it2];
       }

#ifdef ARGWEAVER_MPI
       comm->Bcast(&scale, 1, MPI::DOUBLE, 0);
#endif

       for (int i=0; i < 2*model->ntimes-1; i++)
           model->popsizes[i] = model->popsizes[i] * scale;
       double new_like = sample_popsize_recomb ? calc_arg_prior(model, trees) :
           calc_arg_prior_recomb_integrate(model, trees);
#ifdef ARGWEAVER_MPI
       comm->Reduce(MPI_IN_PLACE, &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
#endif
       double lr = new_like - curr_like;
       double ln_accept = trans_ratio + prior_ratio + lr;
       ln_accept *= heat;
       double pr_accept = (ln_accept > 0 ? 1.0 : exp(ln_accept));
       bool accept = (ln_accept > 0 || frand() < pr_accept);
#ifdef ARGWEAVER_MPI
       comm->Bcast(&accept, 1, MPI::BOOL, 0);
#endif
       if (!accept) {
           for (int i=0; i < 2*model->ntimes-1; i++)
               model->popsizes[i] = model->popsizes[i] / scale;
       } else {
           curr_like = new_like;
       }
       printLog(LOG_LOW, "sample_scale\t%f\t%f\t%f\t%f\t%f\t%s\n",
                popsize_sum, scale, lr, trans_ratio, prior_ratio,
                accept ? "accept" : "reject");

#ifdef ARGWEAVER_MPI
   } else {
       comm->Bcast(&scale, 1, MPI::DOUBLE, 0);
       for (int i=0; i < 2*model->ntimes-1; i++)
           model->popsizes[i]  = model->popsizes[i] * scale;
       double new_like = sample_popsize_recomb ? calc_arg_prior(model, trees) :
           calc_arg_prior_recomb_integrate(model, trees);
       comm->Reduce(&new_like, &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
       bool accept;
       comm->Bcast(&accept, 1, MPI::BOOL, 0);
       if (!accept) {
           for (int i=0; i < 2*model->ntimes-1; i++)
               model->popsizes[i] = model->popsizes[i] / scale;
       } else {
           curr_like = new_like;
       }
   }
#endif
   return curr_like;
}


double resample_single_popsize(ArgModel *model, const LocalTrees *trees,
                               bool sample_popsize_recomb, double heat,
                               const list<PopsizeConfigParam>::iterator &it,
                               double curr_like, int index) {
    list<PopsizeConfigParam> &l = model->popsize_config.params;
    double adjust[2];
    bool accept;
    if (l.size()==1) return curr_like;

#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    if (rank == 0) {
#endif

        adjust[0] = 1000.0 * frand() - 500.0;
        adjust[1] = -adjust[0]/(double)(l.size()-1);

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(&adjust, 2, MPI::DOUBLE, 0);
#endif

    for (list<PopsizeConfigParam>::iterator it2=l.begin();
         it2 != l.end(); it2++) {
        set<int>::iterator it3 = it2->pops.begin();
        if (model->popsizes[*it3] + adjust[it2 == it ? 0 : 1] <= 0.0)
            return curr_like;
    }
    for (list<PopsizeConfigParam>::iterator it2=l.begin();
         it2 != l.end(); it2++) {
        for (set<int>::iterator it3=it2->pops.begin();
             it3 != it2->pops.end(); it3++) {
            model->popsizes[*it3] += adjust[it2 == it ? 0 : 1];
        }
    }
    double *num_coal = (double*)malloc(2 * model->ntimes * sizeof(double));
    double *num_nocoal = (double*)malloc(2 * model->ntimes * sizeof(double));
    double new_like = sample_popsize_recomb ?
        calc_arg_prior(model, trees, num_coal, num_nocoal) :
        calc_arg_prior_recomb_integrate(model, trees, num_coal, num_nocoal);

#ifdef ARGWEAVER_MPI
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : &new_like,
                 &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : num_coal,
                 num_coal, 2 * model->ntimes - 1, MPI::DOUBLE, MPI::SUM, 0);
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : num_nocoal,
                 num_nocoal, 2 * model->ntimes - 1, MPI::DOUBLE, MPI::SUM, 0);
    if (rank == 0) {
#endif

        double lr = new_like - curr_like;  //transition and prior ratios are zero
        lr *= heat;
        accept = (lr > 0 ? true : frand() < exp(lr));

        int pop = *(it->pops.begin());
        printLog(LOG_LOW, "%i\t%.1f\t%.1f\t%f\t%f\t%f\t%s\n",
                 index, num_coal[pop], num_nocoal[pop],
                 model->popsizes[pop], adjust[0], lr,
                 accept ? "accept" : "reject");

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(&accept, 1, MPI::BOOL, 0);
#endif

    if (!accept) {
        for (list<PopsizeConfigParam>::iterator it2=l.begin();
             it2 != l.end(); it2++) {
            for (set<int>::iterator it3=it2->pops.begin();
                 it3 != it2->pops.end(); it3++) {
                model->popsizes[*it3] -= adjust[it2 == it ? 0 : 1];
            }
        }
    } else {
        curr_like = new_like;
    }
    return curr_like;
}


double resample_single_popsize_old(ArgModel *model, const LocalTrees *trees,
                               double heat,
                               const list<PopsizeConfigParam>::iterator &it,
                               double curr_like, int index) {
    list<PopsizeConfigParam> &l = model->popsize_config.params;
    int other_dir;
    double adjust;
    bool accept;
    int other_index;
    int pop, other_pop;
    list<PopsizeConfigParam>::iterator other=it;
    other++;
    bool is_first = (it == l.begin());
    bool is_last = (other == l.end());
    if (is_first && is_last) return curr_like;

#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    if (rank == 0) {
#endif

        other_dir = (frand() < 0.5 ? 1 : -1);

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(&other_dir, 1, MPI::INT, 0);
#endif

    if (other_dir == 1 && other == l.end()) {
        other = l.begin();
        other_index = 0;
    } else if (other_dir == -1 && it == l.begin()) {
        other = l.end();
        other--;
        other_index = l.size()-1;
    } else {
        if (other_dir == -1) {
            other--;
            other--;
        }
        other_index = index + other_dir;
    }

#ifdef ARGWEAVER_MPI
    if (rank == 0) {
#endif

        adjust = 500.0 * frand();

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(&adjust, 1, MPI::DOUBLE, 0);
#endif

    set<int>::iterator it2 = other->pops.begin();
    other_pop = *it2;
    if (model->popsizes[other_pop] <= adjust) return curr_like;
    for ( ; it2 != other->pops.end(); it2++)
        model->popsizes[*it2] -= adjust;
    it2 = it->pops.begin();
    pop = *it2;
    for ( ; it2 != it->pops.end(); it2++)
        model->popsizes[*it2] += adjust;
    double new_like = calc_arg_prior_recomb_integrate(model, trees);

#ifdef ARGWEAVER_MPI
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : &new_like,
                 &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
    if (rank == 0) {
#endif

        double lr = new_like - curr_like;  //transition and prior ratios are zero
        lr *= heat;
        accept = (lr > 0 ? true : frand() < exp(lr));

        printLog(LOG_LOW, "(%i,%i)\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                 index, other_index,
                 model->popsizes[pop] - adjust,
                 model->popsizes[other_pop] + adjust,
                 adjust, model->popsizes[pop] - model->popsizes[other_pop],
                 model->popsizes[pop] - model->popsizes[other_pop] - 2.0 * adjust,
                 lr,
                 accept ? "accept" : "reject");

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(&accept, 1, MPI::BOOL, 0);
#endif

    if (!accept) {
        for (it2 = other->pops.begin(); it2 != other->pops.end(); it2++)
            model->popsizes[*it2] += adjust;
        for (it2 = it->pops.begin(); it2 != it->pops.end(); it2++)
            model->popsizes[*it2] -= adjust;
    } else {
        curr_like = new_like;
    }
    return curr_like;
}


void resample_popsizes(ArgModel *model, const LocalTrees *trees,
                       bool sample_popsize_recomb, double heat) {
    list<PopsizeConfigParam> &l = model->popsize_config.params;
    double curr_like = sample_popsize_recomb ? calc_arg_prior(model, trees) :
        calc_arg_prior_recomb_integrate(model, trees);
#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : &curr_like,
                 &curr_like, 1, MPI::DOUBLE, MPI_SUM, 0);
#endif
    //TEMPORARY: should make this an option if i like it this way
    /*    for (int rep=0; rep < model->popsize_config.numsample; rep++)
        curr_like = resample_popsize_scale(model, trees, sample_popsize_recomb,
        heat, curr_like);*/

    //    for (int rep=0; rep < model->popsize_config.numsample; rep++) {
    for (int rep=0; rep < model->popsize_config.numsample; rep++) {
        int idx=0;
        for (list<PopsizeConfigParam>::iterator it = l.begin();
             it != l.end(); it++) {
            curr_like =
                resample_single_popsize(model, trees, sample_popsize_recomb,
                                        heat, it, curr_like, idx++);
        }
        curr_like = resample_popsize_scale(model, trees, sample_popsize_recomb,
                                           heat, curr_like);
    }

}


void resample_popsizes_old(ArgModel *model, const LocalTrees *trees,
                       double heat) {

#ifdef ARGWEAVER_MPI
    printLog(0, "resample_popsizes %i\t%i\n", MPI::COMM_WORLD.Get_rank(),
             model->mc3.group_comm->Get_rank()); fflush(stdout);
    //    MPI::COMM_WORLD.Barrier();
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    printLog(0, "rank=%i\n", rank);
    if (rank == 0) {
#endif
        int num_accept=0, total=0;
        list<PopsizeConfigParam> &l = model->popsize_config.params;
        double *num_coal, *num_nocoal;
        vector<double> lrs(2 * model->ntimes),
            trans(2 * model->ntimes),
            prior(2 * model->ntimes),
            oldn(2 * model->ntimes),
            newn(2 * model->ntimes),
            praccept(2 * model->ntimes);
        vector<int> accepted(2 * model->ntimes);
        num_coal = (double*)malloc(2 * model->ntimes * sizeof(double));
        num_nocoal = (double*)malloc(2 * model->ntimes * sizeof(double));
        double curr_like = calc_arg_prior_recomb_integrate(model, trees,
                                                           //        double curr_like = calc_arg_prior(model, trees,
                                                           num_coal, num_nocoal);
#ifdef ARGWEAVER_MPI
        comm->Reduce(MPI_IN_PLACE, &curr_like, 1, MPI::DOUBLE, MPI_SUM, 0);
        comm->Reduce(MPI_IN_PLACE, num_coal, 2 * model->ntimes - 1, MPI::DOUBLE,
                     MPI_SUM, 0);
        comm->Reduce(MPI_IN_PLACE, num_nocoal, 2 * model->ntimes - 1, MPI::DOUBLE,
                     MPI_SUM, 0);
#endif
        for (int rep=0; rep < model->popsize_config.numsample; rep++) {
            for (list<PopsizeConfigParam>::iterator it = l.begin();
                          it != l.end(); it++) {
                if (it->sample == false) continue;
                int maxpop=-1;
                for (set<int>::iterator it2 = it->pops.begin();
                     it2 != it->pops.end(); it2++)
                    if ((*it2) > maxpop) maxpop = *it2;
                double old_popsize = model->popsizes[maxpop];
                double s = min(500.0, old_popsize / 2.0);
                s *= s;  //variance of gamma proposal from
                         // old_popsize to new_popsize

                double new_popsize = rand_gamma(old_popsize * old_popsize / s,
                                                s / old_popsize);
#ifdef ARGWEAVER_MPI
                comm->Bcast(&new_popsize, 1, MPI::DOUBLE, 0);
#endif
                double sp = min(500.0, new_popsize / 2.0);
                sp *= sp;  //variance of proposal from new_popsize
                           // to old_popsize
                double logn = log(old_popsize); //log N
                double lognp = log(new_popsize); //log N'
                double nsquare = old_popsize * old_popsize;
                double npsquare = new_popsize * new_popsize;
                double trans_ratio = (npsquare / sp - nsquare / s - 1.0) * logn
                    + (1.0 - nsquare / s + npsquare / sp) * lognp
                    - old_popsize * new_popsize / sp
                    + old_popsize * new_popsize / s
                    - npsquare / sp * log(sp)
                    + nsquare / s* log(s)
                    - lgamma(npsquare / sp)
                    + lgamma(nsquare / s);

                // using an uninformative gamma prior which slowly goes to zero
                // as you move out to infinity, still allows N to be at least a
                // million (maybe this is too big?)
                // has mean of 200000 (k*theta) and sd of 200000, and is pretty
                // flat from 0 to 400k or so
                double prior_ratio;
                double prior_theta = 200000;
                // double prior_k=1.0;  //this is the value but it is
                // commented-out since never used
                if (( ! model->popsize_config.neighbor_prior) ||
                    maxpop >= 2 * model->ntimes - 2)
                    prior_ratio = (old_popsize - new_popsize) / prior_theta;
                else {
                    double prev_popsize = model->popsizes[maxpop + 1];
                    double pneighbor=0.99999;
                    static double neighbor_sigma = 50.0;
                    static double neighbor_sigma22 = 2.0 * 50.0 * 50.0;
                    static double neighbor_scale = 1.0 / (neighbor_sigma *
                                                          sqrt(2.0 * 3.141593));
                    double newprior = (1.0 - pneighbor)*
                        (exp(-new_popsize / prior_theta) / prior_theta) +
                        pneighbor * neighbor_scale *
                        exp(-(new_popsize - prev_popsize) *
                            (new_popsize - prev_popsize) / neighbor_sigma22);
                    double oldprior = (1.0 - pneighbor) *
                        (exp(- old_popsize / prior_theta) / prior_theta) +
                        pneighbor * neighbor_scale *
                        exp(- (old_popsize - prev_popsize) *
                            (old_popsize - prev_popsize) / neighbor_sigma22);
                    prior_ratio = log(newprior / oldprior);
                }


                for (set<int>::iterator it2 = it->pops.begin();
                     it2 != it->pops.end(); it2++)
                    model->popsizes[*it2] = new_popsize;
                double new_like = calc_arg_prior_recomb_integrate(model, trees);
                //                double new_like = calc_arg_prior(model, trees);
#ifdef ARGWEAVER_MPI
                comm->Reduce(MPI_IN_PLACE, &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
#endif
                double lr = new_like - curr_like;
                double ln_accept = trans_ratio + prior_ratio + lr;
                ln_accept *= heat;
                double pr_accept = (ln_accept > 0 ? 1.0 : exp(ln_accept));
                bool accept = (ln_accept > 0 || frand() < pr_accept);
#ifdef ARGWEAVER_MPI
                comm->Bcast(&accept, 1, MPI::BOOL, 0);
#endif
                for (set<int>::iterator it2 = it->pops.begin();
                     it2 != it->pops.end(); it2++) {
                    lrs[*it2] = new_like - curr_like;
                    trans[*it2] = trans_ratio;
                    prior[*it2] = prior_ratio;
                    oldn[*it2] = old_popsize;
                    newn[*it2] = new_popsize;
                    accepted[*it2] = (accept == true);
                    praccept[*it2] = pr_accept;
                }
                if (accept) {
                    num_accept++;
                    curr_like = new_like;
                } else {
                    for (set<int>::iterator it2 = it->pops.begin();
                         it2 != it->pops.end(); it2++) {
                        model->popsizes[*it2] = old_popsize;
                    }
                }
                total++;
            }
        }
        printLog(LOG_LOW, "done resample_popsizes num_accept=%i/%i\n",
                 num_accept, total);
        for (int i=0; i < 2*model->ntimes-1; i++) {
            for (list<PopsizeConfigParam>::iterator it = l.begin();
                 it != l.end(); it++) {
                if (it->pops.find(i) != it->pops.end()) {
                    if (it->sample) {
                        printLog(LOG_LOW,
                                 "%i\t%.1f\t%.1f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                                 i, num_coal[i], num_nocoal[i], oldn[i], newn[i],
                                 lrs[i], trans[i], prior[i],
                                 accepted[i] == 1 ? "accept" : "reject");
                    } else {
                        printLog(LOG_LOW,
                                 "%i\t%.1f\t%.1f\t%f\tnot_sampled\n",
                                 i, num_coal[i], num_nocoal[i],
                                 model->popsizes[i]);
                    }
                }
            }
        }
        fflush(stdout);
        free(num_coal);
        free(num_nocoal);
#ifdef ARGWEAVER_MPI
    } else {
        list<PopsizeConfigParam> l = model->popsize_config.params;
        double *num_coal = (double*)malloc(2 * model->ntimes*sizeof(double));
        double *num_nocoal = (double*)malloc(2 * model->ntimes * sizeof(double));
        double curr_like = calc_arg_prior_recomb_integrate(model, trees,
                                                           //        double curr_like = calc_arg_prior(model, trees,
                                                           num_coal, num_nocoal);
        comm->Reduce(&curr_like, &curr_like, 1, MPI::DOUBLE, MPI_SUM, 0);
        comm->Reduce(num_coal, num_coal, 2 * model->ntimes - 1, MPI::DOUBLE,
                     MPI_SUM, 0);
        comm->Reduce(num_nocoal, num_nocoal, 2 * model->ntimes - 1, MPI::DOUBLE,
                     MPI_SUM, 0);

        for (int rep=0; rep < model->popsize_config.numsample; rep++) {
            for (list<PopsizeConfigParam>::iterator it=l.begin();
                 it != l.end(); it++) {
                if (it->sample == false) continue;
                double old_popsize = model->popsizes[*(it->pops.begin())];
                double new_popsize;
                comm->Bcast(&new_popsize, 1, MPI::DOUBLE, 0);
                for (set<int>::iterator it2 = it->pops.begin();
                     it2 != it->pops.end(); it2++)
                    model->popsizes[*it2] = new_popsize;
                double new_like = calc_arg_prior_recomb_integrate(model, trees);
                //                double new_like = calc_arg_prior(model, trees);
                comm->Reduce(&new_like, &new_like, 1, MPI::DOUBLE, MPI_SUM, 0);
                bool accept;
                comm->Bcast(&accept, 1, MPI::BOOL, 0);
                if (!accept) {
                    for (set<int>::iterator it2 = it->pops.begin();
                         it2 != it->pops.end(); it2++) {
                        model->popsizes[*it2] = old_popsize;
                    }
                }
            }
        }
        free(num_coal);
        free(num_nocoal);
    }
    printLog(0, "done resample popsizes\n");
#endif
}


#define DERIV_EPSILON 1e-6      /* for numerical computation of
                                   derivatives */
#define BOUNDARY_EPS2 1.0e-6     /* use this smaller value in 1d case */
#define ITMAX 200               /* maximum allowed number of
                                   iterations (opt_bfgs) */
#define EPS 1e-15               /* approx machine precision */
#define TOLX_HIGH (4*EPS)       /* convergence criterion for param */
#define ALPHA 1.0e-4            /* threshold for sufficient decrease
                                   in function value (lnsrch) */

#define RHO 0.5                 /* constant for geometric decrease of
                                   step size during backtracking in
                                   one dimensional linesearch */
#define INFTY 999999999


/* compute first and (optionally) second derivative at particular
   abscissa, using numerical approximations to derivatives if
   necessary.  Allows for bounds.  For use in one-dimensional
   optimizers.  Set deriv2 == NULL to skip second derivative.  Note:
   function value fx is only used in the numerical case (avoids a
   redundant function evaluation) */
void opt_derivs_1d(double *deriv, double *deriv2, double x, double fx, 
		   double lb, double ub, double (*f)(double, void*), void *data,
		   double (*compute_deriv)(double x, void *data, double lb, 
					   double ub),
		   double (*compute_deriv2)(double x, void *data, double lb, 
					    double ub),
		   double deriv_epsilon) {
    double fxeps=-1.0, fx2eps=-1.0;
    int at_ub = (ub - x < BOUNDARY_EPS2); /* at upper bound */
    
    if (compute_deriv == NULL) {
	if (at_ub) {           /* use backward method if at upper bound */
	    fxeps = f(x - deriv_epsilon, data);
	    *deriv = (fx - fxeps) / deriv_epsilon;
	}
	else { 
	    fxeps = f(x + deriv_epsilon, data);
	    *deriv = (fxeps - fx) / deriv_epsilon;
	}
    }
    else 
	*deriv = compute_deriv(x, data, lb, ub);
    
    if (deriv2 == NULL) return;
    
    if (compute_deriv2 == NULL) {     /* numerical 2nd deriv */
	if (at_ub) {                    /* at upper bound */
	    if (compute_deriv != NULL)    /* exact 1d available */
		*deriv2 = (*deriv - compute_deriv(x - deriv_epsilon, data, lb, ub)) / 
		    deriv_epsilon;
	    else {                    /* numerical 1st and second derivs */
		fx2eps = f(x - 2*deriv_epsilon, data);
		*deriv2 = (fx2eps + 2*fxeps - fx) / (deriv_epsilon * deriv_epsilon);
	    }
	}
	else {                       /* not at upper bound */
	    if (compute_deriv != NULL) /* exact 1d available */
		*deriv2 = (compute_deriv(x + deriv_epsilon, data, lb, ub) - *deriv) / 
		    deriv_epsilon;
	    else {                    /* numerical 1st and second derivs */
		fx2eps = f(x + 2*deriv_epsilon, data);
		*deriv2 = (fx2eps - 2*fxeps + fx) / (deriv_epsilon * deriv_epsilon);
	    }
	}
    }
    else                          /* exact 2nd deriv */
	*deriv2 = compute_deriv2(x, data, lb, ub);
}
    

/* line search for 1d case.  Sets *x, *fx, *final_lambda on exit */
void opt_lnsrch_1d(double direction, double xold, double fxold, double *x, 
		   double *fx, double deriv, double (*func)(double, void*), 
		   void *data, int *nevals, double *final_lambda, FILE *logf) {
    
    double lambda, slope, test, lambda_min;
    
    /* one-d line search */
    lambda = 1;
    slope = deriv * direction;
    /* Let phi(lambda) = f(xold + direction * lambda).  Then
       phi'(lambda) is the rate of change of the new function value as a
       function of lambda.  phi'(lambda) can be shown to equal to
       direction * f'(xold + direction * lambda).  Thus, slope is
       phi'(lambda) evaluated lambda = 0; in other words, it is the
       downward slope as a function of lambda at xold.  See Nocedal and
       Write, Numerical Optimization, chapter 3, for a reasonably clear
       discussion */
    
    test = fabs(*x)/max(fabs(xold), 1.0);
    lambda_min = TOLX_HIGH/test;
    
    for (;;) {
	(*x) = xold + lambda * direction;
	(*fx) = func(*x, data);
	(*nevals)++;
	
	if (lambda < lambda_min) {
	    (*x) = xold;
	    *final_lambda = lambda;
	    break;
	}
	else if (*fx <= fxold + ALPHA * lambda * slope) {
	    /* the "sufficient decrease" (Armijo) condition has been met--
	       essentially, the slope achieved by the update is less (i.e.,
	       a steeper negative slope) than a constant times the tangent
	       at xold.  This ensures that the update slope is not becoming
	       increasingly smaller, so that the algorithm converges to a
	       suboptimal value  */
	    *final_lambda = lambda;
	    break;
	}
	
	lambda *= RHO;              /* have to backtrack */
	/* Simple geometric decrease in lambda still allows for guaranteed
	   convergence.  See Nocedal and Write, pp 41-42 */
    }
    
    *final_lambda = lambda;
}



/* given two numbers, return number of shared significant figures */
int opt_sigfig(double val1, double val2) {
    double tmp;
    int sf, tv1, tv2;
    if (val1 == val2) return INFTY;
    if ((val1 < 0 && val2 > 0) || (val1 > 0 && val2 < 0)) return 0;
    val1 = fabs(val1); val2 = fabs(val2);
    tmp = pow(10, floor(log10(val1)));
    val1 /= tmp; val2 /= tmp;
    for (sf = 0; sf < 30; sf++) { /* never look at more than 30 digits */
	tv1 = (int)floor(val1); tv2 = (int)floor(val2);
	if (tv1 != tv2 || (val1 < 1e-30 && val2 < 1e-30)) break;
	/* avoid pathological roundoff cases */
	val1 = (val1 - tv1) * 10;   /* avoid overflow */
	val2 = (val2 - tv2) * 10;
    }    
    return sf;
}


/* one-dimensional Newton-Raphson optimizer, with line search and
   bounds.  Allows for analytical computation of both first and second
   derivatives via function pointers.  Set these to NULL to use
   numerical methods.  Abscissa (*x) should be initialized
   appropriately by calling code, and will contain optimized value on
   exit.  Parameter (*fx) will contain minimized value of function on
   exit.  Set sigfigs to desired number of stable significant figures
   for convergence.  This criterion applies both to x and to f(x).
   Function returns 0 on success, 1 if maximum number of iterations is
   reached */
bool opt_newton_1d(double (*f)(double, void*), double (*x), void *data, 
		   double *fx, int sigfigs, double lb, double ub, FILE *logf, 
		   double (*compute_deriv)(double x, void *data, double lb, 
					   double ub),
		   double (*compute_deriv2)(double x, void *data, double lb, 
					    double ub)) {
     
    double xold, fxold, d, d2, direction, lambda = -1;
    int its, nevals = 0;
    bool converged = false;
    struct timeval start_time, end_time;
    
    if (!(*x > lb && *x < ub && ub > lb)) {
	fprintf(stderr, "ERROR opt_newton_1d: x=%e, lb=%e, ub=%e\n", *x, lb, ub);
	assert(0);
    }
    
    if (logf != NULL) {
	gettimeofday(&start_time, NULL);
	fprintf(logf, "%15s %15s %15s %15s %15s\n", "f(x)", "x", "f'(x)", 
		"f''(x)", "lambda");
    }
    
    /* initial function evaluation */
    (*fx) = f(*x, data);
    nevals++;
    
    xold = (*x);                  /* invariant condition at loop start */
    fxold = (*fx);
    
    for (its = 0; !converged && its < ITMAX; its++) { 
	opt_derivs_1d(&d, &d2, *x, *fx, lb, ub, f, data, compute_deriv, 
		      compute_deriv2, DERIV_EPSILON);
	nevals += 2;                /* assume cost of each deriv approx
				       equals that of a functional evaluation */
	
	if (logf != NULL)               
	    fprintf(logf, "%15.6f %15.6f %15.6f %15.6f %15.6f\n", *fx, *x, d, d2, lambda);
	
	if (d2 < 1e-4)
	    d2 = 1;                    /* if second deriv is negative or very
					  close to zero, reduce to simple
					  gradient descent */
	
	direction = -d / d2;
	
	/* truncate for bounds, if necessary */
	if ((*x) + direction - lb < BOUNDARY_EPS2) 
	    direction = lb + BOUNDARY_EPS2 - (*x); 
	else if (ub - ((*x) + direction) < BOUNDARY_EPS2) 
	    direction = ub - BOUNDARY_EPS2 - (*x);
	
	/* line search; function eval occurs here */
	opt_lnsrch_1d(direction, xold, fxold, x, fx, d, f, data, &nevals, 
		      &lambda, logf);
	
	/* test for convergence */
	if (opt_sigfig(*x, xold) >= sigfigs && opt_sigfig(*fx, fxold) >= sigfigs) 
	    converged = true;
	
	fxold = (*fx);
	xold = (*x);
    }  
    
    if (logf != NULL) {
	fprintf(logf, "%15.6f %15.6f %15s %15s %15f\n", *fx, *x, "-", "-", lambda);
	gettimeofday(&end_time, NULL);
	fprintf(logf, "\nNumber of iterations: %d\nNumber of function evaluations: %d\nTotal time: %.4f sec.\n", 
		its, nevals, end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec)/1.0e6);
	
	if (!converged)
	    fprintf(logf, "WARNING: exceeded maximum number of iterations.\n");
    }
    
    return(!converged);
}
    
    /*double popsize_likelihood_func(double log_popsize, void *data0) {
    struct popsize_mle_data *data = (struct popsize_mle_data*)data0;
    int t = data->popsize_idx;
    double old_popsize1, old_popsize2, prob;
    ArgModel *model = data->model;
    double popsize = exp(log_popsize);
    old_popsize1 = model->popsizes[2*t];
    model->popsizes[2*t] = popsize;
    if (t > 0) {
	old_popsize2 = model->popsizes[2*t-1];
	model->popsizes[2*t-1] = popsize;
    }
    //    printf("calling calc_arg_prior_recomb_integrate %i %f\n", t, popsize);
    prob = calc_arg_prior_recomb_integrate(model, data->trees);
    if (t == 0) {
	model->popsizes[0] = old_popsize1;
    } else {
	model->popsizes[2*t-1] = old_popsize1;
	model->popsizes[2*t] = old_popsize2;
    }
    //    printf("popsize_likelihood_func\t%i\t%f\t%f\n", t, popsize, prob);
    return -prob;
    }*/

double mle_one_popsize_likelihood(double log_popsize, struct popsize_mle_data *data);

void numeric_deriv_check(double log_popsize, struct popsize_mle_data *data, double deriv0) {
    double delta=1e-5;
    double like1 = mle_one_popsize_likelihood(log_popsize, data);
    double like2 = mle_one_popsize_likelihood(log_popsize+delta, data);
    double deriv = (like2-like1)/delta;
    printf("check: %f\t%f\t%f\n", deriv, deriv0, fabs(deriv-deriv0));
}


void mle_one_popsize_like_and_dlike(double log_popsize, struct popsize_mle_data *data, 
				    double *likelihood, double *dlikelihood) {
    double popsize = exp(log_popsize);
    double t1 = data->t1;
    double t2 = data->t2;
    int numleaf = data->numleaf;
    int **coal_counts = data->coal_counts[data->popsize_idx];
    int **nocoal_counts = data->nocoal_counts[data->popsize_idx];
    double like=0.0, dlike=0.0;
    /*    if (log_popsize < 0) {
	*likelihood = -INFINITY;
	*dlikelihood = 100.0; //want to go up from here
	}*/
    for (int i=0; i < numleaf; i++) {
	for (int j=0; j < numleaf; j++) {
	    double rate = (t1*i + t2*j)/(2.0*popsize);
	    double erate=exp(-rate);
	    double erate1 = 1.0-erate;
            if (j==0) assert(coal_counts[i][j] == 0 && nocoal_counts[i][j]==0);
	    if (coal_counts[i][j] > 0) {
		like += coal_counts[i][j]*log(erate1);
		dlike -= coal_counts[i][j]/(erate1)*erate*rate;
	    }
	    if (nocoal_counts[i][j] > 0) {
		like -= nocoal_counts[i][j]*rate;
		dlike += nocoal_counts[i][j]*rate;
	    }
	}
    }
    if (0) {
	// these values always seem very close, not checking anymore
	numeric_deriv_check(log_popsize, data, dlike);
    }
    *likelihood = like;
    *dlikelihood = dlike;
}


double mle_one_popsize_likelihood(double log_popsize, struct popsize_mle_data *data) {
    double popsize = exp(log_popsize);
    double t1 = data->t1;
    double t2 = data->t2;
    int numleaf = data->numleaf;
    int **coal_counts = data->coal_counts[data->popsize_idx];
    int **nocoal_counts = data->nocoal_counts[data->popsize_idx];
    double pr=0.0;
    /*    if (log_popsize < 0.0)
	  return -INFINITY;*/
    for (int i=0; i < numleaf; i++) {
	for (int j=0; j < numleaf; j++) {
	    double rate = (t1*i + t2*j)/(2.0*popsize);
            if (j==0) assert(coal_counts[i][j] == 0 && nocoal_counts[i][j]==0);
	    if (coal_counts[i][j] > 0) pr += coal_counts[i][j]*log(1.0-exp(-rate));
	    if (nocoal_counts[i][j] > 0) pr -= nocoal_counts[i][j]*rate;
	}
    }
    return pr;
}


//return first derivative of likelihood
double mle_one_popsize_dlikelihood(double log_popsize, struct popsize_mle_data *data) {
    double popsize = exp(log_popsize);
    double t1 = data->t1;
    double t2 = data->t2;
    int numleaf = data->numleaf;
    int **coal_counts = data->coal_counts[data->popsize_idx];
    int **nocoal_counts = data->nocoal_counts[data->popsize_idx];
    double pr=0.0;
    for (int i=0; i < numleaf; i++) {
	for (int j=0; j < numleaf; j++) {
	    double rate = (t1*i + t2*j)/(2.0*popsize);
	    if (j==0) assert(coal_counts[i][j] == 0 && nocoal_counts[i][j]==0);
	    if (coal_counts[i][j] > 0) pr -= coal_counts[i][j]/(1.0-exp(-rate))*exp(-rate)*rate;
	    if (nocoal_counts[i][j] > 0) pr += nocoal_counts[i][j]*rate;
	    //note: above formulas should be divided by popsize, but then multiplied by popsize
	    // to take into account that we are maximizing x=log(popsize)
	    // log_like(no_coal) = count*t/2N = count*t/(2*exp(x))
            // dlike = -count*t/(2*exp(x))^2 (2 exp(x)) = -count*t/(2*exp(x)) = 
	}
    }
    return pr;
}


double mle_one_popsize_neg_likelihood(double log_popsize, void *data0) {
    struct popsize_mle_data *data = (struct popsize_mle_data*) data0;
    double like = mle_one_popsize_likelihood(log_popsize, data);
    return -like;
}

double mle_one_popsize(double init_popsize, void *data0) {
    double popsize, log_popsize = log(init_popsize);
    double likelihood;
    int sigfigs=4;
    struct popsize_mle_data *data = (struct popsize_mle_data*)data0;
    static double min_popsize = log(100);
    static double max_popsize = log(1e7);
    opt_newton_1d(mle_one_popsize_neg_likelihood, &log_popsize, data0, &likelihood, sigfigs, min_popsize, max_popsize, NULL, NULL, NULL);
    popsize = exp(log_popsize);
    printf("mle_popsize %i\t%f\t%f\t%i\t%i\n", data->popsize_idx, popsize, likelihood, data->coal_totals[data->popsize_idx], data->nocoal_totals[data->popsize_idx]);
    return popsize;
}


void popsize_sufficient_stats(struct popsize_mle_data *data, ArgModel *model, const LocalTrees *trees) {
    int end = trees->start_coord;
    LineageCounts lineages(model->ntimes);
    int numleaf = trees->get_num_leaves();
    //coal_counts[i][j][k] gives number of SPRs which coalesce in time i, with j
    // lineages in the tree interval before time i, and k lineages in the
    // interval after time i. If coalescence happens at the same time as 
    // the recombination (which is implied if coal_time==0), then j is always 0
    //nocoal_counts is the same, but for non-coalescing segments; so counted
    // for each segment from the recomb up until before the coal
    int arr_size = 2*(model->ntimes * numleaf * numleaf + model->ntimes);
    int *arr_alloc = new int[arr_size]();
    int ***coal_counts = new int**[model->ntimes];
    int ***nocoal_counts = new int**[model->ntimes];
    int *coal_totals = &arr_alloc[0];
    int *nocoal_totals = &arr_alloc[model->ntimes];
    int pos = model->ntimes * 2;
    int pseudocount=1;

#ifdef ARGWEAVER_MPI
    //Set pseudocount to zero for all but one MPI, since it will all get combined
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    if (rank > 0) pseudocount = 0;
#endif
    for (int i=0; i < model->ntimes; i++) {
	coal_counts[i] = new int*[numleaf];
	nocoal_counts[i] = new int*[numleaf];
	for (int j = 0; j < numleaf; j++) {
	    coal_counts[i][j] = &(arr_alloc[pos]);
	    pos += numleaf;
	    nocoal_counts[i][j] = &(arr_alloc[pos]);
	    pos += numleaf;
	}
	if (i==0) {
	    coal_counts[i][0][1] = pseudocount;
	    nocoal_counts[i][0][1] = pseudocount;
	} else {
	    coal_counts[i][1][1] = pseudocount;
	    nocoal_counts[i][1][1] = pseudocount;
	}
	coal_totals[i] += pseudocount;
	nocoal_totals[i] += pseudocount;
    }
    if (pos != arr_size) {
	printf("pos=%i arr_size=%i ntimes=%i numleaf=%i\n", pos, arr_size, model->ntimes, numleaf);
	assert(pos == arr_size);
    }
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
	end += it->blocklen;
	LocalTree *tree = it->tree;
	++it;
	if (end >= trees->end_coord) break;
	const Spr *spr = &it->spr;
	lineages.count(tree);
	int broken_age = tree->nodes[tree->nodes[spr->recomb_node].parent].age;
	int nlineage1=0;
	int nlineage2=lineages.nbranches[spr->recomb_time] - int(spr->recomb_time < broken_age);

	if (spr->recomb_time == spr->coal_time) {
	    coal_counts[spr->coal_time][0][nlineage2]++;
	    coal_totals[spr->coal_time]++;
	} else {
	    nocoal_counts[spr->recomb_time][0][nlineage2]++;
	    nocoal_totals[spr->recomb_time]++;
	}
	for (int i=spr->recomb_time + 1; i < spr->coal_time; i++) {
	    nlineage1 = nlineage2;
	    nlineage2 = lineages.nbranches[i] - int(i < broken_age);
	    nocoal_counts[i][nlineage1][nlineage2]++;
	    nocoal_totals[i]++;
	}
	if (spr->recomb_time != spr->coal_time) {
	    nlineage1 = nlineage2;
	    nlineage2 = lineages.nbranches[spr->coal_time] - int(spr->coal_time < broken_age);
	    coal_counts[spr->coal_time][nlineage1][nlineage2]++;
	    coal_totals[spr->coal_time]++;
	}
    }
    data->coal_counts = coal_counts;
    data->nocoal_counts = nocoal_counts;
    data->coal_totals = coal_totals;
    data->nocoal_totals = nocoal_totals;
    data->numleaf = numleaf;
    data->ntimes = model->ntimes;
    data->arr_alloc = arr_alloc;
    data->arr_size = arr_size;
}

void delete_popsize_data(struct popsize_mle_data *data) {
    int ntimes = data->ntimes;
    for (int i=0; i < ntimes; i++) {
	delete data->coal_counts[i];
	delete data->nocoal_counts[i];
    }
    delete data->coal_counts;
    delete data->nocoal_counts;
    delete data->arr_alloc;
}


void set_data_time(struct popsize_mle_data *data, int t, ArgModel *model) {
    data->popsize_idx = t;
    data->t2 = model->coal_time_steps[2*t];
    if (t == 0)
	data->t1 = 0;
    else data->t1 = model->coal_time_steps[2*t-1];
}


//use Hamiltonian MC to update popsize
void update_popsize_hmc(ArgModel *model, const LocalTrees *trees) {
    struct popsize_mle_data data;

    //compute all the coal_counts and nocoal_counts to be used for likelihood calculations
    popsize_sufficient_stats(&data, model, trees);

#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : data.arr_alloc, data.arr_alloc, data.arr_size, MPI::INT, MPI_SUM, 0);
    if (rank == 0) {
#endif

    static double sd=1;
    static double epsilon=0.1;
    double momentum[model->ntimes];
    static int numsteps=100;
    
    // first re-sample the momentums
    double current_K=0.0;
    for (int i=0; i < model->ntimes - 1; i++) {
	momentum[i] = rand_norm(0, sd);
	current_K += momentum[i]*momentum[i]/2.0;
    }

    double log_popsizes[model->ntimes];
    double orig_popsizes[model->ntimes];  // just for printing
    double current_likelihood=0;
    for (int i=0; i < model->ntimes-1; i++) {
	double like, dlike;
	log_popsizes[i] = log(model->popsizes[2*i]);
	orig_popsizes[i] = model->popsizes[2*i];
	set_data_time(&data, i, model);
	mle_one_popsize_like_and_dlike(log_popsizes[i], &data, &like, &dlike);
	current_likelihood -= like;
	momentum[i] += epsilon * dlike / 2.0;
    }

    
    for (int step=1; step <= numsteps; step++) {
	for (int i=0; i < model->ntimes-1; i++) {
	    log_popsizes[i] += epsilon * momentum[i];
	    set_data_time(&data, i, model);
	    if (step != numsteps)
		momentum[i] += epsilon * mle_one_popsize_dlikelihood(log_popsizes[i], &data);
	}
    }

    double proposed_likelihood = 0.0;
    double proposed_K = 0.0;
    for (int i=0; i < model->ntimes-1; i++) {
	double like, dlike;
	set_data_time(&data, i, model);
	mle_one_popsize_like_and_dlike(log_popsizes[i], &data, &like, &dlike);
	proposed_likelihood -= like;
	momentum[i] += epsilon * dlike / 2.0;
	proposed_K += momentum[i] * momentum[i] / 2.0;
    }
    
    double lr = current_likelihood - proposed_likelihood + current_K - proposed_K;
    if (lr > 0 || frand() < exp(lr)) {
	printf("accept HMC update %f (%f %f %f %f)\n", lr, current_likelihood, proposed_likelihood, current_K, proposed_K);
	for (int i=0; i < model->ntimes-1; i++) {
	    model->popsizes[2*i] = exp(log_popsizes[i]);
	    if (i != 0) model->popsizes[2*i-1] = model->popsizes[2*i];
	    printf("%i\t%f\t%f\t%f\taccept\t%i\t%i\n", i, orig_popsizes[i], model->popsizes[2*i], momentum[i], data.coal_totals[i], data.nocoal_totals[i]);
	}
	printf("\n");
    } else {
	printf("reject HMC update %f (%f %f %f %f)\n", lr, current_likelihood, proposed_likelihood, current_K, proposed_K);
	for (int i=0; i < model->ntimes-1; i++)
	    printf("%i\t%f\t%f\t%f\treject\t%i\t%i\n", i, orig_popsizes[i], exp(log_popsizes[i]), momentum[i], data.coal_totals[i], data.nocoal_totals[i]);
	printf("\n");
    }

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(model->popsizes, model->ntimes*2-1, MPI::DOUBLE, 0);
#endif

    delete_popsize_data(&data);
}




void mle_popsize(ArgModel *model, const LocalTrees *trees) {
    struct popsize_mle_data data;
    popsize_sufficient_stats(&data, model, trees);
	
    for (int i=0; i < model->ntimes-1; i++) {
	set_data_time(&data, i, model);
	model->popsizes[2*i] = mle_one_popsize(model->popsizes[2*i], &data);
	if (i > 0) model->popsizes[2*i-1] = model->popsizes[2*i];
    }
    delete_popsize_data(&data);
}



// reset popsize for each time interval to maximum likelihood value
/*void mle_popsize_2(ArgModel *model, const LocalTrees *trees)
{
    double maxlike, new_popsize;
    struct popsize_mle_data data;
    int sigfigs=4;
    static double min_popsize=log(100);
    static double max_popsize=log(1e7);
    data.model = model;
    data.trees = trees;
    
    for (int i=0; i < model->ntimes - 1; i++) {
	data.popsize_idx = i;
	new_popsize = log(model->popsizes[2*i]);
	opt_newton_1d(popsize_likelihood_func, &new_popsize, &data,
		      &maxlike, sigfigs, min_popsize, max_popsize, NULL, NULL, NULL);
	model->popsizes[2*i] = exp(new_popsize);
	if (i > 0) model->popsizes[2*i-1] = model->popsizes[2*i];
	printf("mle_popsize %i\t%f\t%f\n", i, exp(new_popsize), maxlike);
    }
}




double popsize_mle_get_prob(ArgModel *model, const LocalTrees *trees, int t, double popsize) {
    double old_popsize1, old_popsize2, prob;
    if (t == 0) {
	old_popsize1 = model->popsizes[0];
	model->popsizes[0] = popsize;
    } else {
	old_popsize1 = model->popsizes[2*t-1];
	old_popsize2 = model->popsizes[2*t];
	model->popsizes[2*t-1] = popsize;
	model->popsizes[2*t] = popsize;
    }
    prob = calc_arg_prior_recomb_integrate(model, trees);
    if (t == 0) {
	model->popsizes[0] = old_popsize1;
    } else {
	model->popsizes[2*t-1] = old_popsize1;
	model->popsizes[2*t] = old_popsize2;
    }
    return prob;
}


// reset popsize for each time interval to maximum likelihood value
void mle_popsize_old(ArgModel *model, const LocalTrees *trees)
{
    double minx, maxx, midx, minpr, maxpr, midpr;
    double tol=1;

    for (int i=0; i < model->ntimes - 1; i++) {
	minx = 1;
	maxx = 100000;
	midx = 30000;
	minpr = popsize_mle_get_prob(model, trees, i, minx);
	midpr = popsize_mle_get_prob(model, trees, i, midx);
	maxpr = popsize_mle_get_prob(model, trees, i, maxx);
	if (midpr < maxpr || midpr < minpr) {
	    fprintf(stderr, "Error: bad bounds for golden section search in mle_popsize\n");
	    assert(0);
	}
	while (maxx - minx > tol) {
	    double newx, newpr;
 	    static double p = (1.0 + sqrt(5))/2.0;  //golden ratio
	    bool new_is_greater = (maxx - midx > midx - minx);
	    if (new_is_greater) {
		newx = (midx - minx)/p + midx;
		assert(newx > midx && newx < maxx);
	    } else {
		newx = midx - (maxx - midx)/p;
		assert(newx > minx && newx < midx);
	    }
	    newpr = popsize_mle_get_prob(model, trees, i, newx);
	    if (new_is_greater) {
		if (newpr > midpr) {
		    minpr = midpr;
		    minx = midx;
		    midpr = newpr;
		    midx = newx;
		} else {
		    maxpr = newpr;
		    maxx = newx;
		}
	    } else {
		if (newpr > midpr) {
		    maxpr = midpr;
		    maxx = midx;
		    midpr = newpr;
		    midx = newx;
		} else {
		    minpr = newpr;
		    minx = newx;
		}
	    }
	}
	if (i == 0) {
	    model->popsizes[i] = midx;
	} else {
	    model->popsizes[2*i-1] = midx;
	    model->popsizes[2*i] = midx;
	}
    }
    }*/



void est_popsize_trees2(const ArgModel *model, const LocalTree *const *trees,
                        int ntrees, double *popsizes)
{
    assert(ntrees > 0);

    const int ntimes = model->ntimes;
    const int nleaves = trees[0]->get_num_leaves();
    LineageCounts lineages(ntimes);

    int total_ncoals[ntimes];
    int total_pairs[ntimes];
    int total_ncoals_pairs[ntimes];

    fill(total_ncoals, total_ncoals + ntimes, 0);
    fill(total_pairs, total_pairs + ntimes, 0);
    fill(total_ncoals_pairs, total_ncoals_pairs + ntimes, 0);

    printf("ntrees %d\n", ntrees);

    // count lineages
    for (int i=0; i<ntrees; i++) {
        lineages.count(trees[i]);

        for (int j=0; j<ntimes-1; j++) {
            int start = (j == 0) ? nleaves : lineages.nbranches[j-1];
            int end = lineages.nbranches[j];
            int ncoals = start - end;
            int pairs = start * (start-1) / 2;

            total_ncoals[j] += ncoals;
            total_pairs[j] += pairs;
            total_ncoals_pairs[j] += ncoals * pairs;
        }
    }

    for (int j=0; j<ntimes-1; j++) {
        if (total_ncoals[j] == 0)
            popsizes[j] = 0.0;
        else
            popsizes[j] = .5 * model->time_steps[j] *
                (total_ncoals_pairs[j] + total_pairs[j] - total_ncoals[j]) /
                double(total_ncoals[j]);
        printf("> %d %d %d\n", total_ncoals_pairs[j], total_pairs[j],
               total_ncoals[j]);
        printf("popsize %f\n", popsizes[j]);
    }
}



void est_popsize_trees(const ArgModel *model, const LocalTree *const *trees,
                       int ntrees, double *popsizes)
{
    assert(ntrees > 0);

    const int ntimes = model->ntimes;
    const int nleaves = trees[0]->get_num_leaves();
    LineageCounts lineages(ntimes);

    int total_ncoals[ntimes];
    int total_pairs[ntimes];

    fill(total_ncoals, total_ncoals + ntimes, 0);
    fill(total_pairs, total_pairs + ntimes, 0);

    printf("ntrees %d\n", ntrees);

    // count lineages
    for (int i=0; i<ntrees; i++) {
        lineages.count(trees[i]);

        for (int j=0; j<ntimes-1; j++) {
            int start = (j == 0) ? nleaves : lineages.nbranches[j-1];
            int end = lineages.nbranches[j];
            int ncoals = start - end;
            int pairs = start * (start-1) / 2;

            total_ncoals[j] += ncoals;
            total_pairs[j] += pairs;
        }
    }

    for (int j=0; j<ntimes-1; j++) {
        if (total_ncoals[j] == 0)
            popsizes[j] = 0.0;
        else
            popsizes[j] = .5 * model->time_steps[j] * total_pairs[j] /
                double(total_ncoals[j]);
        printf("> %d %d\n", total_pairs[j], total_ncoals[j]);
        printf("popsize %f\n", popsizes[j]);
    }
}


void est_popsize_trees(const ArgModel *model, const LocalTrees *trees,
                       int step, double *popsizes)
{
    vector<LocalTree*> indep_trees;

    int end = trees->start_coord;
    int pos = end;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
        {
            int start = end;
            end += it->blocklen;

            while (start <= pos && pos < end) {
                // record tree
                indep_trees.push_back(it->tree);
                pos += step;
            }
        }

    est_popsize_trees(model, &indep_trees[0], indep_trees.size(), popsizes);
}


//=============================================================================
// C-interface

extern "C" {

// estimate population sizes
void arghmm_est_popsizes_trees(LocalTrees *trees, double *times, int ntimes,
                               int step, double *popsizes)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, 0.0);
    est_popsize_trees(&model, trees, step, popsizes);
}


} // extern C

} // namespace argweaver

