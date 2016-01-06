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

void numeric_deriv_check(double log_popsize, struct popsize_data *data, double deriv0) {
    double delta=1e-5;
    double like1 = one_popsize_likelihood(data->popsize_idx, log_popsize, data);
    double like2 = one_popsize_likelihood(data->popsize_idx, log_popsize+delta, data);
    double deriv = (like2-like1)/delta;
    printf("check: %f\t%f\t%f\n", deriv, deriv0, fabs(deriv-deriv0));
}


void one_popsize_like_and_dlike(int t, double log_popsize, struct popsize_data *data,
				double *likelihood, double *dlikelihood, double *dlikelihood2) {
    if (data->popsize_idx != t) {
	set_data_time(data, t);
    }

    double popsize = exp(log_popsize);
    double t1 = data->t1;
    double t2 = data->t2;
    int numleaf = data->numleaf;
    double **coal_counts = data->coal_counts[data->popsize_idx];
    double **nocoal_counts = data->nocoal_counts[data->popsize_idx];
    double like=0.0, dlike=0.0, dlike2=0.0;
    bool do_like = (likelihood != NULL);
    bool do_dlike = (dlikelihood != NULL);
    bool do_dlike2 = (dlikelihood2 != NULL);
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
		if (do_like) like += coal_counts[i][j]*log(erate1);
		// NOTE: we are optimizing log_popsize so need to take this into consideration for derivative computations
		if (do_dlike) dlike -= coal_counts[i][j]/(erate1)*erate*rate;
		if (do_dlike2) dlike2 -= coal_counts[i][j]*
				   ((erate*rate/erate1)*(erate*rate/erate1)
				    + erate*rate*(rate - 1.0)/erate1);
	    }
	    if (nocoal_counts[i][j] > 0) {
		if (do_like) like -= nocoal_counts[i][j]*rate;
		if (do_dlike) dlike += nocoal_counts[i][j]*rate;
		if (do_dlike2) dlike2 -= nocoal_counts[i][j]*rate;
	    }
	}
    }
    if (0) {
	// these values always seem very close, not checking anymore
	numeric_deriv_check(log_popsize, data, dlike);
    }
    if (do_like) *likelihood = like;
    if (do_dlike) *dlikelihood = dlike;
    if (do_dlike2) *dlikelihood2 = dlike2;
}


double one_popsize_likelihood(int t, double log_popsize, struct popsize_data *data) {
    double like;
    one_popsize_like_and_dlike(t, log_popsize, data, &like, NULL);
    return like;
}


//log_popsize should be a vector of length ntimes-1 with one popsize per whole time interval
double popsize_likelihood(double *log_popsize, struct popsize_data *data) {
    double like=0.0;
    for (int i = 0 ; i < data->model->ntimes - 1; i++)
	like += one_popsize_likelihood(i, log_popsize[i], data);
    return like;
}


//return first derivative of likelihood
double one_popsize_dlikelihood(int t, double log_popsize, struct popsize_data *data) {
     double dlike;
     one_popsize_like_and_dlike(t, log_popsize, data, NULL, &dlike);
     return dlike;
}



//wrapper to send to newton1d program; assume set_data_time has already been called
double one_popsize_neg_likelihood(double log_popsize, void *data0) {
    struct popsize_data *data = (struct popsize_data*) data0;
    double like=0;
    for (int t = data->min_t; t <= data->max_t; t++)
	like += one_popsize_likelihood(t, log_popsize, data);
    return -like;
}


//gets MLE for time interval start_t through end_t inclusive- one parameter for
// whole combined interval
double mle_one_popsize(int start_t, int end_t, double init_popsize, void *data0) {
    double popsize, log_popsize = log(init_popsize);
    double likelihood;
    int sigfigs=4;
    struct popsize_data *data = (struct popsize_data*)data0;
    static double min_popsize = log(100);
    static double max_popsize = log(1e7);
    data->min_t = start_t;
    data->max_t = end_t;
    opt_newton_1d(one_popsize_neg_likelihood, &log_popsize, data0, &likelihood, sigfigs, min_popsize, max_popsize, NULL, NULL, NULL);
    popsize = exp(log_popsize);
    double dlike2=0.0;
    for (int t=start_t; t <= end_t; t++)  {
	double tmp;
	one_popsize_like_and_dlike(t, log_popsize, data, NULL, NULL, &tmp);
	dlike2 += tmp;
    }
    double sd = sqrt(-1.0/dlike2);
    double popsize_min = exp(log_popsize - 2*sd);
    double popsize_max = exp(log_popsize + 2*sd);
    for (int t=start_t; t <= end_t; t++) {
	printLog(LOG_LOW, "mle_popsize %i\t%f\t%f\t%.1f\t%.1f\t%f\t%.1f\t%.1f\n", t, popsize, likelihood, data->coal_totals[t], data->nocoal_totals[t], sqrt(-1.0/dlike2), popsize_min, popsize_max);
    }
   return popsize;
}


    /*void popsize_sufficient_stats_recomb_integrate(struct popsize_data *data, ArgModel *model, const LocalTrees *trees) {
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
    double *arr_alloc = new double[arr_size]();
    double ***coal_counts = new double**[model->ntimes];
    double ***nocoal_counts = new double**[model->ntimes];
    double *coal_totals = &arr_alloc[0];
    double *nocoal_totals = &arr_alloc[model->ntimes];
    int pos = model->ntimes * 2;
    int pseudocount=1;

#ifdef ARGWEAVER_MPI
    //Set pseudocount to zero for all but one MPI, since it will all get combined
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    if (rank > 0) pseudocount = 0;
#endif
    for (int i=0; i < model->ntimes; i++) {
	coal_counts[i] = new double*[numleaf];
	nocoal_counts[i] = new double*[numleaf];
	for (int j = 0; j < numleaf; j++) {
	    coal_counts[i][j] = &(arr_alloc[pos]);
	    pos += numleaf;
	    nocoal_counts[i][j] = &(arr_alloc[pos]);
	    pos += numleaf;
	}
	if (pseudocount > 0) {
	    double pr_nocoal;
	    if (i==0) {
		pr_nocoal = exp(-model->coal_time_steps[0]/20000.0);
		coal_counts[i][0][1] = (1.0 - pr_nocoal) * pseudocount;
		nocoal_counts[i][0][1] = pr_nocoal * pseudocount;
	    } else {
		pr_nocoal = exp(-(model->coal_time_steps[2*i-1] + model->coal_time_steps[2*i])/20000.0);
		coal_counts[i][1][1] = (1.0 - pr_nocoal) * pseudocount;
		nocoal_counts[i][1][1] = pr_nocoal * pseudocount;
	    }
	    coal_totals[i] += (1.0 - pr_nocoal) * pseudocount;
	    nocoal_totals[i] += pr_nocoal * pseudocount;
	}
    }
    if (pos != arr_size) {
	printf("pos=%i arr_size=%i ntimes=%i numleaf=%i\n", pos, arr_size, model->ntimes, numleaf);
	assert(pos == arr_size);
    }
    int rho_idx = 0;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
	end += it->blocklen;
	LocalTree *tree = it->tree;

	int blocklen = it->blocklen;
	double treelen = get_treelen(tree, model->times, model->ntimes, false);
	double recomb_rate = max(model->get_local_rho(trees->start_coord, &rho_idx)*treelen,
				 model->rho);

	//for single site, probability of no recomb
	double pr_no_recomb = exp(-recomb_rate);
	double pr_recomb = 1.0 - pr_no_recomb;
	const int root_age = tree->nodes[tree->root].age;
	lineages.count(tree);
	if (end >= trees->end_coord)
	    blocklen++;

	if (blocklen > 1) {
	    double recomb_sum = 0.0;
	    for (int i=0; i < tree->nnodes; i++) {
		const LocalNode *node = &(tree->nodes[i]);
		if (node->parent != -1) {
		    int maxage = tree->nodes[node->parent].age;
		    assert(node->age <= maxage);
		    for (int age=node->age; age <= maxage; age++) {
			Spr spr(i, age, node->parent, maxage);
			double val = exp(calc_spr_prob(model, tree, spr, lineages,
						       treelen, NULL, NULL, 0, true));

	}

	++it;
	if (end >= trees->end_coord) break;
	const Spr *spr = &it->spr;
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
    data->arr_alloc = arr_alloc;
    data->arr_size = arr_size;
    data->model = model;
    data->popsize_idx = -1;
    data->t1 = -1;
    data->t2 = -1;
    }*/


void popsize_sufficient_stats(struct popsize_data *data, ArgModel *model, const LocalTrees *trees, bool add) {
    int end = trees->start_coord;
    LineageCounts lineages(model->ntimes);
    int numleaf = trees->get_num_leaves();
    //coal_counts[i][j][k] gives number of SPRs which coalesce in time i, with j
    // lineages in the tree interval before time i, and k lineages in the
    // interval after time i. If coalescence happens at the same time as
    // the recombination (which is implied if coal_time==0), then j is always 0
    //nocoal_counts is the same, but for non-coalescing segments; so counted
    // for each segment from the recomb up until before the coal
    double ***coal_counts;
    double ***nocoal_counts;
    double *coal_totals;
    double *nocoal_totals;

    if (!add) {
	int arr_size = 2*(model->ntimes * numleaf * numleaf + model->ntimes);
	double *arr_alloc = new double[arr_size]();
	double pseudocount;
	int pos=model->ntimes * 2;
	coal_counts = new double**[model->ntimes];
	nocoal_counts = new double**[model->ntimes];
	coal_totals = &arr_alloc[0];
	nocoal_totals = &arr_alloc[model->ntimes];
	pseudocount=model->popsize_config.pseudocount;
	//	printf("pseudocount=%f\n", pseudocount);
#ifdef ARGWEAVER_MPI
	//Set pseudocount to zero for all but one MPI, since it will all get combined
	MPI::Intracomm *comm = model->mc3.group_comm;
	int rank = comm->Get_rank();
	if (rank > 0) pseudocount = 0;
#endif
	for (int i=0; i < model->ntimes; i++) {
	    coal_counts[i] = new double*[numleaf];
	    nocoal_counts[i] = new double*[numleaf];
	    for (int j = 0; j < numleaf; j++) {
		coal_counts[i][j] = &(arr_alloc[pos]);
		pos += numleaf;
		nocoal_counts[i][j] = &(arr_alloc[pos]);
		pos += numleaf;
	    }
	    if (pseudocount > 0) {
		double pr_nocoal;
		if (i==0) {
		    pr_nocoal = exp(-model->coal_time_steps[0]/20000.0);
		    coal_counts[i][0][1] = (1.0 - pr_nocoal) * pseudocount;
		    nocoal_counts[i][0][1] = pr_nocoal * pseudocount;
		} else {
		    pr_nocoal = exp(-(model->coal_time_steps[2*i-1] + model->coal_time_steps[2*i])/20000.0);
		    coal_counts[i][1][1] = (1.0 - pr_nocoal) * pseudocount;
		    nocoal_counts[i][1][1] = pr_nocoal * pseudocount;
		}
		coal_totals[i] += (1.0 - pr_nocoal) * pseudocount;
		nocoal_totals[i] += pr_nocoal * pseudocount;
	    }
	}
	if (pos != arr_size) {
	    printf("pos=%i arr_size=%i ntimes=%i numleaf=%i\n", pos, arr_size, model->ntimes, numleaf);
	    assert(pos == arr_size);
	}
	data->arr_alloc = arr_alloc;
	data->arr_size = arr_size;
	data->coal_counts = coal_counts;
	data->nocoal_counts = nocoal_counts;
	data->coal_totals = coal_totals;
	data->nocoal_totals = nocoal_totals;
	data->numleaf = numleaf;
	data->model = model;
	data->popsize_idx = -1;
	data->t1 = -1;
	data->t2 = -1;
    } else {  // add counts to already initialized structure
	coal_counts = data->coal_counts;
	nocoal_counts = data->nocoal_counts;
	coal_totals = data->coal_totals;
	nocoal_totals = data->nocoal_totals;
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
}

void delete_popsize_data(struct popsize_data *data) {
    int ntimes = data->model->ntimes;
    for (int i=0; i < ntimes; i++) {
	delete [] data->coal_counts[i];
	delete [] data->nocoal_counts[i];
    }
    delete [] data->coal_counts;
    delete [] data->nocoal_counts;
    delete [] data->arr_alloc;
}


void set_data_time(struct popsize_data *data, int t) {
    data->popsize_idx = t;
    data->t2 = data->model->coal_time_steps[2*t];
    if (t == 0)
	data->t1 = 0;
    else data->t1 = data->model->coal_time_steps[2*t-1];
}

double dotProduct(double *x, int len) {
    double val=0.0;
    for (int i=0; i < len; i++)
	val += x[i]*x[i];
    return val;
}

double kineticEnergy(double *x, int len) {
    return dotProduct(x, len)/2.0;
}

    // perform L leapfrog steps
void leapFrogL(double *theta, double *r, double epsilon, int L,
	       double *thetaPrime, double *rPrime, struct popsize_data *data) {
    assert(L >= 1);
    int len = data->model->ntimes - 1;
    for (int i=0; i < len; i++) {
	rPrime[i] = r[i] + one_popsize_dlikelihood(i, theta[i], data)*epsilon/2.0;
	thetaPrime[i] = theta[i] + epsilon*rPrime[i];
    }

    for (int l=1; l < L; l++) {
	for (int i=0; i < len; i++) {
	    rPrime[i] += one_popsize_dlikelihood(i, thetaPrime[i], data)*epsilon;
	    thetaPrime[i] += epsilon*rPrime[i];
	}
    }

    for (int i=0; i < len; i++)
	rPrime[i] += one_popsize_dlikelihood(i, thetaPrime[i], data)*epsilon/2.0;
}

// Taken from Algorithm 6, No U-Turn Sampler with Dual Averaging
// www.stat.columbia.edu/~gelman/research/published/nuts.pdf
// "The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo"
// Matthew D Hoffman, Andrew Gelman, Journal of Machine Learning Resaerch 15 (2014) 1351-1381
//return variables are thetaMinus, rMinus, thetaPlus, rPlus, thetaPrime, nprime, sprime
/*void buildTree(double *theta, double *r, double u, double v, int j, double epsilon,
	       double *theta0, double *r0,
	       double *thetaMinus, double *rMinus, double *thetaPlus, double *rPlus,
	       double *thetaPrime, int *nprime, int *sprime, double *alpha, int *nalpha,
	       struct popsize_data *data) {
    static double deltaMax = 1000.0;
    int len = data->model->ntimes - 1;

    if (j == 0) {
	leapFrogL(theta, r, v*epsilon, 1, thetaPrime, rprime);
	double val = popsize_likelihood(thetaPrime, data) - kineticEnergy(rprime, len);
	if (u <= exp(val)) {
	    *nprime = 1;
	} else *nprime = 0;
	if (u < 0 || log(u) < deltaMax + val) {
	    *sprime = 1;
	} else *sprime = 0;
	if (thetaMinus != NULL) memcpy(thetaMinus, thetaPrime, len);
	if (rMinus != NULL) memcpy(rMinus, rPrime, len);
	if (thetaPlus != NULL) memcpy(thetaPlus, thetaPrime, len);
	if (rPlus != NULL) memcpy(rPlus, rPrime, len);
	val -= popsize_likelihood(theta0, data);
	val += kineticEnergy(r0);
	if (val > 0) {
	    *nprime
    }

    }


//use Hamiltonian MC to update popsize
//no u-turn sampler for adaptively choosing stepsize (and also dual averaging to choose epsilon)
void hmc_update_nuts(struct popsize_data *data, ArgModel *model) {
    static double sd=1;
    static double epsilon=model->popsize_config.epsilon;
    double momentum[model->ntimes];
    static int numsteps=100;

    // first re-sample the momentums
    double current_K=0.0;
    for (int i=0; i < model->ntimes - 1; i++) {
	double momentum_scale = (i == 0 ? 1.0 : log((model->coal_time_steps[2*i-1] + model->coal_time_steps[2*i])/model->coal_time_steps[0]) + 1.0);
	momentum[i] = rand_norm(0, sd*sqrt(momentum_scale));
	current_K += momentum[i]*momentum[i]/2.0;
    }

    double log_popsizes[model->ntimes];
    double orig_popsizes[model->ntimes];  // just for printing
    double current_likelihood=0;
    for (int i=0; i < model->ntimes-1; i++) {
	double like, dlike;
	log_popsizes[i] = log(model->popsizes[2*i]);
	orig_popsizes[i] = model->popsizes[2*i];
	one_popsize_like_and_dlike(i, log_popsizes[i], &data, &like, &dlike);
	current_likelihood -= like;
	momentum[i] += epsilon * dlike / 2.0;
    }


    for (int step=1; step <= numsteps; step++) {
	for (int i=0; i < model->ntimes-1; i++) {
	    log_popsizes[i] += epsilon * momentum[i];
	    if (step != numsteps)
		momentum[i] += epsilon * one_popsize_dlikelihood(i, log_popsizes[i], &data);
	}
    }

    double proposed_likelihood = 0.0;
    double proposed_K = 0.0;
    for (int i=0; i < model->ntimes-1; i++) {
	double like, dlike;
	one_popsize_like_and_dlike(i, log_popsizes[i], &data, &like, &dlike);
	proposed_likelihood -= like;
	momentum[i] += epsilon * dlike / 2.0;
	proposed_K += momentum[i] * momentum[i] / 2.0;
    }

    double lr = current_likelihood - proposed_likelihood + current_K - proposed_K;
    if (lr > 0 || frand() < exp(lr)) {
	printLog(LOG_LOW, "accept HMC update %f (%f %f %f %f)\n", lr, current_likelihood, proposed_likelihood, current_K, proposed_K);
	for (int i=0; i < model->ntimes-1; i++) {
	    model->popsizes[2*i] = exp(log_popsizes[i]);
	    if (i != 0) model->popsizes[2*i-1] = model->popsizes[2*i];
	    printLog(LOG_LOW, "%i\t%f\t%f\t%f\taccept\t%.1f\t%.1f\n", i, orig_popsizes[i], model->popsizes[2*i], momentum[i], data.coal_totals[i], data.nocoal_totals[i]);
	}
	printLog(LOG_LOW, "\n");
    } else {
	printLog(LOG_LOW, "reject HMC update %f (%f %f %f %f)\n", lr, current_likelihood, proposed_likelihood, current_K, proposed_K);
	for (int i=0; i < model->ntimes-1; i++)
	    printLog(LOG_LOW, "%i\t%f\t%f\t%f\treject\t%.1f\t%.1f\n", i, orig_popsizes[i], exp(log_popsizes[i]), momentum[i], data.coal_totals[i], data.nocoal_totals[i]);
	printLog(LOG_LOW, "\n");
    }

}
*/


//use Hamiltonian MC to update popsize
void hmc_update(struct popsize_data *data) {
    ArgModel *model = data->model;
    static double sd=1;
    static double epsilon=model->popsize_config.epsilon;
    double momentum[model->ntimes];
    static int numsteps=1000;

    // first re-sample the momentums
    double current_K=0.0;
    for (int i=0; i < model->ntimes - 1; i++) {
	// TODO: should momentum_scale just be 1? should test
	//	double momentum_scale = (i == 0 ? 1.0 : log((model->coal_time_steps[2*i-1] + model->coal_time_steps[2*i])/model->coal_time_steps[0]) + 1.0);
	double momentum_scale = 1.0;
	momentum[i] = rand_norm(0, sd*sqrt(momentum_scale));
	current_K += momentum[i]*momentum[i]/2.0;
    }

    double log_popsizes[model->ntimes];
    double orig_popsizes[model->ntimes];  // just for printing
    double current_likelihood=0;
    for (int i=0; i < model->ntimes-1; i++) {
	double like, dlike;
	log_popsizes[i] = log(model->popsizes[2*i]);
	orig_popsizes[i] = model->popsizes[2*i];
	one_popsize_like_and_dlike(i, log_popsizes[i], data, &like, &dlike);
	current_likelihood -= like;
	momentum[i] += epsilon * dlike / 2.0;
    }

    for (int step=1; step <= numsteps; step++) {
	for (int i=0; i < model->ntimes-1; i++) {
	    log_popsizes[i] += epsilon * momentum[i];
	    if (step != numsteps)
		momentum[i] += epsilon * one_popsize_dlikelihood(i, log_popsizes[i], data);
	}
    }

    double proposed_likelihood = 0.0;
    double proposed_K = 0.0;
    for (int i=0; i < model->ntimes-1; i++) {
	double like, dlike;
	one_popsize_like_and_dlike(i, log_popsizes[i], data, &like, &dlike);
	proposed_likelihood -= like;
	momentum[i] += epsilon * dlike / 2.0;
	proposed_K += momentum[i] * momentum[i] / 2.0;
    }

    double lr = current_likelihood - proposed_likelihood + current_K - proposed_K;
    if (lr > 0 || frand() < exp(lr)) {
	printLog(LOG_LOW, "accept HMC update %f (%f %f %f %f)\n", lr, current_likelihood, proposed_likelihood, current_K, proposed_K);
	for (int i=0; i < model->ntimes-1; i++) {
	    model->popsizes[2*i] = exp(log_popsizes[i]);
	    if (i != 0) model->popsizes[2*i-1] = model->popsizes[2*i];
	    printLog(LOG_LOW, "%i\t%f\t%f\t%f\taccept\t%.1f\t%.1f\n", i, orig_popsizes[i], model->popsizes[2*i], momentum[i], data->coal_totals[i], data->nocoal_totals[i]);
	}
	printLog(LOG_LOW, "\n");
    } else {
	printLog(LOG_LOW, "reject HMC update %f (%f %f %f %f)\n", lr, current_likelihood, proposed_likelihood, current_K, proposed_K);
	for (int i=0; i < model->ntimes-1; i++)
	    printLog(LOG_LOW, "%i\t%f\t%f\t%f\treject\t%.1f\t%.1f\n", i, orig_popsizes[i], exp(log_popsizes[i]), momentum[i], data->coal_totals[i], data->nocoal_totals[i]);
	printLog(LOG_LOW, "\n");
    }

}


void no_update_popsize(ArgModel *model, const LocalTrees *trees) {
    struct popsize_data data;

    //compute all the coal_counts and nocoal_counts to be used for likelihood calculations
    popsize_sufficient_stats(&data, model, trees);

#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : data.arr_alloc, data.arr_alloc, data.arr_size, MPI::DOUBLE, MPI_SUM, 0);
    if (rank == 0) {
#endif
	for (int i=0; i< model->ntimes-1; i++)
	    printLog(LOG_LOW, "%i\t%f\t%.f\t%.1f\n", i, model->popsizes[2*i], data.coal_totals[i], data.nocoal_totals[i]);

#ifdef ARGWEAVER_MPI
    }
#endif
    delete_popsize_data(&data);
}

//use Hamiltonian MC to update popsize
void update_popsize_hmc(ArgModel *model, const LocalTrees *trees) {
    struct popsize_data data;

    //compute all the coal_counts and nocoal_counts to be used for likelihood calculations
    popsize_sufficient_stats(&data, model, trees);

#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : data.arr_alloc, data.arr_alloc, data.arr_size, MPI::DOUBLE, MPI_SUM, 0);
    if (rank == 0) {
#endif
	hmc_update(&data);

#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(model->popsizes, model->ntimes*2-1, MPI::DOUBLE, 0);
#endif

    delete_popsize_data(&data);
}



void mle_popsize(ArgModel *model, const struct popsize_data *data, double min_total) {
    int start_time = 0;
    double curr_total = 0.0;
    for (int i=0; i < model->ntimes-1; i++) {
	curr_total += data->coal_totals[i] + data->nocoal_totals[i];
	if (curr_total < min_total && i < model->ntimes - 2) continue;
	double popsize = mle_one_popsize(start_time, i, model->popsizes[2*i], (void*)data);
	for (int j = start_time; j <= i; j++) {
	    model->popsizes[2*j] = popsize;
	    if (j > 0) model->popsizes[2*j-1] = popsize;
	}
	start_time = i+1;
	curr_total = 0.0;
    }
}


void mle_popsize(ArgModel *model, const LocalTrees *trees, double min_total) {
    struct popsize_data data;
    popsize_sufficient_stats(&data, model, trees);
#ifdef ARGWEAVER_MPI
    MPI::Intracomm *comm = model->mc3.group_comm;
    int rank = comm->Get_rank();
    comm->Reduce(rank == 0 ? MPI_IN_PLACE : data.arr_alloc, data.arr_alloc, data.arr_size, MPI::DOUBLE, MPI_SUM, 0);
    if (rank == 0) {
#endif
	mle_popsize(model, &data, min_total);
#ifdef ARGWEAVER_MPI
    }
    comm->Bcast(model->popsizes, model->ntimes*2-1, MPI::DOUBLE, 0);
#endif
    delete_popsize_data(&data);
}



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

