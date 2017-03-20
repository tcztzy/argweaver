// headers c++
#include "common.h"


namespace argweaver {

/* make a draw from a gamma distribution with parameters 'a' and
 * 'b'. Be sure to call srandom externally.  If a == 1, exp_draw is
 * called.  If a > 1, Best's (1978) rejection algorithm is used, and
 * if a < 1, rejection sampling from the Weibull distribution is
 * performed, both as described in "Non-Uniform Random Variate
 * Generation" by Luc Devroye, available online at
 * http://cgm.cs.mcgill.ca/~luc/rnbookindex.html */
// NOTE: "a" and "b" refer to parameters "k" and "theta"
// as described on Wikipedia, so that the expected value is a*b
// and variance is a*b*b
 // NOTE: This is buggy for small a! when c=1/a is big,
 // then pow(Z, c) often hits infinity
double rand_gamma(double a, double b) {
  double retval = -1;

  if (a <= 0) {
    fprintf(stderr, "ERROR gamma_draw got a=%f\n", a);
    exit(1);
  }
  if (a == 1) return rand_exp(b);

  else if (a > 1) {
    while (retval == -1) {
      double U, V, W, X, Y, Z;
      double d = a - 1, c = 3 * a - 0.75;

      U = frand();          /* uniform on [0, 1]; used for draw */
      V = frand();          /* also uniform on [0, 1]; used for
                                       rejection/acceptance */
      W = U * (1 - U);
      Y = sqrt(c / W) * (U - 0.5);
      X = d + Y;
      /* Y is a scaled version of a random variate from a t distribution
         with 2 dof; X is a random deviate from a shifted distribution
         whose ratio with a gamma(a, 1) can be bounded by a constant */

      if (X < 0) continue;      /* truncate because of nonnegativity
                                   of gamma */

      Z = 64 * W * W * W * V * V;
      if (log(Z) <= 2 * (d * log(X / d) - Y))
        retval = X;
      /* there's some algebra behind this, but underneath it's just
         the standard rejection sampling idea of accepting with
         probability p(x)/(M*q(x)), where p(x) is the density of the
         desired distrib, q(x) is the density of the distrib from which
         you have sampled, and M is a constant such that p(x) / q(x) <=
         M for all x */
    }
  }
  else {
      double c = 1/a;
      double ld = log(1-a) + a/(1-a)*log(a);
      while (retval == -1) {
          double lE = log(rand_exp(1));
          double lZ = log(rand_exp(1));
          double lX = c*lZ;
          if (logadd(lZ, lE) >= logadd(ld, lX))
              retval = exp(lX);
      }
  }
  /*  else {                        // use Weibull
    double c = 1/a, d = pow(a, a/(1-a)) * (1-a);
    while (retval == -1) {
      double E, Z, X;
      E = rand_exp(1);
      Z = rand_exp(1);
      X = pow(Z, c);            // X is Weibull(a)
      if (Z + E >= d + X)       // note: wrong in book, correct
                               // formula in errata
        retval = X;
    }
    } */

  /* so far we only have a draw from gamma(a, 1); multiply by b to
 *      obtain a draw from gamma(a, b) */
  return retval * b;
}


double rand_beta(double a, double b, bool do_test) {
    double x = rand_gamma(a, 1);
    double y = rand_gamma(b, 1);
    double rv = x/(x+y);
    if (do_test) {
        int numtest=100000;
        double vals[numtest];
        double mean=0,var=0;
        for (int l=0; l < numtest; l++) {
            vals[l] = rand_beta(a, b, false);
            mean += vals[l];
        }
        mean /= (double)numtest;
        for (int l=0; l < numtest; l++) {
            var += (vals[l] - mean)*(vals[l]-mean);
        }
        var /= (double)(numtest-1);
        double expmean = a/(a+b);
        double expvar = a*b/((a+b)*(a+b)*(a+b+1));
        printf("mean %e %f %f var=%e %f %f\n",
               (mean-expmean)/expmean, mean, expmean,
               (var-expvar)/expvar, var, expvar);
    }
    return rv;
}
}
