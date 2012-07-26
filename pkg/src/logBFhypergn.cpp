/*
 * logMargLikHypergn.cpp
 *
 *  Created on: 28.03.2012
 *      Author: daniel
 *
 *  Uses the C code from the R-package BAS
 */


#include <rcppExport.h>

#include <logBFhypergn.h>

using namespace Rcpp;


void
posroot(double a,
        double b,
        double c,
        double *root,
        double *status)
{ /* this computes the real roots of a cubic polynomial; in the end, if
 status==1, root stores the nonegative root; if status is not one,
 then status is the total number of nonegative roots and root is
 useless
 */

    int i;
    double Q, R, disc, Q3, A, B, aux, x[3];

    * root = 0.;
    * status = 0.;

    Q = (R_pow_di(a,
                  2) - 3. * b) / 9.;
    R = (2 * R_pow_di(a,
                      3) - 9 * a * b + 27. * c) / 54.;
    Q3 = R_pow_di(Q,
                  3);

    disc = R_pow_di(R,
                    2) - Q3;

    if (disc >= 0.)
    {
        if (R >= 0)
            A = - cbrt(R + sqrt(disc));
        else
            A = - cbrt(R - sqrt(disc));

        if (A == 0.)
            B = 0.;
        else
            B = Q / (A);
        * root = (A + B) - a / 3.;
        if (* root >= 0)
            * status = 1.;
    }
    else
    {
        A = acos(R / sqrt(Q3));
        aux = 2. * sqrt(Q);
        x[0] = - aux * cos(A / 3.);
        x[1] = - aux * cos((A + 4. * asin(1.)) / 3.);
        x[2] = - aux * cos((A - 4. * asin(1.)) / 3.);
        aux = a / 3.;
        for (i = 0; i < 3; i++)
            x[i] = x[i] - aux;
        for (i = 0; i < 3; i++)
        {
            if (x[i] >= 0.)
            {
                * status = * status + 1.;
                * root = x[i];
            }
        }
    }
}



double
lik_null_HG(double g,
            double R2,
            int n,
            int k,
            double alpha,
            int gpower)
{
    /* this computes log(likelihood  x prior), where the likelihood is marginal
     on the intercept, regression coefficients and the variance
     */
    double aux;

    aux = ((double) n - 1. - (double) k) * log(1. + g) - ((double) n - 1.)
            * log(1. + (1. - R2) * g) + 2. * (double) gpower * log(g) - alpha
            * log(1. + g / (double) n);
    aux = aux / 2.;
    aux = aux - log((double) n) + log(alpha / 2. - 1.);

    return (aux);
}


double
info_null_HG(double g,
             double R2,
             int n,
             int k,
             double alpha)
{/* This computs the second derivative of LogLik(tau) which is
 equal to
 (1/2)*(- alpha ng/(n+g)^2-(n-1) eg/(1+eg)^2+(n-1-p)g/(1+g)^2)
 where g=e^{\tau}
 */

    double aux;

    aux = ((double) n - 1. - (double) k) * g / R_pow_di(1. + g,
                                                        2);
    aux = aux - ((double) n - 1.) * (1. - R2) * g
            / R_pow_di(1. + (1. - R2) * g,
                       2);
    aux = aux - alpha * (double) n * g / R_pow_di((double) n + g,
                                                  2);
    aux = aux / 2.;
    return (aux);
}


// 15/6/12: use warnings instead of prints
// 18/6/12: if status != 1, return NaN and not NA
double
LogBF_Hg_null(double R2,
              int n,
              int d,
              double alpha,
              int gpower)
{
    /* this computes a Laplace approximation to the log of the Bayes factor
     with respect to the null model (intercept only), log(m_k)-log(m_0)

     R2 = 1-SSE/SST; the coefficient of determination
     e = 1 - R2;
     n  = sample size;
     k  = number of covariates of the current model (exclusing the intercept)

     The prior under consideration is Hyper-g with prior (1+g/n)^{-alpha/2}

     The cubic equation for loglik'(g)=0 is
     -e (alpha - 2 gpower + k) g^3
     - {[e (k - 2 gpower) - (1 - e)] n + (1 - e) + k + (1 + e) (alpha - 2 gpower)} g^2
     + [(1 - e) (n - 1) n - k n - alpha + 2 gpower (1 + (1 + e) n)] g + 2 gpower n
     */

    /* this version: April 11 2005 */

    double status, root, logmarg;
    double a, b, c, e, aux;
    int k;

    k = d - 1;
    e = 1. - R2;
    aux = - e * (alpha - (double) gpower * 2. + (double) k);
    a
            = - ((e * ((double) k - 2. * (double) gpower) - (1. - e))
                    * (double) n + (1. - e) + (double) k + (1. + e) * (alpha
                    - 2. * (double) gpower));
    b = ((1. - e) * ((double) n - 1.) * (double) n - (double) k * (double) n
            - alpha + 2. * (double) gpower * (1. + (1. + e) * (double) n));
    c = (double) n * 2. * (double) gpower;

    a = a / aux;
    b = b / aux;
    c = c / aux;
    posroot(a,
            b,
            c,
            & root,
            & status);
    if (status != 1.)
    {
        if (status == 0.)
            Rf_warning("\n No positive roots\n");
        else
            Rf_warning("\n More than one positive root\n");

        return R_NaN;
    }
    else
    {
        if (k == 0)
        {
            logmarg = 0.0;
        }
        else
        {
            logmarg = lik_null_HG(root,
                                  R2,
                                  n,
                                  k,
                                  alpha,
                                  gpower) + (log(4. * asin(1.))
                    - log(- info_null_HG(root,
                                         R2,
                                         n,
                                         k,
                                         alpha))) / 2.;
        }
        return (logmarg);
    }
}


// compute the log BF against the null model, either with the Appell
// function or with the Laplace approximation
double
logBFhypergn(int n,
             int p,
             double R2)
{
    // the return value will be put in here:
    double ret = 0.0;

    // if n is small, try the Appell function via the R interface
    if(n < 100)
    {
        Environment ns = Environment::namespace_env("hypergsplines");
        Function fn = ns["logBFhypergn"];
        ret = as<double>(fn(Named("n", n),
                            Named("p", p),
                            Named("R2", R2)));

        // if it is OK, return that
        if(! R_IsNA(ret))
        {
            return ret;
        }
    }

    // else or if the Appell function failed, use the Laplace approximation
    ret = LogBF_Hg_null(R2, n, p + 1, 4.0, 1);
    return ret;
}


// compute the Laplace approximation of the log marginal likelihood
// under the hyper-g/n prior.
SEXP
cpp_logMargLikHypergn(SEXP R_n, // number of observations
                      SEXP R_p, // number of covariates (excluding the intercept)
                      SEXP R_SSE0, // SSE in the null model
                      SEXP R_SSEm) // SSE in the current model
{
    // convert inputs
    int n = as<int>(R_n);
    int p = as<int>(R_p);
    double SSE0 = as<double>(R_SSE0);
    double SSEm = as<double>(R_SSEm);

    // compute R2
    double R2 = 1.0 - SSEm / SSE0;

    // compute log BF against the null model
    double ret = LogBF_Hg_null(R2, n, p + 1, 4.0, 1);

    // add log marginal likelihood of the null model
    double nm1h = (static_cast<double>(n) - 1.0) / 2.0;
    ret += - M_LN2 - nm1h * log(M_PI * SSE0) + lgamma(nm1h) - 0.5 * log(static_cast<double>(n));

    // so the result is the log marginal likelihood of the current model
    return wrap(ret);
}
