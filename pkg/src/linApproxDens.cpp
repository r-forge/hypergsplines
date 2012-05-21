/*
 * linApproxDens.cpp
 *
 *  Created on: 18.03.2011
 *      Author: daniel
 */

#include <linApproxDens.h>

#include <rcppExport.h>
// #include <cassert>
#include <random.h>

using namespace Rcpp;

// for testing from R
SEXP
cpp_linApproxDens(SEXP args,
                  SEXP logDens,
                  SEXP grid,
                  SEXP nSamples)
{
    NumericVector x = args;
    NumericVector ly = logDens;
    NumericVector g = grid;
    PosInt n = as<PosInt>(nSamples);

    LinApproxDens obj(x, ly);

    NumericVector samples(n);
    for(PosInt i = 0; i < n; ++i)
    {
        samples[i] = obj.sample();
    }

    NumericVector vals(g.size());
    for(int i = 0; i < g.size(); ++i)
    {
        vals[i] = obj.dens(g[i]);
    }

    return List::create(_["samples"] = samples,
                        _["vals"] = vals);

}


// ctr
LinApproxDens::LinApproxDens(const AVector& args,
                             const AVector& logDens)
{
    // assert(args.n_elem == logDens.n_elem);

    // get the (ordered) finite values
    arma::uvec order = arma::sort_index(args, 0);

    for(PosInt i = 0; i < args.n_elem; ++i)
    {
        PosInt index = order[i];

        // only take finite values (this also ensures that the
        // density values are positive, otherwise the log density is infinite)
        if(R_finite(logDens[index]))
        {
            ord_args.push_back(args[index]);
            ord_dens.push_back(exp(logDens[index]));
        }
    }

    // decide on the range of the spline function
    double argMin = ord_args.front();
    double argMax = ord_args.back();

    double constant = 0.2 * (argMax - argMin);

    ord_args.push_front(argMin - constant);
    ord_args.push_back(argMax + constant);

    // add zeroes to both ends of the density values
    // to smooth out the approximation to zero.
    ord_dens.push_front(0.0);
    ord_dens.push_back(0.0);

    nPoints = ord_args.size();

    // calculate the cdf points
    cdf.resize(nPoints);

    cdf[0] = 0.0;
    for(PosInt i = 1; i < nPoints; ++i)
    {
        cdf[i] = cdf[i - 1] + (ord_dens[i] + ord_dens[i - 1]) / 2.0 * (ord_args[i] - ord_args[i - 1]);
    }

    // normalize cdf and density values:
    double normConst = cdf[nPoints - 1];

    for(PosInt i = 0; i < nPoints; ++i)
    {
        cdf[i] /= normConst;
        ord_dens[i] /= normConst;
    }

    // check the finiteness of the approximation on a grid
    const PosInt gridLength = 400;
    const double increment = (ord_args[nPoints - 1] - ord_args[0]) / (gridLength - 1);
    double testArg = ord_args[0];
    for(PosInt i = 0; i < gridLength; ++i, testArg += increment)
    {
        if(! R_finite(dens(testArg)))
        {
            Rf_error("LinApproxDens: not finite approximation detected");
        }
    }
}

// inverse cdf sampling from the linear density approximation
double
LinApproxDens::sample() const
{
    // get a uniform random variate
    double p = unif();

    // check the boundary cases
    if(p == 0.0)
        return ord_args[0];
    if(p == 1.0)
        return ord_args[nPoints - 1];

    // so now we are inside.

    // find the right interval
    PosInt i = 0;
    while(p > cdf[i])
    {
        ++i;
    }
    // so now cdf[i - 1] < p <= cdf[i]

    double slope = (ord_dens[i] - ord_dens[i - 1]) / (ord_args[i] - ord_args[i - 1]);
    double intercept = ord_dens[i - 1];
    double constant = cdf[i - 1] - p;

    double solution = (sqrt(intercept * intercept - 2.0 * slope * constant) - intercept) / slope;

    return solution + ord_args[i - 1];
}

// return the density approximation at one argument value
double
LinApproxDens::dens(double arg) const
{
    int i = 0;
    int j = nPoints - 1;

    // handle out-of-domain points
    if (arg < ord_args[i])
        return ord_dens[i];
    if (arg > ord_args[j])
        return ord_dens[j];

    // find the correct interval by bisection:
    while (i < j - 1)
    { /* ord_args[i] <= arg <= ord_args[j] */
        int ij = (i + j) / 2; /* i+1 <= ij <= j-1 */

        if (arg < ord_args[ij])
            j = ij;
        else
            i = ij;
        /* still i < j */
    }
    /* probably have i == j-1 */

    /* interpolation */
    if (arg == ord_args[j])
        return ord_dens[j];
    if (arg == ord_args[i])
        return ord_dens[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

    return ord_dens[i] + (ord_dens[j] - ord_dens[i]) * ((arg - ord_args[i]) / (ord_args[j] - ord_args[i]));
}

