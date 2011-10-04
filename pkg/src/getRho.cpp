/*
 * getRho.cpp
 *
 *  Created on: 28.07.2011
 *      Author: daniel
 */

#include <getRho.h>

#include <rcppExport.h>
#include <R_ext/Applic.h>



using namespace Rcpp;

struct TargetData {
    const NumericVector& lambdas;
    double degree;

    TargetData(const NumericVector& lambdas,
               double degree) :
                   lambdas(lambdas),
                   degree(degree)
    {
    }
};

double
target(double rho, TargetData* data)
{
    double tmp = sum(data->lambdas / (data->lambdas + 1.0 / rho));
    return tmp - data->degree;
}

double
getRho(const NumericVector& lambdas,
       double degree)
{
    // generate arguments for R_zeroin:
    TargetData data(lambdas,
                    degree);
    double tol = pow(DOUBLE_EPS,
                     0.25);
    int maxIter = 10000;

    // do the "uniroot" call
    double ret = R_zeroin(DOUBLE_XMIN,
                          DOUBLE_XMAX,
                          (double (*)(double, void *)) target,
                          (void *)(&data),
                          &tol,
                          &maxIter);

    // check convergence
    if(maxIter < 0)
    {
        Rf_warning("root-finder in getRho did not converge in %d iterations!",
                   maxIter);
    }

    // if OK, return
    return ret;
}



