/*
 * random.cpp
 *
 *  Created on: 17.03.2011
 *      Author: daniel
 */

#include <random.h>

#include <rcppExport.h>
#include <linalgInterface.h>

// draw a single random normal vector from N(mean, (precisionCholeskyFactor * t(precisionCholeskyFactor))^(-1))
AVector
drawNormalVector(const AVector& mean,
                 const AMatrix& precisionCholeskyFactor)
{
    GetRNGstate();

    // get vector from N(0, I)
    AVector w = Rcpp::rnorm(mean.n_rows, // as many normal variates as required by the dimension.
                            0,
                            1);

    PutRNGstate();

    // then solve L' * ret = w, and overwrite w with the result:
    trs(false,
        true,
        precisionCholeskyFactor,
        w);

    // return the shifted vector
    return (w + mean);
}


// draw a single uniform random variable
double
unif()
{
    GetRNGstate();

    double ret = unif_rand();

    PutRNGstate();

    return ret;
}

