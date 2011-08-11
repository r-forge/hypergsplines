/*
 * random.h
 *
 *  Created on: 17.03.2011
 *      Author: daniel
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <types.h>

// draw a single random normal vector from N(mean, (precisionCholeskyFactor * t(precisionCholeskyFactor))^(-1))
AVector
drawNormalVector(const AVector& mean,
                 const AMatrix& precisionCholeskyFactor);

// draw a single uniform random variable
double
unif();

// get random int x with lower <= x < upper
template<class INT>
INT
discreteUniform(const INT& lower, const INT& upper)
{
    if (lower >= upper)
    {
        Rf_error("\nlower = %d >= %d = upper in discreteUniform call\n", lower,
                 upper);
    }

    double u = unif();

    INT size = upper - lower;
    INT ret = lower;

    while (u > 1.0 / size * (ret - lower + 1))
    {
        ret++;
    }

    return ret;
}


#endif /* RANDOM_H_ */
