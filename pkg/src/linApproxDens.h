/*
 * linApproxDens.h
 *
 *  Created on: 17.03.2011
 *      Author: daniel
 */

#ifndef LINAPPROXDENS_H_
#define LINAPPROXDENS_H_

#include <types.h>

// LinApproxDens:
// This class gets vectors x and log f(x) where
// f is the unnormalised density. It then produces
// a linear spline interpolation density from that.
class LinApproxDens {
public:

    // ctr
    LinApproxDens(const AVector& args,
                  const AVector& logDens);

    // default ctr: only for convenience, object not to be used in reality!
    LinApproxDens(){}

    // inverse cdf sampling from the linear density approximation
    double
    sample() const;

    // return the density approximation at one argument value
    double
    dens(double arg) const;

    // getters for the internal stuff
    DoubleDeque
    getArgs() const
    {
        return ord_args;
    }

    DoubleDeque
    getDens() const
    {
        return ord_dens;
    }

private:

    // sorted argument values:
    DoubleDeque ord_args;

    // normalised density values and cumulative distribution function:
    DoubleDeque ord_dens;
    DoubleVector cdf;

    // number of points, i.e. the dimension of ord_args, ord_dens and cdf:
    PosInt nPoints;
};






#endif /* LINAPPROXDENS_H_ */
