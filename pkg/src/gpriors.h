/*
 * gpriors.h
 *
 *  Created on: 15.03.2010
 *      Author: daniel
 *
 * Classes for priors on the shrinkage paramater g.
 *
 * 30/03/2012: add the hyper-g/n prior
 *
 */

#ifndef GPRIORS_H_
#define GPRIORS_H_

#include <rcppExport.h>

#include <functionWraps.h>

// ***************************************************************************************************//

// Virtual base class for all g priors.
struct GPrior
{
    // Log prior density
    virtual double
    logDens(double g) const = 0;

    // we need a virtual destructor here,
    // cf. Accelerated C++ pp. 242 ff.
    virtual ~GPrior(){}
};

// ***************************************************************************************************//

// Inverse gamma g prior
class InvGammaGPrior : public GPrior
{
public:
    // ctr
    InvGammaGPrior(double a, double b) :
        a(a),
        b(b)
        {
        }

    // Log prior density
    double
    logDens(double g) const
    {
        return - (a + 1.0) * log(g) - b / g + a * log(b) - Rf_lgammafn(a);
    }

private:
    // the parameters for the inverse gamma density
    const double a;
    const double b;
};

// ***************************************************************************************************//

// Hyper-g prior
class HypergPrior : public GPrior
{
public:
    // ctr
    HypergPrior(double a) :
        a(a)
        {
        }

    // Log prior density
    double
    logDens(double g) const
    {
        return log(a - 2.0) - M_LN2 - (a / 2.0) * log1p(g);
    }

private:
    // the hyperparameter
    const double a;
};

// ***************************************************************************************************//

// Hyper-g/n prior
class HypergnPrior : public GPrior
{
public:
    // ctr
    HypergnPrior(double a,
                 int nObs) :
        a(a),
        n(static_cast<double>(nObs))
        {
        }

    // Log prior density
    double
    logDens(double g) const
    {
        return log(a - 2.0) - M_LN2 - log(n) - (a / 2.0) * log1p(g/n);
    }

private:
    // the hyperparameter
    const double a;

    // the sample size
    const double n;
};

// ***************************************************************************************************//

// Custom g-prior
class CustomGPrior : public GPrior
{
public:
    // ctr
    CustomGPrior(SEXP R_function) :
        wrappedRfunction(R_function)
        {
        }

    // Log prior density
    double
    logDens(double g) const
    {
        return wrappedRfunction(g);
    }

private:
    // the wrapped R function
    const RFunction wrappedRfunction;
};

// ***************************************************************************************************//


#endif /* GPRIORS_H_ */
