/*
 * zdensity.h
 *
 *  Created on: 16.11.2009
 *      Author: daniel
 */

#ifndef ZDENSITY_H_
#define ZDENSITY_H_

#include <iwls.h>
#include <types.h>
#include <linApproxDens.h>
#include <rcppExport.h>

// 19/5/2011: generalise "binaryLogisticCorrection" to "higherOrderCorrection"
class NegLogUnnormZDens {
public:

    // call the function object:
    // z is the argument,
    double
    operator()(double z);

    // constructor
    NegLogUnnormZDens(const ModelPar &mod,
                      const GlmModelData& modelData,
                      // return the approximate *conditional* density f(y | z, mod) by operator()?
                      // otherwise return the approximate unnormalized *joint* density f(y, z | mod).
                      bool conditional,
                      bool verbose,
                      bool higherOrderCorrection,
                      PosInt nIter=40) :
                          mod(mod),
                          modelData(modelData),
                          iwlsObject(mod,
                                     modelData,
                                     conditional, // fixed z if conditional density of y given z is wished
                                     EPS, // take EPS as the convergence epsilon
                                     verbose), // and debug if verbose is wished.
                          verbose(verbose),
                          higherOrderCorrection(higherOrderCorrection),
                          nIter(nIter)
    {
    }

private:
    // save the model reference and fp info, so we can write a nice warning message if the
    // IWLS fails for some z
    const ModelPar& mod;
    const GlmModelData& modelData;

    // the IWLS object
    Iwls iwlsObject;

    // be verbose?
    const bool verbose;

    // make the Laplace correction for canonical models?
    const bool higherOrderCorrection;

    // number of IWLS iterations
    PosInt nIter;
};


// compute the log marginal likelihood of a specific model
double
glmGetLogMargLik(const ModelPar& modPar,
                 const GlmModelData& modelData,
                 const Computation& computation);

// compute the log marginal likelihood of a specific model
// another interface with other rhos
double
glmGetLogMargLik(const GlmModelData& modelData,
                 const Rcpp::IntegerVector& whichLinear,
                 const Rcpp::IntegerVector& whichOptimize,
                 const Rcpp::NumericVector& rhos,
                 const Computation& computation);

// ZdensApprox: class encapsulating the approximation of the
// marginal posterior density of z = log(g).
class ZdensApprox {
public:

    // ctr
    ZdensApprox(const ModelPar& modPar,
                const GlmModelData& modelData,
                const Computation& computation);

    // get the log density value at a specific z
    double
    logDens(double z) const;

    // get one sample from the approximate distribution
    double
    sample() const;

    // getter for the mode of the density approximation
    double
    getzMode() const
    {
        return zMode;
    }

    // getter for logMargLik
    double
    getLogMargLik() const
    {
        return logMargLik;
    }

    // getter for the linear approximation
    Rcpp::List
    getLinApproxDens() const;

private:

    // the mode, is computed inside ctr.
    double zMode;

    // the ILA estimate (is a byproduct)
    double logMargLik;

    // is this the null model without any z?
    const bool isNullModel;

    // the internal approximation object
    LinApproxDens linApproxDens;
};


#endif /* ZDENSITY_H_ */
