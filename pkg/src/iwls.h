/*
 * zdensity.h
 *
 *  Created on: 16.11.2009
 *      Author: daniel
 */

#ifndef IWLS_H_
#define IWLS_H_

#include <types.h>
#include <dataStructure.h>


class Iwls {

public:
    // constructor: constructs the Iwls object for given model and data and
    // start linear predictor
    Iwls(const ModelPar &mod,
         const GlmModelData& modelData,
         bool useFixedZ,
         double epsilon,
         bool verbose);

    // do the Iwls algorithm for a given covariance factor g and current linear predictor
    // linPred,
    // until convergence or until the maximum number of iterations is reached.
    // so also only one iwls step can be performed with this function.
    // returns the number of iterations.
    PosInt
    startWithLastLinPred(PosInt maxIter,
                         double g);

    // do the Iwls algorithm for a given covariance factor g and new start linear predictor
    // linPredStart.
    PosInt
    startWithNewLinPred(PosInt maxIter,
                        double g,
                        const AVector& linPredStart);

    // do the Iwls algorithm for a given covariance factor g and new start coefficients vector
    // coefsStart.
    PosInt
    startWithNewCoefs(PosInt maxIter,
                      double g,
                      const AVector& coefsStart);


    // compute the log of the (unnormalized)
    // posterior density for a given parameter consisting of the coefficients vector and z
    double
    computeLogUnPosteriorDens(const Parameter& sample) const;

    // getter for results
    IwlsResults
    getResults() const
    {
        return results;
    }


    // this can be public:

    // is this the null model?
    const bool isNullModel;

    // use a fixed z?
    const bool useFixedZ;

    // number of observations
    const PosInt nObs;

    // get the design matrix (by const reference):
    // (This shall avoid that we have to make (the constant) "design" public,
    // but can still access it.)
    // AMatrix getDesign()
    const AMatrix& getDesign() const
    {
        return design;
    }

    // getter for nCoefs
    PosInt getnCoefs() const
    {
        return nCoefs;
    }

private:

    // needed in function calls (mean, variance and other glm functions):
    const GlmModelData& modelData;

    // the convergence epsilon
    const double epsilon;

    // status messages??
    const bool verbose;


    // further intermediate results which are needed again internally:

    // This is diag(dispersions)^(-1/2), use diagmat to interpret it as diagonal matrix!
    const AVector invSqrtDispersions;

    // full design matrix for this model
    AMatrix design;

    // dimensions:
    int dimLinear; // number of linear coefficients
    int dimSpline; // number of spline terms
    int dimZ; // number of spline coefficients (dimension of Z)

    // dimension including the intercept ("ncol(design)" in R syntax)
    // this is 1 + dimLinear + dimZ
    PosInt nCoefs;

    // container for the results of the iwls computation:
    IwlsResults results;

    // unscaled prior precision matrix R^-1 = blockDiag(0, fixed part, diagonal spline coef part)
    // (without g^-1 which is later multiplied to the fixed part only!)
    // Saved in lower triangular storage.
    AMatrix unscaledPriorPrec;

    // the log of the determinant of the unscaled prior precision (obviously without
    // considering the intercept), so this is
    // log(det(blockDiag(fixed part, spline part))) ==
    // log(det(fixed part) * det(spline part)) ==
    // log(det(fixed part)) + log(det(spline part))
    double logUnscaledPriorPrecDeterminant;

    // lower-triangular Cholesky factor of the linear part of the prior precision matrix:
    AMatrix linPartCholFactor;

    // and here basically the spline coefficients equivalent, as the diagonal vector of that matrix.
    AVector splinePartPrecSqrtEntries;
};


#endif /* IWLS_H_ */
