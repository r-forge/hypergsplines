/*
 * getLogMargLikEst.cpp
 *
 *  Created on: 07.04.2011
 *      Author: daniel
 */

#include <getLogMargLikEst.h>

#include <rcppExport.h>

// hopefully fast kernel function
static inline void
addKron(arma::mat22& target,
        const arma::subview_col<double>& a,
        const arma::subview_col<double>& b)
{
    target.at(0, 0) += a[0] * b[0];
    target.at(0, 1) += a[0] * b[1];
    target.at(1, 0) += a[1] * b[0];
    target.at(1, 1) += a[1] * b[1];
}

// get one Omega
static arma::mat22
getOmega(const PosInt s, const AMatrix& hDiffHat)
{
    arma::mat22 ret;
    ret.zeros();

    for(PosInt g = s; g < hDiffHat.n_cols; ++g)
    {
        // this is equivalent to:
        // ret += hDiffHat.col(g) * arma::trans(hDiffHat.col(g - s));
        addKron(ret, hDiffHat.col(g), hDiffHat.col(g - s));
    }

    ret /= static_cast<double>(hDiffHat.n_cols);
    return ret;
}

double
getLogMargLikEst(const MyDoubleVector& numeratorTerms,
                 const MyDoubleVector& denominatorTerms,
                 double highDensityPointLogUnPosterior,
                 double& se,
                 PosInt endLag)
{
    // note that here it is assumed that there are equally many numerator and
    // denominator numbers, therefore
    PosInt nSamples = numeratorTerms.size();
    // assert(nSamples == denominatorTerms.size());

    // convert to rowvectors and compute respective means
    arma::rowvec numerator = arma::conv_to<arma::rowvec>::from(numeratorTerms);
    arma::rowvec denominator = arma::conv_to<arma::rowvec>::from(denominatorTerms);

    double hHat1 = arma::mean(numerator);
    double hHat2 = arma::mean(denominator);

    // so the estimate for the log posterior density ordinate is:
    double logPostDensityEstimate = log(hHat1) - log(hHat2);

    // compute the standard error, if that is possible
    se = R_NaReal;
    if(nSamples >= endLag)
    {
        // combine the corresponding numerator and denominator (centered) samples
        AMatrix hDiffHat = arma::join_cols(numerator - hHat1,
                                           denominator - hHat2);

        // initialize the 2x2 covariance matrix with Omega_0:
        arma::mat22 covMat = getOmega(0, hDiffHat);

        // then add the other terms
        for(PosInt s = 1; s <= endLag; ++s)
        {
            arma::mat22 thisOmega = getOmega(s, hDiffHat);
            covMat += (1.0 - s / (endLag + 1.0)) * (thisOmega + arma::trans(thisOmega));
        }

        // now covMat is ready to get delta ruled:
        arma::vec2 aHat;
        aHat.at(0) = 1.0 / hHat1;
        aHat.at(1) = - 1.0 / hHat2;

        se = sqrt(arma::as_scalar(arma::trans(aHat) * covMat * aHat / nSamples));
    }

    return highDensityPointLogUnPosterior - logPostDensityEstimate;
}
