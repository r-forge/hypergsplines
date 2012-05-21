/*
 * calculateModel.cpp
 *
 *  Created on: 13.10.2010
 *      Author: daniel
 */

#include <rcppExport.h>
// #include <backsolve.h>
#include <linalgInterface.h>
#include <getRho.h>

#include <algorithm>

using namespace Rcpp;


SEXP cpp_calculateModel(SEXP R_config,
                        SEXP R_modelData)
{
    // convert to Rcpp types
    NumericVector config(R_config);
    List modelData(R_modelData);

    // extract dimensions
    int nObs = as<int>(modelData["nObs"]);
    int nCovs = as<int>(modelData["nCovs"]);

    // translate config into linear and splines part
    IntegerVector whichLinear;
    IntegerVector whichSpline;

    for(int i=0; i < nCovs; ++i)
    {
        const double d = config(i);

        if(d > 0.0)
        {
            whichLinear.push_back(i);

            if(d > 1.0)
            {
                whichSpline.push_back(i);
            }
        }
    }

    // is this the null model / only linear effects?
    const int dimLinear = whichLinear.length();
    const bool isNullModel = (dimLinear == 0);

    const bool hasOnlyLinear = (whichSpline.length() == 0);

    // the original response
    NumericVector ytemp = modelData["y"];
    arma::colvec y(ytemp.begin(), nObs);

    // the mean response vector
    arma::colvec meanY(nObs);
    meanY.fill(arma::mean(y));

    // construct the linear design matrix Xlin (without the intercept):
    NumericMatrix Xfulltemp = modelData["X"];
    const double * Xfull_ptr = Xfulltemp.begin();

    arma::mat Xlin(nObs, dimLinear);
    double * Xlin_ptr = Xlin.memptr();

    for(int i = 0; i < dimLinear; ++i)
    {
        const int start = whichLinear(i) * nObs;

        std::copy(Xfull_ptr + start, // start source pointer
                  Xfull_ptr + start + nObs, // one beyond the end source pointer
                  Xlin_ptr + (i * nObs)); // start destination pointer
    }

    // possible root of covariance matrix
    arma::mat Vroot;

    if(! hasOnlyLinear)
    {
        List rhoList = modelData["rho.list"];
        List lambdasList = modelData["lambdas.list"];
        List ztcrossprodList = modelData["Z.tcrossprod.list"];

//        // get also the getRhos R-function
//        Environment package("package:hypergsplines");
//        Function getRhos(package.find("getRhos"));

        // construct the covariance matrix of the marginal model:
        arma::mat V(nObs, nObs);
        V.eye();

        for(IntegerVector::iterator
                s = whichSpline.begin();
                s != whichSpline.end();
                ++s)
        {
            // get the correct rho
            const double d = config(*s);
            double rho;

            if(d == floor(d))
            {
                // dNumber <- match(d - 1,
                //                  modelData$splineDegrees)
                const int dNumber = round(d - 2);

                const NumericVector rhos = rhoList[*s];
                rho = rhos[dNumber];
            }
            else
            {
                const NumericVector lambdas = lambdasList[*s];
                rho = getRho(lambdas, d - 1.0);
            }

            // and add to V with that factor
            const NumericMatrix Ztcrossprod = ztcrossprodList[*s];
            const arma::mat Zt(Ztcrossprod.begin(), nObs, nObs, false);

            V += rho * Zt;
        }

        // do a Cholesky decomposition of V
        // and scale everything with that:
        Vroot = arma::chol(V);

        trs(true, true, Vroot, y);
        trs(true, true, Vroot, meanY);
        trs(true, true, Vroot, Xlin);

//        y = backsolve(Vroot,
//                      y,
//                      true,
//                      true);
//
//        meanY = backsolve(Vroot,
//                          meanY,
//                          true,
//                          true);
//
//        Xlin = backsolve(Vroot,
//                         Xlin,
//                         true,
//                         true);
    }

    double coefR2 = 0;
    arma::colvec betaOLS;
    arma::colvec olsFit;

    if(! isNullModel)
    {
        // compute the OLS solution and get the coefficient of determination (R^2):
        betaOLS = arma::solve(Xlin, y);
        olsFit = meanY + Xlin * betaOLS;

        coefR2 = arma::as_scalar(arma::sum(arma::square(olsFit - meanY))) /
                arma::as_scalar(arma::sum(arma::square(y - meanY)));
    }

    // return the list with necessary info
    IntegerVector whichLinearOne = whichLinear + 1;
    IntegerVector whichSplineOne = whichSpline + 1;
    return List::create(_["isNullModel"] = isNullModel,
                        _["hasOnlyLinear"] = hasOnlyLinear,
                        _["whichLinear"] = whichLinearOne,
                        _["whichSpline"] = whichSplineOne,
                        _["X.lin"] = Xlin,
                        _["dim.lin"] = dimLinear,
                        _["y"] = y,
                        _["meanY"] = meanY,
                        _["Vroot"] = Vroot,
                        _["betaOLS"] = betaOLS,
                        _["coefR2"] = coefR2);
}
