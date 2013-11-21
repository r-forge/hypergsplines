/*
 * exhaustive.cpp
 *
 *  Created on: 14.10.2010
 *      Author: daniel
 *
 *  30/03/2012: remove hyperparameter a, instead use string for g-prior
 */

#include <rcppExport.h>
#include <backsolve.h>
#include <vector>
#include <string>
#include <types.h>
#include <sum.h>

#include <dataStructure.h>
#include <logBFhypergn.h>
#include <logBFhyperg.h>

using namespace Rcpp;

// 31/01/2013: add R_modelPrior option. Inclusion probabilities are computed in R.

SEXP
cpp_exhaustive(SEXP R_modelData,
               SEXP R_modelPrior,
               SEXP R_modelConfigs)
{
    // convert and extract the modelData object:
    List modelData(R_modelData);

    List rhoList = modelData["rho.list"];
    List ztcrossprodList = modelData["Z.tcrossprod.list"];

    NumericMatrix Xfulltemp = modelData["X"];
    const double * Xfull_ptr = Xfulltemp.begin();

    const int nObs(as<int>(modelData["nObs"]));
    const int nCovs(as<int>(modelData["nCovs"]));

    const std::string gPriorString(as<std::string>(modelData["gPrior"]));

    const NumericVector ytemp = modelData["y"];
    const arma::colvec y_orig(ytemp.begin(), nObs, false);

    arma::colvec meanY_orig(nObs);
    meanY_orig.fill(arma::mean(y_orig));

    // the prior on the model space:
    ModelPrior modelPrior(R_modelPrior,
                          as<int>(modelData["nCovs"]),
                          as<IntVector>(modelData["degrees"]).size(),
                          as<LogicalVector>(modelData["continuous"]));

    // convert and extract the model configs:
    IntegerMatrix modelConfigs(R_modelConfigs);
    const int nModels(modelConfigs.nrow());

    // setup containers for results
    MyDoubleVector R2;
    MyDoubleVector logMargLik;MyDoubleVector logPrior;
    MyDoubleVector logPost;

    // this is needed for computing the normalising constant
    // of the posterior probabilities:
    SafeSum vec;

    // progress bar setup:
    StatusBar statusBar(true, 20, 5, nModels);

    // now process each model configuration:
    for(int i_mod = 0; i_mod < nModels; statusBar.showProgress(++i_mod))
    {
        // check if any interrupt signals have been entered
        R_CheckUserInterrupt();

        // get this model config
        // IntegerMatrix::Row config = modelConfigs.row(i_mod);
        IntegerVector config = modelConfigs.row(i_mod);

        // translate config into ModelPar
        ModelPar modelPar(config);

        // translate config into linear and splines part
        IntVector whichLinear;
        IntVector whichSpline;

        for(int i_cov = 0; i_cov < nCovs; ++i_cov)
        {
            const int d = config[i_cov];

            if(d > 0)
            {
                whichLinear.push_back(i_cov);

                if(d > 1)
                {
                    whichSpline.push_back(i_cov);
                }
            }
        }

        // is this the null model / only linear effects?
        const int dimLinear = whichLinear.size();
        const bool isNullModel = (dimLinear == 0);

        const bool hasOnlyLinear = (whichSpline.size() == 0);

        // the original response
        arma::colvec y = y_orig;

        // the mean response vector
        arma::colvec meanY = meanY_orig;

        // construct the linear design matrix Xlin (without the intercept):
        arma::mat Xlin(nObs, dimLinear);
        double * Xlin_ptr = Xlin.memptr();

        for(int i_lin = 0; i_lin < dimLinear; ++i_lin)
        {
            const int start = whichLinear.at(i_lin) * nObs;

            std::copy(Xfull_ptr + start, // start source pointer
                      Xfull_ptr + start + nObs, // one beyond the end source pointer
                      Xlin_ptr + (i_lin * nObs)); // start destination pointer
        }

        // possible root of covariance matrix
        arma::mat Vroot;

        if(! hasOnlyLinear)
        {
            // construct the covariance matrix of the marginal model:
            arma::mat V = arma::eye(nObs, nObs);

            for(std::vector<int>::const_iterator
                    s = whichSpline.begin();
                    s != whichSpline.end();
                    ++s)
            {
                // get the correct rho
                const int d = config[*s];
                // dNumber <- match(d - 1,
                //                  modelData$splineDegrees)
                const int dNumber = d - 2;

                NumericVector rhos = rhoList[*s];
                const double rho = rhos[dNumber];

                // and add to V with that factor
                const NumericMatrix Ztcrossprod = ztcrossprodList[*s];
                const arma::mat Zt(Ztcrossprod.begin(), nObs, nObs, false);

                V += rho * Zt;
            }

            // do a Cholesky decomposition of V
            // and scale everything with that:
            Vroot = arma::chol(V);

            y = backsolve(Vroot,
                          y,
                          true,
                          true);

            meanY = backsolve(Vroot,
                              meanY,
                              true,
                              true);

            Xlin = backsolve(Vroot,
                             Xlin,
                             true,
                             true);
        }

        // sum of squares total:
        const double sst = arma::as_scalar(arma::sum(arma::square(y - meanY)));

        // start with R2 = 0, which is correct in the null model.
        double thisR2 = 0;
        if(! isNullModel)
        {
            // compute the OLS solution and get the coefficient of determination (R^2):
            arma::colvec betaOLS = arma::solve(Xlin, y);
            arma::colvec olsFit = meanY + Xlin * betaOLS;
            thisR2 = arma::as_scalar(arma::sum(arma::square(olsFit - meanY))) / sst;
        }
        R2.push_back(thisR2);

        // start the log marginal likelihood
        double thisLogMargLik = - (nObs - 1.0) / 2.0 * log(sst);

        // additional terms depend on the model structure:
        if(! isNullModel)
        {
            // and then add the log BF, depending on the hyperprior on g:
            if(gPriorString == "hyper-g/n")
            {
                thisLogMargLik += logBFhypergn(nObs, dimLinear, thisR2);
            }
            else // gPriorString == "hyper-g"
            {
                thisLogMargLik += logBFhyperg(nObs, dimLinear, thisR2);
            }

            // if splines were used, we must also take into account the
            // covariance matrix V for the marginal likelihood:
            if(! hasOnlyLinear)
            {
                thisLogMargLik -= arma::as_scalar(arma::sum(arma::log(Vroot.diag())));
            }
        }

        // save the log marginal likelihood
        logMargLik.push_back(thisLogMargLik);

        // and the log prior
        double lp = modelPrior.getLogPrior(modelPar);
        logPrior.push_back(lp);

        // so the unnormalised log posterior is
        double ulp = thisLogMargLik + lp;
        logPost.push_back(ulp);
        // add that to the SafeSum
        vec.add(ulp);

    }

    // log normalising constant for the posterior probs:
    double logNormConst = vec.logSumExp();

    // so we overwrite the unnormalised log posterior values
    // with the normalised posterior values:
    for(MyDoubleVector::iterator i = logPost.begin();
            i != logPost.end();
            ++i)
    {
        *i = exp(*i - logNormConst);
    }


    // return the results
    return DataFrame::create(_["R2"] = R2,
                             _["logMargLik"] = logMargLik,
                             _["logPrior"] = logPrior,
                             _["post"] = logPost);
}
