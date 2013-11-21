/*
 *  glmExhaustive.cpp
 *
 *  Created on: 10.03.2011
 *      Author: daniel
 *
 *  This is the GLM version for exhaustive model space evaluation.
 *
 */

#include <rcppExport.h>
#include <dataStructure.h>
#include <zdensity.h>
#include <sum.h>

using namespace Rcpp;

// 31/01/2013: add R_modelPrior option. Inclusion probabilities are computed in R.

// R-interface for exhaustive evaluation of all models provided
SEXP
cpp_glmExhaustive(SEXP R_modelData,
                  SEXP R_modelPrior,
                  SEXP R_modelConfigs,
                  SEXP R_computation)
{
    // convert and extract the modelData object:
    GlmModelData modelData(R_modelData);

    // the prior on the model space:
    ModelPrior modelPrior(R_modelPrior, modelData.nCovs, modelData.nDegrees, modelData.continuous);

    // convert the model configs:
    IntegerMatrix modelConfigs(R_modelConfigs);
    const int nModels(modelConfigs.nrow());

    // and convert the computation options
    Computation computation(R_computation);

    // setup containers for results
    MyDoubleVector logMargLik;
    MyDoubleVector logPrior;
    MyDoubleVector logPost;

    // this is needed for computing the normalising constant
    // of the posterior probabilities:
    SafeSum vec;

    // progress bar setup:
    StatusBar statusBar(computation.verbose, 20, 5, nModels);

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

        // and compute log marginal likelihood
        double lm = glmGetLogMargLik(modelPar,
                                     modelData,
                                     computation);
        logMargLik.push_back(lm);

        // and the log prior
        double lp = modelPrior.getLogPrior(modelPar);
        logPrior.push_back(lp);

        // so the unnormalised log posterior is
        double ulp = lm + lp;
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
    return DataFrame::create(_["logMargLik"] = logMargLik,
                             _["logPrior"] = logPrior,
                             _["post"] = logPost);
}
