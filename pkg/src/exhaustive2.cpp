/*
 * exhaustive2.cpp
 *
 *  Created on: 15.10.2010
 *      Author: daniel
 *
 *  This is a second version of cpp_exhaustive, which shall be more
 *  efficient by avoiding factorization of an n x n matrix.
 *
 *  The results however, should be identical!!
 *
 *  4/4/2011: Depend on the ModelData class in dataStructure.h
 */

#include <rcppExport.h>
#include <vector>
#include <types.h>
#include <sum.h>

#include <dataStructure.h>

using namespace Rcpp;

// 31/01/2013: add R_modelPrior option. Inclusion probabilities are computed in R.

SEXP
cpp_exhaustive2(SEXP R_modelData,
                SEXP R_modelPrior,
                SEXP R_modelConfigs)
{
    // convert the modelData object:
    ModelData modelData(R_modelData);

    // the prior on the model space:
    ModelPrior modelPrior(R_modelPrior, modelData.nCovs, modelData.nDegrees, modelData.continuous);

     // convert and extract the model configs:
    IntegerMatrix modelConfigs(R_modelConfigs);
    const int nModels(modelConfigs.nrow());

    // setup containers for results
    MyDoubleVector R2;
    MyDoubleVector logMargLik;
    MyDoubleVector logPrior;
    MyDoubleVector logPost;

    // this is needed for computing the normalising constant
    // of the posterior probabilities:
    SafeSum vec;

    // progress bar setup:
    StatusBar statusBar(true, 20, 5, nModels);

    double thisR2 = 0;

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

        // compute log marginal likelihood and R2 and save them
        double lm = modelData.getLogMargLik(modelPar, thisR2);
        logMargLik.push_back(lm);
        R2.push_back(thisR2);

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
    return DataFrame::create(_["R2"] = R2,
                             _["logMargLik"] = logMargLik,
                             _["logPrior"] = logPrior,
                             _["post"] = logPost);
}
