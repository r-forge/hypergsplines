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


using namespace Rcpp;


// R-interface for exhaustive evaluation of all models provided
SEXP
cpp_glmExhaustive(SEXP R_modelData,
                  SEXP R_modelConfigs,
                  SEXP R_computation)
{
    // convert and extract the modelData object:
    GlmModelData modelData(R_modelData);

    // convert the model configs:
    IntegerMatrix modelConfigs(R_modelConfigs);
    const int nModels(modelConfigs.nrow());

    // and convert the computation options
    Computation computation(R_computation);

    // setup containers for results
    DoubleVector logMargLik;

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
        logMargLik.push_back(glmGetLogMargLik(modelPar,
                                              modelData,
                                              computation));
    }

    // return the results
    return DataFrame::create(_["logMargLik"] = logMargLik);
}
