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

#include <dataStructure.h>

using namespace Rcpp;

SEXP
cpp_exhaustive2(SEXP R_modelData,
                SEXP R_modelConfigs)
{
    // convert the modelData object:
    ModelData modelData(R_modelData);

     // convert and extract the model configs:
    IntegerMatrix modelConfigs(R_modelConfigs);
    const int nModels(modelConfigs.nrow());

    // setup containers for results
    std::vector<double> R2;
    std::vector<double> logMargLik;

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
        logMargLik.push_back(modelData.getLogMargLik(modelPar, thisR2));
        R2.push_back(thisR2);
    }

    // return the results
    return DataFrame::create(_["R2"] = R2,
                             _["logMargLik"] = logMargLik);
}
