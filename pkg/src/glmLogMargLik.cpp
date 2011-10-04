/*
 * logMargLik.cpp
 *
 *  Created on: 21.06.2011
 *      Author: daniel
 */

#include <rcppExport.h>
#include <dataStructure.h>
#include <zdensity.h>

using namespace Rcpp;

SEXP
cpp_glmLogMargLik(SEXP R_modelConfig,
                  SEXP R_modelData,
                  SEXP R_computation)
{
    // translate config into ModelPar
    ModelPar modelPar(as<NumericVector>(R_modelConfig));

    // convert the modelData object:
    GlmModelData modelData(R_modelData);

    // and convert the computation options
    Computation computation(R_computation);

    // compute log marginal likelihood:
    const double ret = glmGetLogMargLik(modelPar,
                                        modelData,
                                        computation);

    return wrap(ret);
}

