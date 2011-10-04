/*
 * rcppExport.h
 *
 *  Created on: 13.10.2010
 *      Author: daniel
 */

#ifndef RCPPEXPORT_H_
#define RCPPEXPORT_H_

#include <RcppArmadillo.h>

RcppExport SEXP cpp_calculateModel(SEXP R_config,
                                   SEXP R_modelData);

RcppExport SEXP cpp_exhaustive(SEXP R_modelData,
                               SEXP R_modelConfigs);

RcppExport SEXP cpp_exhaustive2(SEXP R_modelData,
                                SEXP R_modelConfigs);

RcppExport SEXP cpp_stochSearch(SEXP R_modelData,
                                SEXP R_modelPrior,
                                SEXP R_searchSettings);

RcppExport SEXP cpp_glmStochSearch(SEXP R_modelData,
                                   SEXP R_modelPrior,
                                   SEXP R_searchSettings,
                                   SEXP R_computation);

RcppExport SEXP cpp_glmExhaustive(SEXP R_modelData,
                                  SEXP R_modelConfigs,
                                  SEXP R_computation);

RcppExport SEXP cpp_glmGetSamples(SEXP R_modelConfig,
                                  SEXP R_modelData,
                                  SEXP R_mcmc,
                                  SEXP R_computation);

RcppExport SEXP cpp_linApproxDens(SEXP args,
                                  SEXP logDens,
                                  SEXP grid,
                                  SEXP nSamples);

RcppExport SEXP cpp_hyp2f1(SEXP R_interface);

RcppExport SEXP cpp_log_hyp2f1_laplace(SEXP R_interface);

RcppExport SEXP cpp_logMargLik(SEXP R_modelConfig,
                               SEXP R_modelData);

RcppExport SEXP cpp_glmLogMargLik(SEXP R_modelConfig,
                                  SEXP R_modelData,
                                  SEXP R_computation);

RcppExport SEXP cpp_aggregateModelsTable(SEXP R_modelsTable,
                                  SEXP R_posterior,
                                  SEXP R_cut);

#endif /* RCPPEXPORT_H_ */
