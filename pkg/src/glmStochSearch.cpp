/*
 * glmStochSearch.cpp
 *
 *  Created on: 04.04.2011
 *      Author: daniel
 */

#include <rcppExport.h>
#include <vector>
#include <algorithm>

#include <types.h>
#include <dataStructure.h>
#include <random.h>
#include <zdensity.h>


using namespace Rcpp;


// todo: Note that this function is almost identical to the normal response version....
// This is clear but not used at the moment.
// 02/05/2012: it is no longer mandatory to have at least
//             two continuous covariates. For 0 or 1 continuous
//             covariates the move step is always chosen.
// 10/05/2012: allow flexible starting model
SEXP
cpp_glmStochSearch(SEXP R_modelData,
                   SEXP R_modelPrior,
                   SEXP R_searchSettings,
                   SEXP R_computation)
{
    // convert and extract the modelData object:
    GlmModelData modelData(R_modelData);

    // determine probability of a "move" step
    const double moveProb = (modelData.nContCovs < 2) ? 1.01 : 0.75;

    // the prior on the model space:
    ModelPrior modelPrior(R_modelPrior, modelData.nCovs, modelData.nDegrees, modelData.continuous);

    // the search settings in this space:
    List searchSettings(R_searchSettings);
    ModelPar startModPar(as<NumericVector>(searchSettings["startModel"]));
    startModPar.compDegIndex(modelData.degrees);
    const int chainlength = as<int>(searchSettings["chainlength"]);
    const int nCache = as<int>(searchSettings["nCache"]);
    const int nModels = as<int>(searchSettings["nModels"]);

    // and convert the computation options
    Computation computation(R_computation);

    // setup cache infrastructure
    ModelCache modelCache(nCache);

    // compute start model properties
    const double logMargLik = glmGetLogMargLik(startModPar, modelData, computation);
    const double logPrior = modelPrior.getLogPrior(startModPar);

    ModelInfo modInfo(logMargLik, logPrior);

    // insert into the cache
    modelCache.insert(startModPar, modInfo);

    // put the start model in the mcmc containers
    ModelMcmc old(startModPar, logMargLik, logPrior);
    ModelMcmc now(old);

    // progress bar setup
    StatusBar statusBar(computation.verbose, 20, 5, chainlength);

    // Start MCMC sampler***********************************************************//

    // use R's random number generator
    GetRNGstate();

    for(int t = 0; t < chainlength; statusBar.showProgress(++t))
    {
        // check if any interrupt signals have been entered
        R_CheckUserInterrupt();

        // log(proposal ratio)
        double logPropRatio;

        // randomly select move type: move or switch.
        if (unif_rand() < moveProb) // this may also be a search option in the future
        {
            // ****************************** //
            // implement a move

            // first randomly select a covariate index
            int covIndex = discreteUniform(0, modelData.nCovs);

            // differentiate between continuous and "linear" covariate
            if(modelData.continuous[covIndex])
            {
                // where are we in the degrees table?
                int& degIndex = now.modPar.degIndex.at(covIndex);
                int oldDegIndex = degIndex;

                // decide the move
                bool goDown;
                if(degIndex == modelData.nDegrees - 1)
                {
                    goDown = true;
                }
                else if (degIndex == 0)
                {
                    goDown = false;
                }
                else
                {
                    goDown = (unif_rand() < 0.5);
                }

                // do the move
                if(goDown)
                {
                    now.modPar.config.at(covIndex) = modelData.degrees.at(--degIndex);
                }
                else // go up
                {
                    now.modPar.config.at(covIndex) = modelData.degrees.at(++degIndex);
                }

                // calculate the log proposal ratio
                if(oldDegIndex == (modelData.nDegrees - 1) ||
                   oldDegIndex == 0)
                {
                    logPropRatio = M_LN2;
                }
                else if(degIndex == (modelData.nDegrees - 1) ||
                        degIndex == 0)
                {
                    logPropRatio = - M_LN2;
                }
                else
                {
                    logPropRatio = 0.0;
                }
            }
            else // "linear" covariate
            {
                // move from 0 to 1 or from 1 to 0 via:
                now.modPar.config[covIndex] = 1 - now.modPar.config[covIndex];
                now.modPar.degIndex[covIndex] = 1 - now.modPar.degIndex[covIndex];

                logPropRatio = 0.0;
            }
        }
        else
        {
            // ****************************** //
            // implement a switch
            // (in the continuous covariates only!)

            // determine the two covariates where we want to switch:

            int firstCov = discreteUniform(0, modelData.nContCovs);

            int secondCov = discreteUniform(0, modelData.nContCovs - 1);
            if(secondCov >= firstCov)
                ++secondCov;

            firstCov = modelData.contCovs[firstCov];
            secondCov = modelData.contCovs[secondCov];

            // do the switch:

            now.modPar.config.at(firstCov) = old.modPar.config.at(secondCov);
            now.modPar.degIndex.at(firstCov) = old.modPar.degIndex.at(secondCov);

            now.modPar.config.at(secondCov) = old.modPar.config.at(firstCov);
            now.modPar.degIndex.at(secondCov) = old.modPar.degIndex.at(firstCov);

            // and the log proposal ratio is easy
            logPropRatio = 0.0;
        }

        // search for info of the proposed model
        bool nowWasFound = modelCache.getModelInfo(now.modPar, modInfo);

        if (nowWasFound)
        {
            // "now" is an old model, so get the necessary info from there.
            now.logMargLik = modInfo.logMargLik;
            now.logPrior = modInfo.logPrior;
        }
        else
        { // "now" is a new model

            // so we must compute the log marg lik now.
            now.logMargLik = glmGetLogMargLik(now.modPar, modelData, computation);

            // check if the new model is OK
            if (R_IsNaN(now.logMargLik))
            {
                // we do not save this model in the model cache
                modelCache.nanCounter++;
            }
            else
            { // OK => then compute the rest, and insert into model cache

                now.logPrior = modelPrior.getLogPrior(now.modPar);

                // insert the model parameter/info into the model cache

                // problem: this could erase the old model from the model cache,
                // and invalidate the iterator old.mapPos!
                // ==> so we cannot work with the iterators here.
                modelCache.insert(now.modPar,
                                  ModelInfo(now.logMargLik,
                                            now.logPrior));
            }
        }

        // decide acceptance:
        // for acceptance, the new model must be valid and the acceptance must be sampled
        if ((R_IsNaN(now.logMargLik) == FALSE) &&
                (log(unif_rand()) <= (now.logMargLik - old.logMargLik +
                                      now.logPrior - old.logPrior +
                                      logPropRatio)))
        { // acceptance
            old = now;
        }
        else
        { // rejection
            now = old;
        }

        // so now definitely old == now, and we can
        // increment the associated sampling frequency
        // (if the model is not in the cache because it was too bad this
        // is not a problem)
        modelCache.incrementFrequency(now.modPar);
    }
    Rprintf("\n");

    PutRNGstate(); // no RNs required anymore

    // Finished MCMC sampler********************************************************//

    // normalize posterior probabilities
    const long double logNormConst = modelCache.getLogNormConstant();

    // get the nModels best models from the cache as an R list
    List ret = modelCache.getListOfBestModels(nModels, logNormConst);

    // set the attributes
    ret.attr("numVisited") = modelCache.size();
    ret.attr("inclusionProbs") = modelCache.getInclusionProbs(logNormConst, modelData.nCovs);
    ret.attr("logNormConst") = logNormConst;

    // Print statistics
    Rprintf("\nNumber of non-identifiable model proposals:     %d", modelCache.nanCounter);
    Rprintf("\nNumber of total cached models:                  %d", modelCache.size());
    Rprintf("\nNumber of returned models:                      %d\n", ret.length());

    // return
    return ret;
}



