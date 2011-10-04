/*
 * glmGetSamples.cpp
 *
 *  Created on: 17.03.2011
 *      Author: daniel
 *
 * MCMC sampling for a specific GLM configuration.
 *
 */

#include <rcppExport.h>

#include <dataStructure.h>
#include <types.h>
#include <iwls.h>
#include <random.h>
#include <zdensity.h>
#include <linalgInterface.h>
#include <getLogMargLikEst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// Helper class encapsulating an MCMC state
class McmcState
{
public:
    // ctr
    McmcState(const ZdensApprox& marginalz,
              PosInt nObs,
              PosInt nCoefs) :
                  sample(nCoefs),
                  proposalInfo(nObs, nCoefs),
                  marginalz(marginalz)
    {
    }

    // the current parameter sample
    Parameter sample;

    // the unnormalized log posterior of this sample
    double logUnPosterior;

    // info about the normal proposal distribution given the sampled z
    IwlsResults proposalInfo;

    // compute the log of the normalized proposal density when the z log density is provided
    // normalize correctly in order not to get problems (perhaps) with the Chib-Jeliazkov estimate
    // computation.
    double
    computeLogProposalDens() const
    {
        // Be careful: qFactor is in fact lower-triangular, so a simple multiplication would fail!
        // use instead directly a BLAS routine for this multiplication.
        AVector tmp = sample.coefs - proposalInfo.coefs;
        trmv(false, true, proposalInfo.qFactor, tmp);

        return 0.5 * (proposalInfo.logPrecisionDeterminant - arma::dot(tmp, tmp)) -
               M_LN_SQRT_2PI * proposalInfo.qFactor.n_rows +
               marginalz.logDens(sample.z);
    }

    // non-default assignment operator
    McmcState&
    operator=(const McmcState& rhs)
    {
        if(this == &rhs)
        {
            return *this;
        }
        else
        {
            sample = rhs.sample;
            logUnPosterior = rhs.logUnPosterior;
            proposalInfo = rhs.proposalInfo;

            return *this;
        }
    }


private:
    // the marginal z info: same for all Mcmc objects,
    // therefore it is not assigned by the assignment operator
    const ZdensApprox marginalz;
    // important: copy the object, because otherwise (reference/pointer)
    // we are not sure that the functions are still available if we do not use "new"
};


// Store the MCMC samples
struct Samples
{
    // constructor: basically allocates beta matrix
    Samples(const ModelPar& mod,
            const GlmModelData& modelData,
            PosInt nSamples) :
                nSamples(nSamples),
                nSaved(0)
    {
        // translate config into linear and splines part:
        for(int i_cov = 0; i_cov < modelData.nCovs; ++i_cov)
        {
            const double d = mod.config[i_cov];

            if(d > 0.0)
            {
                whichLinear.push_back(i_cov);

                if(d > 1.0)
                {
                    whichSpline.push_back(i_cov);
                }
            }
        }

        // corresponding number of linear and spline effects
        PosInt dimLinear = whichLinear.size();
        PosInt dimSpline = whichSpline.size();
        PosInt dimZ = dimSpline * modelData.dimSplineBasis;

        // so the number of coefficients altogether is:
        PosInt nCoefs = 1 + dimLinear + dimZ;

        // allocate the matrix
        coefsSamples.set_size(nCoefs, nSamples);
    }

    // save a sample consisting of coefs and z
    void
    storeParameters(const Parameter& sample)
    {
        coefsSamples.col(nSaved++) = sample.coefs;
        zSamples.push_back(sample.z);
    }

    // save terms for marginal likelihood estimate
    void
    storeMargLikTerms(double num, double denom)
    {
        numerator.push_back(num);
        denominator.push_back(denom);
    }

    // output everything to an R list
    List
    convert2list(const GlmModelData& modelData) const
    {
        // compute t = g / (g + 1) samples:
        DoubleVector tSamples;
        for(DoubleVector::const_iterator
                z = zSamples.begin();
                z != zSamples.end();
                ++z)
        {
            tSamples.push_back(Rf_plogis(*z, 0.0, 1.0, TRUE, FALSE));
        }

        // now process the big samples matrix.
        // this variable records which row to extract next:
        PosInt nextRow = 0;

        // extract intercept samples, the first row of the big samples matrix:
        NumericVector interceptSamples = wrap(arma::rowvec(coefsSamples.row(nextRow++)));
        interceptSamples.attr("dim") = Dimension(nSamples);

        // for the linear and spline coefs list, we need the names
        // of the covariates
        List dimnames = modelData.Xfull.attr("dimnames");
        StringVector allNames = dimnames[1]; // colnames of Xfull

        // extract linear coefficients:
        // one numeric vector for each in a list.
        List linearCoefsSamples;
        StringVector linearCoefsSamplesNames;

        for(PosInt i = 0; i < whichLinear.size(); ++i)
        {
            NumericVector thisVec = wrap(arma::rowvec(coefsSamples.row(nextRow++)));
            thisVec.attr("dim") = Dimension(nSamples);

            linearCoefsSamples.push_back(thisVec);
            linearCoefsSamplesNames.push_back(allNames[whichLinear[i]]);
        }
        linearCoefsSamples.names() = linearCoefsSamplesNames;

        // extract linear coefficients:
        // one numeric matrix for each in a list.
        List splineCoefsSamples;
        StringVector splineCoefsSamplesNames;

        for(PosInt i = 0; i < whichSpline.size(); ++i)
        {
            // offset 1 + dimLinear bc. int. and lin coefs are in front:
            AMatrix thisMat = coefsSamples.rows(nextRow, nextRow + modelData.dimSplineBasis - 1);
            nextRow += modelData.dimSplineBasis;

            splineCoefsSamples.push_back(thisMat);
            splineCoefsSamplesNames.push_back(allNames[whichSpline[i]]);
        }
        splineCoefsSamples.names() = splineCoefsSamplesNames;

        // return the list
        return List::create(_["t"] = tSamples,
                            _["intercept"] = interceptSamples,
                            _["linearCoefs"] = linearCoefsSamples,
                            _["splineCoefs"] = splineCoefsSamples,
                            _["z"] = zSamples);
    }

    // nCoefs x nSamples:
    PosInt nSamples;
    AMatrix coefsSamples;

    // counts the number of saved parameters, so that we know where to store the
    // next coefficients vector
    PosInt nSaved;

    // is gradually extended:
    DoubleVector zSamples;

    // possibly stays empty if not required by the user:
    // the numerator and denominator terms for the marginal likelihood estimate
    DoubleVector numerator;
    DoubleVector denominator;

    // which covariates are included as linear and spline effects?
    // this also translates to the ordering of the design matrix!
    IntVector whichLinear;
    IntVector whichSpline;
};


// R-interface to get the MCMC samples

// 19/5/2011 checked that no changes are necessary for offsets
SEXP
cpp_glmGetSamples(SEXP R_modelConfig,
                  SEXP R_modelData,
                  SEXP R_mcmc,
                  SEXP R_computation)
{
    // convert the arguments:
    ModelPar modPar(as<NumericVector>(R_modelConfig));
    GlmModelData modelData(R_modelData);
    McmcOptions mcmc(R_mcmc);
    Computation computation(R_computation);

    // use only one thread if we do not want to use openMP.
#ifdef _OPENMP
    if(! computation.useOpenMP)
    {
        omp_set_num_threads(1);
    } else {
        omp_set_num_threads(omp_get_num_procs());
    }
#endif


    // ----------------------------------------------------------------------------------
    // prepare the sampling
    // ----------------------------------------------------------------------------------

    // get the marginal z density approximation for this model
    ZdensApprox marginalz(modPar,
                          modelData,
                          computation);

    // construct IWLS object, which can be used for all IWLS stuff,
    // and also contains the design matrix etc
    Iwls iwlsObject(modPar,
                    modelData,
                    false,
                    EPS,
                    computation.debug);

    // allocate sample container
    Samples samples(modPar,
                    modelData,
                    mcmc.samples);

    // count how many proposals we have accepted:
    PosInt nAccepted(0);

    // at what z do we start?
    double startZ = marginalz.getzMode();

    // get the mode for beta given the mode of the approximated marginal posterior as z
    PosInt iwlsIterations = iwlsObject.startWithNewLinPred(30,
                                                           // this is the corresponding g
                                                           exp(startZ),
                                                           // and the start value for the linear predictor is taken from the Glm model config
                                                           modelData.linPredStart);

    // echo debug-level message?
    if(computation.debug)
    {
        Rprintf("\ncpp_sampleGlm: Initial IWLS for high density point finished after %d iterations",
                iwlsIterations);
    }

    // start container with current things
    McmcState now(marginalz,
                  modelData.nObs,
                  iwlsObject.getnCoefs());

    // this is the current proposal info:
    now.proposalInfo = iwlsObject.getResults();

    // and this is the current parameters sample:
    now.sample = Parameter(now.proposalInfo.coefs,
                           startZ);

    // compute the (unnormalized) log posterior of the proposal
    now.logUnPosterior = iwlsObject.computeLogUnPosteriorDens(now.sample);

    // so the parameter object "now" is then also the high density point
    // required for the marginal likelihood estimate:
    const McmcState highDensityPoint(now);

    // we accept this starting value, so initialize "old" with the same ones
    McmcState old(now);

    // ----------------------------------------------------------------------------------
    // start sampling
    // ----------------------------------------------------------------------------------

    // progress bar setup
    StatusBar statusBar(computation.verbose, 20, 5, mcmc.iterations);


    // echo debug-level message?
    if(computation.debug)
    {
        Rprintf("\ncpp_sampleGlm: Starting MCMC loop");
    }


    // loop
    for(PosInt i_iter = 0; i_iter < mcmc.iterations; statusBar.showProgress(++i_iter))
    {
        // echo debug-level message?
        if(computation.debug)
        {
            Rprintf("\ncpp_sampleGlm: Starting iteration no. %d", i_iter);
        }

        // ----------------------------------------------------------------------------------
        // store the proposal
        // ----------------------------------------------------------------------------------

        try
        {

        // sample one new log covariance factor z
        now.sample.z = marginalz.sample();

        // then do 1 IWLS step, starting from the last linear predictor and the new z
        // (here the return value is not very interesting, as it must be 1)
        iwlsObject.startWithNewCoefs(mcmc.nIwlsIterations,
                                     exp(now.sample.z),
                                     now.sample.coefs);

        // get the results
        now.proposalInfo = iwlsObject.getResults();

        // draw the proposal coefs:
        now.sample.coefs = drawNormalVector(now.proposalInfo.coefs,
                                            now.proposalInfo.qFactor);

        // compute the (unnormalized) log posterior of the proposal
        now.logUnPosterior = iwlsObject.computeLogUnPosteriorDens(now.sample);

        // ----------------------------------------------------------------------------------
        // get the reverse jump normal density
        // ----------------------------------------------------------------------------------

        // copy the old Mcmc object
        McmcState reverse(old);

        // do again 1 IWLS step, starting from the sampled linear predictor and the old z
        iwlsObject.startWithNewCoefs(mcmc.nIwlsIterations,
                                     exp(reverse.sample.z),
                                     now.sample.coefs);

        // get the results for the reverse jump Gaussian:
        // only the proposal has changed in contrast to the old container,
        // the sample stays the same!
        reverse.proposalInfo = iwlsObject.getResults();


        // ----------------------------------------------------------------------------------
        // compute the proposal density ratio
        // ----------------------------------------------------------------------------------

        // first the log of the numerator, i.e. log(f(old | new)):
        double logProposalRatioNumerator = reverse.computeLogProposalDens();

        // second the log of the denominator, i.e. log(f(new | old)):
        double logProposalRatioDenominator = now.computeLogProposalDens();

        // so the log proposal density ratio is
        double logProposalRatio = logProposalRatioNumerator - logProposalRatioDenominator;

        // ----------------------------------------------------------------------------------
        // compute the posterior density ratio
        // ----------------------------------------------------------------------------------

        double logPosteriorRatio = now.logUnPosterior - old.logUnPosterior;

        // ----------------------------------------------------------------------------------
        // accept or reject proposal
        // ----------------------------------------------------------------------------------

        double acceptanceProb = exp(logPosteriorRatio + logProposalRatio);

        if(unif() < acceptanceProb)
        {
            old = now;

            ++nAccepted;
        }
        else
        {
            now = old;
        }

        // ----------------------------------------------------------------------------------
        // store the sample?
        // ----------------------------------------------------------------------------------

        // if the burnin was passed and we are at a multiple of step beyond that, then store
        // the sample.
        if((i_iter >= mcmc.burnin) &&
           (((i_iter + 1 - mcmc.burnin) % mcmc.step) == 0))
        {
            // echo debug-level message
            if(computation.debug)
            {
                Rprintf("\ncpp_sampleGlm: Storing samples of iteration no. %d", i_iter);
            }

            // store the current parameter sample
            samples.storeParameters(now.sample);

            // ----------------------------------------------------------------------------------
            // compute marginal likelihood terms
            // ----------------------------------------------------------------------------------

            // echo debug-level message?
            if(computation.debug)
            {
                Rprintf("\ncpp_sampleGlm: Compute marginal likelihood estimation terms");
            }

            // ----------------------------------------------------------------------------------
            // compute next term for the denominator
            // ----------------------------------------------------------------------------------

            // draw from the high density point proposal distribution
            McmcState denominator(highDensityPoint);
            denominator.sample.z = marginalz.sample();

            iwlsObject.startWithNewLinPred(mcmc.nIwlsIterations,
                                           exp(denominator.sample.z),
                                           highDensityPoint.proposalInfo.linPred);

            denominator.proposalInfo = iwlsObject.getResults();

            denominator.sample.coefs = drawNormalVector(denominator.proposalInfo.coefs,
                                                        denominator.proposalInfo.qFactor);

            // get posterior density of the sample
            denominator.logUnPosterior = iwlsObject.computeLogUnPosteriorDens(denominator.sample);

            // get the proposal density at the sample
            double denominator_logProposalDensity = denominator.computeLogProposalDens();

            // then the reverse stuff:
            // first we copy again the high density point
            McmcState revDenom(highDensityPoint);

            // but choose the new sampled coefficients as starting point
            iwlsObject.startWithNewCoefs(mcmc.nIwlsIterations,
                                         exp(revDenom.sample.z),
                                         denominator.sample.coefs);
            revDenom.proposalInfo = iwlsObject.getResults();

            // so the reverse proposal density is
            double revDenom_logProposalDensity = revDenom.computeLogProposalDens();


            // so altogether the next term for the denominator is the following acceptance probability
            double denominatorTerm = denominator.logUnPosterior - highDensityPoint.logUnPosterior +
                    revDenom_logProposalDensity - denominator_logProposalDensity;
            denominatorTerm = exp(fmin(0.0, denominatorTerm));

            // ----------------------------------------------------------------------------------
            // compute next term for the numerator
            // ----------------------------------------------------------------------------------

            // compute the proposal density of the current sample starting from the high density point
            McmcState numerator(now);

            iwlsObject.startWithNewLinPred(mcmc.nIwlsIterations,
                                           exp(numerator.sample.z),
                                           highDensityPoint.proposalInfo.linPred);
            numerator.proposalInfo = iwlsObject.getResults();

            double numerator_logProposalDensity = numerator.computeLogProposalDens();

            // then compute the reverse proposal density of the high density point when we start from the current
            // sample
            McmcState revNum(highDensityPoint);

            iwlsObject.startWithNewCoefs(mcmc.nIwlsIterations,
                                         exp(revNum.sample.z),
                                         now.sample.coefs);
            revNum.proposalInfo = iwlsObject.getResults();

            double revNum_logProposalDensity = revNum.computeLogProposalDens();

            // so altogether the next term for the numerator is the following guy:
            double numeratorTerm = exp(fmin(revNum_logProposalDensity,
                                            highDensityPoint.logUnPosterior - now.logUnPosterior +
                                            numerator_logProposalDensity));

            // ----------------------------------------------------------------------------------
            // finally store both terms
            // ----------------------------------------------------------------------------------

            samples.storeMargLikTerms(numeratorTerm, denominatorTerm);

        }

        } // try
        catch (std::domain_error error)
        {
            Rf_error("Error in MCMC iteration %d because\n%s",
                     i_iter, error.what());
        }

        // echo debug-level message?
        if(computation.debug)
        {
            Rprintf("\ncpp_sampleGlm: Finished iteration no. %d", i_iter);
        }

    } // end MCMC loop


    // echo debug-level message?
    if(computation.debug)
    {
        Rprintf("\ncpp_sampleGlm: Finished MCMC loop");
    }


    // ----------------------------------------------------------------------------------
    // build up return list for R and return that.
    // ----------------------------------------------------------------------------------

    double se;
    double logMargLikEst = getLogMargLikEst(samples.numerator,
                                            samples.denominator,
                                            highDensityPoint.logUnPosterior,
                                            se);

    return List::create(_["samples"] = samples.convert2list(modelData),
                        _["mcmc"] = List::create(_["nAccepted"] = nAccepted,
                                                 _["acceptanceRatio"] = nAccepted * 1.0 / mcmc.iterations),
                        _["logMargLik"] = List::create(_["ilaEstimate"] = marginalz.getLogMargLik(),
                                                       _["mcmcEstimate"] = logMargLikEst,
                                                       _["mcmcSe"] = se,
                                                       _["margApproxZdens"] = marginalz.getLinApproxDens()));
}


