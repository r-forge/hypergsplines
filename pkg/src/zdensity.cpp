#include <zdensity.h>

#include <rcppExport.h>
#include <stdexcept>

#include <linalgInterface.h>
#include <optimize.h>
#include <sum.h>
#include <linApproxDens.h>

using namespace Rcpp;


// *************************************************************************************
// Integrate a function and return the log of the approximation.
// You provide the negative log of the positive target function and computation
// options, and also get the mode and the variance estimate for the target function.
template <class Fun>
double
logIntegral(Fun& negLogFun,
            double& zMode,
            double& zVar,
            const Computation& computation)
{
    // construct an appropriate object for using the optimize routine
    Brent<Fun> brent(negLogFun,
                     -100.0, // log(DBL_MIN) + 40.0,
                     +200.0, // log(DBL_MAX) - 40.0,
                     sqrt(EPS));

    // and get the mode from that.
    zMode = brent.minimize();

    // here we have to compute the inverse hessian afterwards:
    // use the same epsilon here as for the minimization routine.
    AccurateNumericInvHessian<Fun> invHess(negLogFun);
    zVar = invHess(zMode);

    // be sure that we get a positive variance estimate
    if(zVar <= 0.0)
    {
        // warn
        Rf_warning("\nlogIntegral: Non-positive variance estimate %f encountered!\nProbably the optimization did not converge.\nResetting to default variance",
                   zVar);

        // set large default to explore some space
        zVar = 15.0;
    }

    // echo detailed progress in debug mode
    if(computation.debug)
    {
        Rprintf("\nlogIntegral: finished minimization at mode %f with variance %f", zMode, zVar);
    }


    // then compute the Gauss-Hermite quadrature, using the supplied standard nodes and
    // weights from R
    DoubleVector nodes;
    DoubleVector logWeights;

    // get the nodes and log weights for this mode and variance:
    computation.gaussHermite.getNodesAndLogWeights(zMode, zVar, nodes, logWeights);

    // the log contributions which will be stored here:
    SafeSum logContributions;

    // compute them now
    DoubleVector::const_iterator n = nodes.begin();
    for(DoubleVector::const_iterator
            w = logWeights.begin();
            w != logWeights.end();
            ++w, ++n)
    {
        logContributions.add((*w) - negLogFun(*n));
    }

    // the result is the log of the sum of the exp'ed values
    double ret = logContributions.logSumExp();

    // check finiteness
    if (! R_finite(ret))
    {
        std::ostringstream stream;
        stream << "logIntegral: got non-finite log integral approximation " << ret;
        throw std::domain_error(stream.str().c_str());
    }

    // return the result
    return ret;
}


// *************************************************************************************
// NegLogUnnormZDens

// call the function object

// 19/5/2011: generalise "binaryLogisticCorrection" to "higherOrderCorrection" and
//            also treat the canonical Poisson case.
double
NegLogUnnormZDens::operator()(double z)
{
    // map back to the original covariance factor scale
    double const g = exp(z);

    // optional status message
    if(verbose)
    {
        Rprintf("\nNegLogUnnormZDens: Computing function call NegLogUnnormZDens(%f) ...", z);
    }

    // if g is numerically zero, we cannot proceed
    if(g == 0.0)
    {
        std::ostringstream stream;
        stream << "g numerically zero";
        throw std::domain_error(stream.str().c_str());
    }

    // protect against errors in IWLS and non-finite result
    try
    {
        // compute Gaussian approximation:
        PosInt requiredIter = iwlsObject.startWithLastLinPred(nIter, // iterate at most nIter times (because we want convergence)
                                                              g);    // for this covariance factor, for the current linear predictor start.

        // if we did not converge in the maximum number of iterations, then start again from the original linear predictor,
        // and allow even more iterations.
        if(requiredIter > nIter)
        {
            PosInt higherIter = 2 * nIter;
            requiredIter = iwlsObject.startWithNewLinPred(higherIter,
                                                          g,
                                                          modelData.linPredStart);

            // if we still do not have converged, print a warning message if we are verbose.
            if(requiredIter > higherIter && verbose)
            {
                Rprintf("\nnNegLogUnnormZDens: IWLS did not converge in %d iterations (even after restart)!", higherIter);
            }

            // do not warn the user, because this need not be a serious problem.
        }

        if(verbose)
        {
            Rprintf("\nNegLogUnnormZDens: finished after %d IWLS iterations", requiredIter);
        }

        // get iwls results
        const IwlsResults iwlsResults = iwlsObject.getResults();

        // then the return value is:
        double ret = 0.5 * iwlsResults.logPrecisionDeterminant -
                     iwlsObject.getnCoefs() * M_LN_SQRT_2PI -    // This is dependent on the model dimension! (we had to add it here because
                                                            // computeLogUnPosteriorDens also contains now the analogous normalization constant.
                     iwlsObject.computeLogUnPosteriorDens(Parameter(iwlsResults.coefs, z));


        // the Raudenbush & Yang & Yosef correction to the Laplace approximation
        // in the canonical link case

        // Warning: this may not work yet for some reason for some data sets,
        // strange negative correctionFactor values can occur!
        // (both for the removed Fortran and for the actual C++ code...)
        if(higherOrderCorrection && modelData.canonicalLink)
        {

            // initialize E(T_4), E(T_6):
            // (excluding constant factors which are multiplied at the end)
            double exp_T4 = 0.0;
            double exp_T6 = 0.0;

            // for the sum over m_i^(3) * x_i * B_i:
            AVector k = arma::zeros(iwlsObject.getnCoefs());

            // design matrix is also needed
            const AMatrix& design = iwlsObject.getDesign();

            // iterate through the observations:
            for(PosInt i = 0; i < iwlsObject.nObs; ++i)
            {
                // calculate the fitted probability and resulting derivatives m_i^(x)
                // along the way.
                const double mu = modelData.link->linkinv(iwlsResults.linPred(i));
                double m3, m4, m6;

                // see page 147 in Raudenbush et al. for the formulas of m3, m4, m6
                if(modelData.familyString == "binomial")
                {
                    const double m2 = mu * (1.0 - mu);
                    m3 = m2 * (1.0 - 2.0 * mu);
                    m4 = m2 * (1.0 - 6.0 * m2);
                    m6 = m4 * (1.0 - 12.0 * m2) - 12.0 * m3 * m3;
                }
                else if (modelData.familyString == "poisson")
                {
                    m3 = mu;
                    m4 = mu;
                    m6 = mu;
                }
                else
                {
                    Rf_error("Higher order correction not implemented for this family.");
                }

                // add to the sums:
                AVector tmp = arma::trans(design.row(i));
                trs(false, false, iwlsResults.qFactor, tmp);
                const double B = arma::dot(tmp, tmp);
                const double B2 = B * B;

                exp_T4 += m4 * B2;
                exp_T6 += m6 * B2 * B;
                k += (m3 * B) * arma::trans(design.row(i));
            }

            // calculate k^T * Cov * k, and overwrite k with the intermediate result:
            trs(false, false, iwlsResults.qFactor, k);
            double exp_T3T3 = arma::dot(k, k);

            // So the correction factor for the conditional marginal likelihood is (with the new terms):
            double correctionFactor = - 1.0 / 8.0 * exp_T4 -
                    1.0 / 48.0 * exp_T6 +
                    5.0 / 24.0 * exp_T3T3;

            if(verbose)
            {
                Rprintf("\nNegLogUnnormZDens: Higher-order correction factor is %f", 1.0 + correctionFactor);
            }

            if(correctionFactor > -1.0)
            {
                // since we return here the negative log of the conditional marg lik:
                ret -= log1p(correctionFactor);
            } else {
                Rf_warning("negative value for correction factor! We are not using the higher order correction here.");
            }


        } // end if(higherOrderCorrection)

        // check finiteness of return value
        if(! R_finite(ret))
        {
            std::ostringstream stream;
            stream << "NegLogUnnormZDens() got non-finite result " << ret << " for z=" << z;
            throw std::domain_error(stream.str().c_str());
        }

        // and only if it is finite we return
        return ret;

    }
    // catch errors in IWLS and non-finite result
    catch (std::domain_error& error)
    {
        if(verbose)
        {
            Rprintf("For z=%f, the density value could not be computed because\n%s",
                    z, error.what());
        }
        Rf_warning("for z=%f, the density value could not be computed for the following model:\n%s\nCheck for near-collinearity of covariates.",
                   z, mod.print().c_str());

        // return NaN. This can be handled by the Brent optimize routine! It apparently replaces it (implicitely) by
        // the highest positive number.
        return R_NaN;
    }
}


// *************************************************************************************
// compute the log marginal likelihood of a specific model

// 10/02/2012: Also protect the creation of the NegLogUnnormZDens object
// against numerical problems.
double
glmGetLogMargLik(const ModelPar& modPar,
                 const GlmModelData& modelData,
                 const Computation& computation)
{
    // first the easy null model case
    if(modPar.isNullModel())
    {
        return modelData.nullModelLogMargLik;
    }
    else
    { // so now there is at least a linear term.

        // the return value will be placed in here:
        double ret = 0.0;

        // protect everything for problems in IWLS etc.
        try
        {
            // get negative log unnormalized z density: a function object.
            NegLogUnnormZDens negLogUnnormZDens(modPar,
                                                modelData,
                                                false,
                                                computation.debug,
                                                computation.higherOrderCorrection);

            // cache the function, because we do not want to evaluate it more often
            // than necessary.
            CachedFunction<NegLogUnnormZDens> cachedNegLogUnnormZDens(negLogUnnormZDens);

            // allocate space for zMode and zVar
            double zMode = 0.0;
            double zVar = 0.0;

            ret = logIntegral<CachedFunction<NegLogUnnormZDens> >(cachedNegLogUnnormZDens,
                                                                  zMode,
                                                                  zVar,
                                                                  computation);
        }
        catch (std::domain_error& error)
        {
            if(computation.debug)
            {
                Rprintf("\ngetLogMargLik: Model can not be included, because\n%s", error.what());
            }

            ret = R_NaN;
        }

        // return the estimate
        return ret;
    }
}

// *************************************************************************************
// ZdensApprox

// ctr
ZdensApprox::ZdensApprox(const ModelPar& modPar,
                         const GlmModelData& modelData,
                         const Computation& computation) :
                         zMode(0.0),
                         logMargLik(modelData.nullModelLogMargLik),
                         isNullModel(modPar.isNullModel()),
                         linApproxDens()
{
    if(! isNullModel)
    {
        // get negative log unnormalized z density: a function object.
        NegLogUnnormZDens negLogUnnormZDens(modPar,
                                            modelData,
                                            false,
                                            computation.debug,
                                            computation.higherOrderCorrection);

        // cache the function, because we do not want to evaluate it more often
        // than necessary.
        CachedFunction<NegLogUnnormZDens> cachedNegLogUnnormZDens(negLogUnnormZDens);

        // protect everything for problems
        try
        {
            // allocate space for zVar
            double zVar = 0.0;

            // compute log marginal likelihood, which also explores nicely
            // the unnormalised z density along the way:
            logMargLik = logIntegral<CachedFunction<NegLogUnnormZDens> >(cachedNegLogUnnormZDens,
                                                                         zMode,
                                                                         zVar,
                                                                         computation);

            // get the cached values
            Cache funCache = cachedNegLogUnnormZDens.getCache();

            AVector zVals = arma::conv_to<AVector>::from(funCache.getArgs());
            AVector logDensVals = arma::conv_to<AVector>::from(funCache.getVals());

            logDensVals = - logDensVals - logMargLik;

            // and then set up the approximation object
            linApproxDens = LinApproxDens(zVals,
                                          logDensVals);
        }
        catch (std::domain_error& error)
        {
            Rf_error("\nZdensApprox: Marginal z density approximation could not be constructed, because\n%s",
                     error.what());
        }
    }
}

// get the log density value at a specific z
double
ZdensApprox::logDens(double z) const
{
    if(isNullModel)
        return 0.0;
    else
    {
        return log(linApproxDens.dens(z));
    }
}

// get one sample from the approximate distribution
double
ZdensApprox::sample() const
{
    if(isNullModel)
        return 0.0;
    else
    {
        return linApproxDens.sample();
    }
}

// getter for the linear approximation
List
ZdensApprox::getLinApproxDens() const
{
    return List::create(_["args"] = linApproxDens.getArgs(),
                        _["dens"] = linApproxDens.getDens());
}
// *************************************************************************************
