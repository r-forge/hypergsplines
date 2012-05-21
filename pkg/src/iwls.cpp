#include <iwls.h>
#include <types.h>
#include <rcppExport.h>
// #include <cassert>
#include <stdexcept>
#include <sstream>
#include <linalgInterface.h>
#include <getRho.h>

using namespace Rcpp;

// criterion for comparison of two Column vectors of the same size
// max_j (abs(a_j - b_j) / abs(b_j) + 0.01)
// this is similar to the criterion used by R's glm routine on the deviance scale.
// However, we want to avoid the posterior scale because it would slow down the algorithm...
// (computing this criterion is easier than computing the likelihood * prior)
static double
criterion(const AVector& a, const AVector& b)
{
    // check lengths
    // assert(a.n_elem == b.n_elem);

    // this will be the returned value
    double ret = 0.0;

    // now iterate over the elements
#pragma omp parallel for
    for (PosInt j = 0; j < a.n_elem; ++j)
    {
        double tmp = fabs(a(j) - b(j)) / (fabs(b(j)) + 0.01);

        // note that the "critical" directive is vital here!!
        // Otherwise two executions of the same code can lead to different answers.
#pragma omp critical
        ret = (tmp > ret) ? tmp : ret; /* fmax(ret, tmp); */
    }

    // return the criterion value.
    return ret;
}


// This function checked 7/4/2011:

// constructor: constructs the Iwls object for given model and data.
Iwls::Iwls(const ModelPar& mod,
           const GlmModelData& modelData,
           bool useFixedZ,
           double epsilon,
           bool verbose) :
           isNullModel(mod.isNullModel()),
           useFixedZ(useFixedZ),
           nObs(modelData.nObs),
           modelData(modelData),
           epsilon(epsilon),
           verbose(verbose),
           invSqrtDispersions(1.0 / arma::sqrt(modelData.dispersions)),
           dimLinear(0),
           dimSpline(0),
           dimZ(0),
           logUnscaledPriorPrecDeterminant(0.0)
{
    // translate config into linear and splines part:
    IntVector whichLinear;
    IntVector whichSpline;

    for(int i_cov = 0; i_cov < modelData.nCovs; ++i_cov)
    {
        const double d = mod.config[i_cov];

        if(d > 0)
        {
            whichLinear.push_back(i_cov);

            if(d > 1)
            {
                whichSpline.push_back(i_cov);
            }
        }
    }

    // corresponding number of linear and spline effects
    dimLinear = whichLinear.size();
    dimSpline = whichSpline.size();
    dimZ = dimSpline * modelData.dimSplineBasis;

    // so the number of coefficients altogether is:
    nCoefs = 1 + dimLinear + dimZ;

    // now that we finally know what nCoefs is, set up the iwls results container:
    results = IwlsResults(modelData.linPredStart, nCoefs);

    // and we can initialize the unscaled prior precision matrix
    // and design matrix accordingly:
    unscaledPriorPrec.zeros(nCoefs, nCoefs);

    design.set_size(nObs, nCoefs);
    design.col(0) = arma::ones(nObs);

    // both commands above have now already covered the intercept portion
    // (zero prior precision and vector of ones in the design matrix)

    // if there is at least one linear effect:
    if(! isNullModel)
    {
        // construct the linear design matrix Xlin (without the intercept):
        arma::mat Xlin(nObs, dimLinear);
        double * Xlin_ptr = Xlin.memptr();

        for(int i_lin = 0; i_lin < dimLinear; ++i_lin)
        {
            const int start = whichLinear.at(i_lin) * nObs;

            std::copy(modelData.Xfull_ptr + start, // start source pointer
                      modelData.Xfull_ptr + start + nObs, // one beyond the end source pointer
                      Xlin_ptr + (i_lin * nObs)); // start destination pointer
        }

        // input the linear part into the design matrix
        design.cols(1, dimLinear) = Xlin;

        // scale the linear part with the sqrt of W_0:
        Xlin = arma::diagmat(modelData.sqrtWeightMatrixEntries) * Xlin;

        // and compute XtWX = Xlin^T * weightMatrix * Xlin:
        // (we use the storage which will later be occupied by the Cholesky
        // factor)
        linPartCholFactor.set_size(dimLinear, dimLinear);
        syrk(false, // so the result is in lower triangular storage!
             true,
             1.0,
             Xlin,
             0.0,
             linPartCholFactor);
        // todo: similarly to the assembly of M below, we could do it here
        // for XtWX, that could save a bit time (but probably not much since
        // the dimension of X is mostly small)

        // if there is at least one spline term:
        if(dimSpline > 0)
        {
            // this is at first very similar to the normal distribution code,
            // because we also needed to construct the Fisher information there.
            // The new thing first here is that the weight matrix W_0 comes into play.

            // first we compute
            // A = Z M^-1 Z^T,

            // where Z is the grand spline basis matrix of dimension
            // n x pK, when p is the number of spline covariates and
            // K is the number of knots,

            // and M = Z^T * W_0 * Z + diag(rhos)^-1 is of dimension pK x pK,
            // and thus can be Cholesky factorized independent of the
            // number of observations.

            // first assemble M and Z^T in parallel,
            arma::mat M(dimZ, dimZ);
            arma::mat Zt(dimZ, nObs);

            // also the spline coef part of the unscaled precision matrix can be filled in
            splinePartPrecSqrtEntries.set_size(dimZ);

            // and the spline portion of the determinant can be computed in parallel:
            double logSplineDet = 0.0;

            // process all spline terms
            for(int i = 0; i < dimSpline; ++i)
            {
                // get the covariate index:
                const int s = whichSpline.at(i);

                // input transpose of this Z into Zt:
                NumericMatrix thisZ_temp = modelData.zList[s];
                const arma::mat thisZ(thisZ_temp.begin(), nObs, modelData.dimSplineBasis, false);

                // syntax:
                // submat(first row, first col, last row, last col)
                Zt.submat(i * modelData.dimSplineBasis,
                          0,
                          (i + 1) * modelData.dimSplineBasis - 1,
                          nObs - 1) = arma::trans(thisZ);

                // start with filling this "line" no. i in M:
                // (Note that we only fill the *upper*-triangular part)
                for(int j = i; j < dimSpline; ++j)
                {
                    // get the covariate index:
                    const int t = whichSpline.at(j);

                    // then copy Z_s^T * W_0 * Z_t into M:
                    M.submat(i * modelData.dimSplineBasis,
                             j * modelData.dimSplineBasis,
                             (i + 1) * modelData.dimSplineBasis - 1,
                             (j + 1) * modelData.dimSplineBasis - 1) = modelData.ZtZarray.at(s).at(t - s);
                }

                // get the correct rho:
                const double d = mod.config[s];
                double rho;

                if(d == floor(d))
                {
                    // dNumber <- match(d - 1,
                    //                  modelData$splineDegrees)
                    const int dNumber = round(d - 2);

                    const NumericVector rhos = modelData.rhoList[s];
                    rho = rhos[dNumber];
                }
                else
                {
                    const NumericVector lambdas = modelData.lambdasList[s];
                    rho = getRho(lambdas, d - 1.0);
                }

                // add the inverse of this rho to the diagonal block
                // of M:
                // (and in parallel also fill it into unscaledPriorPrec)
                for(int k = i * modelData.dimSplineBasis;
                        k < (i + 1) * modelData.dimSplineBasis;
                        ++k)
                {
                    double rhoInv = 1.0 / rho;

                    M.at(k, k) += rhoInv;
                    splinePartPrecSqrtEntries.at(k) = sqrt(rhoInv);

                    int l = k + 1 + dimLinear;
                    unscaledPriorPrec.at(l, l) = rhoInv;
                }

                // update the log determinant:
                logSplineDet += log(rho);
            }
            logSplineDet *= (- 1.0 * modelData.dimSplineBasis);
            // now the spline part of the log determinant is finished,
            // and can be copied into the storage:
            logUnscaledPriorPrecDeterminant = logSplineDet;

            // add the spline part to the design matrix:
            design.cols(1 + dimLinear, nCoefs - 1) = arma::trans(Zt);

            // OK, now we replace M by its upper triangular Cholesky factor:
            potrf(true,
                  M);
            // so now A = Z (M^T M)^-1 Z^T
            //          = Z M^-1 (M^-1)^T Z^T
            //          = W W^T
            // where W^T = (M^-1)^T Z^T = (M^T)^-1 Z^T,
            // or M^T W^T = Z^T.
            // Since we only need W^T in what follows, we can use the storage in Zt for it.

            trs(true,
                true,
                M,
                Zt);
            // so now we have A = Zt^T Zt !

            // The Fisher information is
            // X^T W_0 X - X^T W_0 Zt^T Zt W_0 X.
            // We already have the first part (X^T W_0 X) in "linPartCholFactor",
            // so we only need to update that with the second part,
            // X^T W_0 Zt^T Zt W_0 X == crossprod(ZtWX):
            Zt = Zt * arma::diagmat(modelData.sqrtWeightMatrixEntries);
            AMatrix ZtWX = Zt * Xlin;
            syrk(false, // linPartCholFactor is in lower triangular storage!
                 true,
                 - 1.0,
                 ZtWX,
                 1.0,
                 linPartCholFactor);
            // now linPartCholFactor contains the part of the unscaled prior precision matrix
            // corresponding to the linear (or "fixed") effects.
        }

        // now "linPartCholFactor" is complete in both cases, and we can fill it into
        // the grand precision matrix:
        unscaledPriorPrec.submat(1, 1, dimLinear, dimLinear) = linPartCholFactor;
        // [ actually for this copying we would not need to copy the upper triangle
        // of "linPartCholFactor" (because it is in lower triangular storage), but ok for now. ]

        // use the cholesky routine in order to get the log determinant
        potrf(false,
              linPartCholFactor);
        // now "linPartCholFactor" finally contains the lower-triangular cholesky factor
        // of the linear part of the prior precision matrix

        // so we can update the log determinant
        logUnscaledPriorPrecDeterminant +=
                2.0 * arma::as_scalar(arma::sum(arma::log(arma::diagvec(linPartCholFactor))));

    }
}




// do the Iwls algorithm for a given covariance factor g and start linear predictor linPred,
// until convergence or until the maximum number of iterations is reached.
// Note that the linear predictor is the sum of X^T * beta and the vector of offsets.
// so also only one iwls step can be performed with this function.
// returns the number of iterations.

// 19/5/2011 modified to accomodate offsets
PosInt
Iwls::startWithLastLinPred(PosInt maxIter,
                           double g)
{
    // initialize iteration counter and stopping criterion
    PosInt iter = 0;
    bool converged = false;

    // do IWLS for at most 30 iterations and unfulfilled stopping criterion
    while ((iter++ < maxIter) && (! converged))
    {
        // compute the pseudo-observations and corresponding sqrt(weights) from the linear predictor
        AVector pseudoObs(nObs);
        AVector sqrtWeights(invSqrtDispersions);

#pragma omp parallel for
        for(PosInt i = 0; i < nObs; ++i)
        {
            double mu = modelData.link->linkinv(results.linPred(i));
            double dmudEta = modelData.link->mu_eta(results.linPred(i));

            pseudoObs(i) = results.linPred(i) - modelData.offsets(i) +
                           (modelData.y(i) - mu) / dmudEta;
            sqrtWeights(i) *= dmudEta / sqrt(modelData.distribution->variance(mu));
        }

        // calculate sqrt(W)X, which is needed twice
        AMatrix sqrtWX = arma::diagmat(sqrtWeights) * design;

        // now calculate the precision matrix Q by doing a rank update:
        // Q = crossprod(sqrt(W)X) + scaled PriorPrec
        results.qFactor = unscaledPriorPrec;

        if(! isNullModel)
        {
            // scale the linear part of the prior precision with 1/g
            results.qFactor.submat(1, 1, dimLinear, dimLinear) /= g;
        }

        // then do the "crossprod" update:
        syrk(false,
             true,
             1.0,
             sqrtWX,
             1.0,
             results.qFactor);

        // decompose into Cholesky factor, Q = LL':
        potrf(false,
              results.qFactor);

        // save the old coefficients vector
        AVector coefs_old = results.coefs;

        // the rhs of the equation Q * m = rhs   or    R'R * m = rhs
        pseudoObs = arma::diagmat(sqrtWeights) * pseudoObs;
        results.coefs = arma::trans(sqrtWX) * pseudoObs;

        // forward-backward solve LL' * v = rhs
        potrs(false,
              results.qFactor,
              results.coefs);
        // now results.coefs is finished.

        // so the new linear predictor is
        results.linPred = design * results.coefs + modelData.offsets;

        // compare on the coefficients scale, but not in the first iteration where
        // it is not clear from where coefs_old came. Be safe and always
        // decide for non-convergence in this case.
        converged = (iter > 1) ? (criterion(coefs_old, results.coefs) < epsilon) : false;
    }

    // do not (!)
    // warn if IWLS did not converge within the maximum number of iterations
    // because the maximum number of iterations can be set by user of this function.

    // compute log precision determinant
    results.logPrecisionDeterminant = 2.0 * arma::as_scalar(arma::sum(arma::log(arma::diagvec(results.qFactor))));

    // last but not least return the number of iterations
    // Note that we must subtract 1 here to give the correct number!
    return iter;
}

// do the Iwls algorithm for a given covariance factor g and new start linear predictor
// linPredStart.
PosInt
Iwls::startWithNewLinPred(PosInt maxIter,
                          double g,
                          const AVector& linPredStart)
{
    // copy start value into linear predictor of the iwls object
    results.linPred = linPredStart;

    // then run the iwls algorithm
    return startWithLastLinPred(maxIter, g);
}


// do the Iwls algorithm for a given covariance factor g and new start coefficients vector
// coefsStart.

// 19/5/2011 modified to accomodate offsets
PosInt
Iwls::startWithNewCoefs(PosInt maxIter,
                  double g,
                  const AVector& coefsStart)
{
    // start with new linpred deriving from the coefs
    return startWithNewLinPred(maxIter, g, design * coefsStart + modelData.offsets);
}




// compute the log of the (unnormalized)
// posterior density for a given parameter consisting of the coefficients vector and z.
//
// Note that it is important to incorporate all model-depending constants here,
// because this function is also used in the Chib-Jeliazkov marginal likelihood estimation,
// comparing different models!!
//
// useFixedZ: is the log-covariance factor z fixed? Then the conditional posterior
// density of the coefficients vector is returned (so the prior of z is not included).

// 7/4/2011 checked the function
// 19/5/2011 modified to accomodate offsets
// 20/4/2012 other computation of scaledSplineCoefs to avoid Armadillo problems
double
Iwls::computeLogUnPosteriorDens(const Parameter& sample) const
{
    // compute the sample of the linear predictor:
    AVector linPredSample = design * sample.coefs + modelData.offsets;

    // compute the resulting mean vector from the linear predictor via the response function
    AVector meansSample(linPredSample.n_elem);

#pragma omp parallel for
    for(PosInt i = 0; i < meansSample.n_elem; ++i)
    {
        meansSample(i) = modelData.link->linkinv(linPredSample(i));
    }

    // start with the log likelihood of this coefficients, it is always included.
    // this part is included in both cases because it does not depend on
    // the prior on the (non-intercept) effects:

    double ret = modelData.distribution->loglik(meansSample.memptr());

    // now it depends again on null model or not.

    if(! isNullModel)
    {
        // map z sample on the original g scale
        double g = exp(sample.z);

        // calculate the quadratic forms:

        // first for the linear coefficients
        AVector scaledLinCoefs = sample.coefs.subvec(1, dimLinear);
        trmv(false,
             true,
             linPartCholFactor,
             scaledLinCoefs);
        double linearQuadForm = arma::dot(scaledLinCoefs, scaledLinCoefs);

        // then (optionally) for the spline coefficients
        double splineQuadForm = 0.0;
        if(dimZ > 0)
        {
            AVector scaledSplineCoefs = splinePartPrecSqrtEntries % sample.coefs.subvec(dimLinear + 1, nCoefs - 1);
            splineQuadForm = arma::dot(scaledSplineCoefs, scaledSplineCoefs);
        }

        // now add the non-null model specific part, which comes from the prior on
        // the coefficients
        ret += 0.5 * (logUnscaledPriorPrecDeterminant -
                     (dimLinear + dimZ) * 2.0 * M_LN_SQRT_2PI -
                     dimLinear * sample.z -
                     linearQuadForm / g -
                     splineQuadForm);

        if(! useFixedZ)
        {
            // and the log prior of this g
            double logGPrior = modelData.gPrior->logDens(g);

            // add the log prior of z
            ret += logGPrior + sample.z;

            // note that the sample.z has its origin in the density transformation from g to z.
            // if the prior on g was discrete and not continuous, this part would not be included in general.
        }
    }

    // return the correct value
    return ret;
}
