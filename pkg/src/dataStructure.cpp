/*
 * dataStructure.cpp
 *
 *  Created on: 14.03.2011
 *      Author: daniel
 *
 *  30/03/2012: remove hyperparameter a, instead use string for g-prior
 */

#include <dataStructure.h>
#include <string>
#include <sstream>
#include <sum.h>
#include <getRho.h>
#include <logBFhypergn.h>
#include <logBFhyperg.h>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


using namespace Rcpp;

// *************************************************************************************

// GaussHermite: encapsulate the Gauss-Hermite integration stuff

// compute nodes and log weights for given mode and sd of target unnormalized density
void
GaussHermite::getNodesAndLogWeights(double mode,
                                    double var,
                                    MyDoubleVector& nodes,
                                    MyDoubleVector& logWeights) const
{
    // logarithm of square root of (2 * var).
    double logSqrt2Var = 0.5 * (M_LN2 + log(var));

    MyDoubleVector::const_iterator t = tVec.begin();
    for (MyDoubleVector::const_iterator w = wVec.begin(); w != wVec.end(); ++w, ++t)
    {
        nodes.push_back(mode + exp(logSqrt2Var) * (* t));
        logWeights.push_back(log(* w) + (* t) * (* t) + logSqrt2Var);
    }
}


// *************************************************************************************

// ModelData:

// ctr from R list
ModelData::ModelData(SEXP R_modelData) :
    modelData(R_modelData),
    rhoList(as<List>(modelData["rho.list"])),
    zList(as<List>(modelData["Z.list"])),
    lambdasList(as<List>(modelData["lambdas.list"])),
    Xfull(as<NumericMatrix>(modelData["X"])),
    Xfull_ptr(Xfull.begin()),
    nObs(as<int>(modelData["nObs"])),
    nCovs(as<int>(modelData["nCovs"])),
    dimSplineBasis(as<int>(modelData["dimSplineBasis"])),
    continuous(as<LogicalVector>(modelData["continuous"])),
    gPriorString(as<std::string>(modelData["gPrior"])),
    degrees(as<IntVector>(modelData["degrees"])),
    nDegrees(degrees.size()),
    y(as<NumericVector>(modelData["y"])),
    // meanY(arma::mean(y) * arma::ones<arma::colvec>(nObs, 1)),
    y_centered(y - arma::mean(y)),
    yCenterNormSq(arma::dot(y_centered, y_centered))
{
    // indices of continuous covariates:
    for(int i = 0; i < nCovs; ++i)
    {
        if(continuous[i])
            contCovs.push_back(i);
    }
    nContCovs = contCovs.size();

    // check the degrees vector:
    if(degrees[0] != 0 || degrees[1] != 1)
    {
        Rf_error("Degree vector must start with 0, 1 !");
    }

    // compute Z_i^T Z_j for each pair 0 <= i <= j < nCovs
    // (so we are interested in the *upper* triangle of the grand matrix Z^T Z)
    for(int i = 0; i < nCovs; ++i)
    {
        MatVector line_i;

        NumericMatrix Z_i_temp = zList[i];
        const arma::mat Z_i(Z_i_temp.begin(), nObs, dimSplineBasis, false);

        for(int j = i; j < nCovs; ++j)
        {
            NumericMatrix Z_j_temp = zList[j];
            const arma::mat Z_j(Z_j_temp.begin(), nObs, dimSplineBasis, false);

            line_i.push_back(arma::trans(Z_i) * Z_j);
        }

        ZtZarray.push_back(line_i);
        // so ZtZarray.at(i).at(j - i) equals Z_i^T * Z_j
    }
}



// 15/6/2012: correct criterion for non-identifiability
// compute the log marginal likelihood (and as a byproduct the R2)
// of a specific model
double
ModelData::getLogMargLik(const ModelPar& modPar,
                         double& R2)
{
    const MyDoubleVector& config = modPar.config;

    // translate config into linear and splines part
    std::vector<int> whichLinear;
    std::vector<int> whichSpline;

    for(int i_cov = 0; i_cov < nCovs; ++i_cov)
    {
        const double d = config[i_cov];

        if(d > 0.0)
        {
            whichLinear.push_back(i_cov);

            if(d > 1.0)
            {
                whichSpline.push_back(i_cov);
            }
        }
    }

    // is this the null model / only linear effects?
    const int dimLinear = whichLinear.size();
    const bool isNullModel = (dimLinear == 0);

    const int dimSpline = whichSpline.size();
    const bool hasOnlyLinear = (dimSpline == 0);

    // return NaN if the number of covariates + intercept is not smaller than the number of observations
    // (because the arma::solve below could still work in that case...)
    // Note that in the border case of dimLinear + 1 = nObs we have R2 = 1 and then
    // the log marginal likelihood is not finite.
    if(dimLinear + 1 >= nObs)
    {
        // if(verbose)
        // {
        //      Rprintf("getLogMargLik: p > n for model:\n%s",
        //              modPar.print().c_str());
        // }

        return R_NaN;
    }

    // first the easy null model case
    if(isNullModel)
    {
        R2 = 0.0;
        return - (nObs - 1.0) / 2.0 * log(yCenterNormSq);
    }
    else
    { // so now there is at least a linear term.

        // construct the linear design matrix Xlin (without the intercept):
        arma::mat Xlin(nObs, dimLinear);
        double * Xlin_ptr = Xlin.memptr();

        for(int i_lin = 0; i_lin < dimLinear; ++i_lin)
        {
            const int start = whichLinear.at(i_lin) * nObs;

            std::copy(Xfull_ptr + start, // start source pointer
                      Xfull_ptr + start + nObs, // one beyond the end source pointer
                      Xlin_ptr + (i_lin * nObs)); // start destination pointer
        }

        if(hasOnlyLinear)
        {
            // protect against errors in the betaOLS computation, coming from
            // collinear columns in Xlin, e.g.
            arma::colvec betaOLS;
            try
            {
                // compute the OLS solution and get the coefficient of determination (R^2):
                betaOLS = arma::solve(Xlin, y);
            }
            catch (std::runtime_error& error)
            {
//                if(verbose)
//                {
//                    Rprintf("getLogMargLik: betaOLS could not be computed for model:\n%s",
//                            modPar.print().c_str());
//                }

                return R_NaN;
            }
            arma::colvec y_fitted = Xlin * betaOLS;
            R2 = arma::dot(y_fitted, y_fitted) / yCenterNormSq;

            // R2 must not be larger than 1
            R2 = std::min(R2, 1.0);

            // compute the log marginal likelihood:

            // start with the log marginal likelihood of the null model
            double ret = - (nObs - 1.0) / 2.0 * log(yCenterNormSq);

            // and then add the log BF, depending on the hyperprior on g:
            if(gPriorString == "hyper-g/n")
            {
                ret += logBFhypergn(nObs, dimLinear, R2);
            }
            else // gPriorString == "hyper-g"
            {
                ret += logBFhyperg(nObs, dimLinear, R2);
            }

            // return the result
            return ret;
        }
        else // the most interesting case with spline-modelled covariate effects!
        {
            // we now construct (components of) the *precision* (or weight) matrix
            // of the marginal model.

            // this is Vinv = I - Z M^-1 Z^T,

            // where Z is the grand spline basis matrix of dimension
            // n x pK, when p is the number of spline covariates and
            // K is the number of knots.

            // M = Z^T Z + diag(rhos)^-1 is of dimension pK x pK,
            // and thus can be Cholesky factorized independent of the
            // number of observations.

            // first assemble M and Z^T in parallel:
            const int dimZ = dimSpline * dimSplineBasis;
            arma::mat M(dimZ, dimZ);
            arma::mat Zt(dimZ, nObs);

            for(int i = 0; i < dimSpline; ++i)
            {
                // get the covariate index:
                const int s = whichSpline.at(i);

                // input transpose of this Z into Zt:
                NumericMatrix thisZ_temp = zList[s];
                const arma::mat thisZ(thisZ_temp.begin(), nObs, dimSplineBasis, false);

                // syntax:
                // submat(first row, first col, last row, last col)
                Zt.submat(i * dimSplineBasis,
                          0,
                          (i + 1) * dimSplineBasis - 1,
                          nObs - 1) = arma::trans(thisZ);

                // start with filling this "line" no. i in M:
                for(int j = i; j < dimSpline; ++j)
                {
                    // get the covariate index:
                    const int t = whichSpline.at(j);

                    // then copy Z_s^T Z_t into M:
                    M.submat(i * dimSplineBasis,
                             j * dimSplineBasis,
                             (i + 1) * dimSplineBasis - 1,
                             (j + 1) * dimSplineBasis - 1) = ZtZarray.at(s).at(t - s);
                }

                // get the correct rho
                const double d = config[s];
                double rho;

                if(d == floor(d))
                {
                    // dNumber <- match(d - 1,
                    //                  modelData$splineDegrees)
                    const int dNumber = round(d - 2);

                    const NumericVector rhos = rhoList[s];
                    rho = rhos[dNumber];
                }
                else
                {
                    const NumericVector lambdas = lambdasList[s];
                    rho = getRho(lambdas, d - 1.0);
                }

                // add the inverse of this rho to the diagonal block
                // of M:
                for(int k = i * dimSplineBasis; k < (i + 1) * dimSplineBasis; ++k)
                {
                    M.at(k, k) += 1.0 / rho;
                }
            }

            // OK, now we replace M by its upper triangular Cholesky factor:
            int info = 0;
            F77_CALL(dpotrf)("U",
                    &dimZ,
                    M.memptr(),
                    &dimZ,
                    &info);
            if(info != 0)
            {
//                    Rf_error("Error in Cholesky factorization of M, code %d", info);
                return R_NaN;
            }

            // so now Vinv = I - Z (M^T M)^-1 Z^T
            //             = I - Z M^-1 (M^-1)^T Z^T
            //             = I - W W^T
            // where W^T = (M^-1)^T Z^T = (M^T)^-1 Z^T,
            // or M^T W^T = Z^T.
            // Since we only need W^T in what follows, we can use the storage in Zt for it.

            // we use here the level 3 blas routine
            // DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            const double doubleOne = 1.0;
            F77_CALL(dtrsm)("L", // side : left
                    "U", // uplo : M is upper-triangular
                    "T", // transa: use M^T in the equation
                    "N", // diag : not unit triangular
                    &dimZ, // M: number of rows of Zt
                    &nObs, // N: number of columns of Zt
                    &doubleOne, // alpha: do not scale Zt
                    M.memptr(), // the Cholesky factor M
                    &dimZ, // number of rows of M
                    Zt.memptr(), // the rhs
                    &dimZ); // number of rows of Zt

            // so now we have Vinv = I - Zt^T Zt !

            // we can compute the sum of squares total
            arma::colvec Zty_centered = Zt * y_centered;
            const double sst = yCenterNormSq - arma::dot(Zty_centered, Zty_centered);

            // for the sum of squares model:
            arma::mat left = trans(Xlin) - trans(Zt * Xlin) * Zt;
            arma::colvec rhs = left * y;
            left *= Xlin;

            info = 0;
            F77_CALL(dpotrf)("U",
                    &dimLinear,
                    left.memptr(),
                    &dimLinear,
                    &info);
            if(info != 0)
            {
//                    Rf_error("Error in Cholesky factorization of left, code %d", info);
                return R_NaN;
            }
            // now left contains the Cholesky factor R of X^T V^-1 X

            // now compute v = R * betaOLS as the solution of R^T v = rhs
            // and overwrite rhs with that:
            const int intOne = 1;
            F77_CALL(dtrsm)("L", // side : left
                    "U", // uplo : M is upper-triangular
                    "T", // transa: use M^T in the equation
                    "N", // diag : not unit triangular
                    &dimLinear, // M: number of rows of rhs
                    &intOne, // N: number of columns of rhs
                    &doubleOne, // alpha: do not scale rhs
                    left.memptr(), // the Cholesky factor
                    &dimLinear, // number of rows of the factor
                    rhs.memptr(), // the rhs
                    &dimLinear); // number of rows of rhs

            // so we have the ssm:
            const double ssm = arma::dot(rhs, rhs);

            // this is more efficient than:
            // const arma::colvec Xt_Vinv_y = Xt_Vinv * y;
            // const arma::colvec betaOLS = arma::solve(Xt_Vinv * Xlin,
            //                                          Xt_Vinv_y);
            // const double ssm = arma::dot(betaOLS, Xt_Vinv_y);

            // so we have the R2:
            R2 = ssm / sst;

            // R2 must not be larger than 1
            R2 = std::min(R2, 1.0);

            // for the marginal likelihood, we need the determinant of Vinv.
            // luckily, we can compute that as
            // the determinant of I - Zt Zt^T. Use the storage of M for this.
            M = arma::eye(dimZ, dimZ) - Zt * trans(Zt);
            double Vinv_log_det;
            double Vinv_log_det_sign;
            arma::log_det(Vinv_log_det,
                          Vinv_log_det_sign,
                          M);
            if(Vinv_log_det_sign < 0)
            {
//                    Rf_error("Determinant not positive");
                return R_NaN;
            }

            // now we can compute the log marginal likelihood:

            // start with the log marginal likelihood of the null model,
            // and the Jacobian
            double ret = - (nObs - 1.0) / 2.0 * log(sst) + 0.5 * Vinv_log_det;

            // and then add the log BF, depending on the hyperprior on g:
            if(gPriorString == "hyper-g/n")
            {
                ret += logBFhypergn(nObs, dimLinear, R2);
            }
            else // gPriorString == "hyper-g"
            {
                ret += logBFhyperg(nObs, dimLinear, R2);
            }

            // return the result
            return ret;
        }
    }
}


// *************************************************************************************

// GlmModelData

// constructor

// 19/5/2011: names of distribution and link are saved in the object
// 20/9/2011: move to g prior classes in R_modelData object, so do no longer save
//            the gPriorString
GlmModelData::GlmModelData(SEXP R_modelData) :
    rcpp_modelData(R_modelData),
    rhoList(as<List>(rcpp_modelData["rho.list"])),
    zList(as<List>(rcpp_modelData["Z.list"])),
    lambdasList(as<List>(rcpp_modelData["lambdas.list"])),
    familyList(as<List>(rcpp_modelData["family"])),
    Xfull(as<NumericMatrix> (rcpp_modelData["X"])),
    Xfull_ptr(Xfull.begin()),
    nObs(as<int> (rcpp_modelData["nObs"])),
    nCovs(as<int> (rcpp_modelData["nCovs"])),
    dimSplineBasis(as<int> (rcpp_modelData["dimSplineBasis"])),
    continuous(as<LogicalVector> (rcpp_modelData["continuous"])),
    degrees(as<IntVector> (rcpp_modelData["degrees"])),
    nDegrees(degrees.size()),
    y(as<NumericVector> (rcpp_modelData["Y"])),
    y_centered(y - arma::mean(y)),
    yCenterNormSq(arma::dot(y_centered, y_centered)),
    weightMatrixEntries(as<NumericVector> (rcpp_modelData["weightMatrixEntries"])),
    sqrtWeightMatrixEntries(arma::sqrt(weightMatrixEntries)),
    linPredStart(as<NumericVector>(familyList["linPredStart"])),
    dispersions(as<NumericVector>(familyList["dispersions"])),
    offsets(as<NumericVector>(familyList["offsets"])),
    nullModelLogMargLik(as<double> (rcpp_modelData["nullModelLogMargLik"])),
    familyString(as<std::string>(familyList["family"])),
    linkString(as<std::string>(familyList["link"])),
    canonicalLink((familyString == "binomial" && linkString == "logit") ||
                  (familyString == "poisson" && linkString == "log"))
{
    // indices of continuous covariates:
    for (int i = 0; i < nCovs; ++i)
    {
        if (continuous[i])
            contCovs.push_back(i);
    }
    nContCovs = contCovs.size();

    // check the degrees vector:
    if (degrees[0] != 0 || degrees[1] != 1)
    {
        Rf_error("Degree vector must start with 0, 1 !");
    }

    // todo: can we accelerate this computation??
    // (Better linear algebra, OpenMP, ...)
    // Perhaps it is cheaper to really compute Z^T * W_0 * Z as a whole
    // as a diagonal matrix, and then select the required subblocks appropriately
    // in iwls.cpp
    // One idea to test this: Just have exactly this code in a function
    // to be called from R, and a switch which selects the approach.

    // compute Z_i^T * W_0 * Z_j for each pair 0 <= i <= j < nCovs
    // (so we are interested in the *upper* triangle of the grand matrix Z^T * W_0 * Z)
    for (int i = 0; i < nCovs; ++i)
    {
        MatVector line_i;

        NumericMatrix Z_i_temp = zList[i];
        const arma::mat Z_i(Z_i_temp.begin(),
                            nObs,
                            dimSplineBasis,
                            false);

        for (int j = i; j < nCovs; ++j)
        {
            NumericMatrix Z_j_temp = zList[j];
            const arma::mat Z_j(Z_j_temp.begin(),
                                nObs,
                                dimSplineBasis,
                                false);

            line_i.push_back(arma::trans(Z_i)
                    * arma::diagmat(weightMatrixEntries) * Z_j);
        }

        ZtZarray.push_back(line_i);
        // so ZtZarray.at(i).at(j - i) equals Z_i^T * W_0 * Z_j
    }

    // first get the class name of the S4 g-prior object
    Rcpp::S4 rcpp_gPrior = as<S4>(rcpp_modelData["gPrior"]);
    std::string gPriorString = rcpp_gPrior.attr("class");

    // get the g prior object
    if(gPriorString == "HypergnPrior")
    {
        gPrior = new HypergnPrior(as<double>(rcpp_gPrior.slot("a")),
                                  as<int>(rcpp_gPrior.slot("n")));
    }
    else if (gPriorString == "HypergPrior")
    {
        gPrior = new HypergPrior(as<double>(rcpp_gPrior.slot("a")));
    }
    else if (gPriorString == "InvGammaGPrior")
    {
        gPrior = new InvGammaGPrior(as<double>(rcpp_gPrior.slot("a")),
                                    as<double>(rcpp_gPrior.slot("b")));
    }
    else if (gPriorString == "CustomGPrior")
    {
        gPrior = new CustomGPrior(as<SEXP>(rcpp_gPrior.slot("logDens")));
    }
    else
    {
        Rf_error("g-prior not implemented!");
    }

    // and now to the family business:

    // generate the distribution:

    const AVector weights = as<NumericVector> (familyList["weights"]);
    const double phi = familyList["phi"];

    if (familyString == "binomial")
    {
        distribution = new Binomial(y,
                                    weights);
    }
    else if (familyString == "gaussian")
    {
        distribution = new Gaussian(y,
                                    weights,
                                    phi);
    }
    else if (familyString == "poisson")
    {
        distribution = new Poisson(y,
                                   weights);
    }
    else
    {
        Rf_error("Distribution not implemented");
    }

    // generate the link:

    if (linkString == "logit")
    {
        link = new LogitLink();
    }
    else if (linkString == "probit")
    {
        link = new ProbitLink();
    }
    else if (linkString == "cloglog")
    {
        link = new CloglogLink();
    }
    else if (linkString == "inverse")
    {
        link = new InverseLink();
    }
    else if (linkString == "log")
    {
        link = new LogLink();
    }
    else if (linkString == "identity")
    {
        link = new IdentityLink();
    }
    else
    {
        Rf_error("Link not implemented!");
    }

}

// *************************************************************************************

// ModelPar

// compute degree index vector, given a degree vector
// note that if one element in "config" is not in "degrees",
// then the resulting degree index will be the length of
// "degrees", so one past the last element of "degrees"
void
ModelPar::compDegIndex(const IntVector& degrees)
{
    degIndex.clear();

    for(MyDoubleVector::const_iterator
            i = config.begin();
            i != config.end();
            ++i)
    {
        PosInt t = 0;
        while((t < degrees.size()) & (static_cast<double>(degrees[t]) != *i))
        {
            ++t;
        }
        degIndex.push_back(t);
    }
}

// return a textual description of this model configuration
std::string
ModelPar::print() const
{
    // put everything into a stringstream, and return the string conversion at the end.
    std::ostringstream stream;

    // Start
    stream << "\nmodel with the following degrees of freedom:\n";

    // Iterate over elements
    for(MyDoubleVector::const_iterator i = config.begin(); i != config.end(); ++i)
    {
        stream << *i << " ";
    }

    // End
    stream << "\n";
    return stream.str();
}

// *************************************************************************************

// StatusBar

// 19/5/2011: modified to accomodate small totaliters (< nBars) and really use numberGap!
StatusBar::StatusBar(bool p,
                     int nB,
                     int nG,
                     int ti) :
                     print(p),
                     nBars(std::min(nB, ti)),
                     numberGap(nG),
                     totaliters(ti),
                     barInterval(totaliters / nBars)
{
    if(print)
    {
        IntegerVector steps = seq_len(nBars + 1);
        steps = floor((steps - 1.0) * 100.0 / nBars);

        // first line with the numbers
        bool lastWasTwoDigits = false;
        for(int k = 1; k <= nBars; ++k)
        {
            if(k % numberGap == 0)
            {
                Rprintf("%d", steps[k]);

                lastWasTwoDigits = steps[k] >= 10;
            }
            else
            {
                if(lastWasTwoDigits)
                    lastWasTwoDigits = false;
                else
                    Rprintf("-");
            }
        }
        Rprintf("\n");

        // second line with the ticks
        for(int k = 1; k <= nBars; ++k)
        {
            if(k % numberGap == 0)
                Rprintf("|");
            else
                Rprintf("-");
        }
        Rprintf("\n");
    }
}

void
StatusBar::showProgress(int iteration)
{
    if(print)
    {
        if((iteration + 1) % barInterval == 0)
        {
            Rprintf("=");
        }
    }
}

// *************************************************************************************

// ModelPrior

// 29/6/2011: revise model probabilities
// 13/10/2011: include "dep.linear" model prior, which is the 'old' "dependent" model prior.

double
ModelPrior::getLogPrior(const ModelPar& mod) const
{
    if (type == "exponential")
    {
        return - 2.0 * std::accumulate(mod.config.begin(), mod.config.end(), 0);
    }
    else if (type == "independent")
    {
        const double exclLogProb = log(0.5);
        const double inclDegreeLogProb = log(0.5 / nInclDegrees);

        double ret = 0.0;
        LogicalVector::iterator j = continuous.begin();
        for(IntVector::const_iterator i = mod.degIndex.begin();
                i != mod.degIndex.end();
                ++i, ++j)
        {
            double inc;

            if(*j) // this is a continuous covariate
            {
                if(*i == 0)
                    inc = exclLogProb;
                else
                    inc = inclDegreeLogProb;
            }
            else
            {
                inc = exclLogProb; // only because this is log(0.5)!
            }

            ret += inc;
        }
        return ret;
    }
    else if (type == "dependent")
    {
        // determine number of included covariates and which are nonlinear:
        int nIncluded = 0;
        int nSplinePossible = 0;

        LogicalVector::iterator j = continuous.begin();
        for(IntVector::const_iterator i = mod.degIndex.begin();
                i != mod.degIndex.end();
                ++i, ++j)
        {
            if (*i > 0)
            {
                ++nIncluded;

                if(*j)
                {
                    ++nSplinePossible;
                }
            }
        }

        double result = - Rf_lchoose(nCovs, nIncluded) - log1p(nCovs) -
                nSplinePossible * log(static_cast<double>(nInclDegrees));

        return result;
    }
    else if (type == "dep.linear")
    {
        // determine number of included covariates and which are nonlinear:
        int nIncluded = 0;
        int nSplinePossible = 0;
        int nSpline = 0;
        LogicalVector::iterator j = continuous.begin();
        for(IntVector::const_iterator i = mod.degIndex.begin();
                i != mod.degIndex.end();
                ++i, ++j)
        {
            if (*i > 0)
            {
                ++nIncluded;

                if(*j)
                {
                    ++nSplinePossible;

                    if (*i > 1)
                    {
                        ++nSpline;
                    }
                }
            }
        }

        double result = - Rf_lchoose(nCovs, nIncluded) - log1p(nCovs) -
                Rf_lchoose(nSplinePossible, nSpline) - log1p(nSplinePossible) -
                nSpline * log(static_cast<double>(nSplineDegrees));

        return result;
    }
    else // if (type == "flat")
    {
        return 0.0;
    }
}

// *************************************************************************************

// ModelInfo
// convert to an R list
List
ModelInfo::convert2list(long double logNormConst) const
{
    return List::create(_["logMargLik"] = logMargLik,
                        _["logPrior"] = logPrior,
                        _["post"] = exp(logPost - logNormConst),
                        _["hits"] = hits);
}

// *************************************************************************************

// ModelCache

// insert model parameter and belonging model info into the cache.
// returns false if not inserted (e.g. because the par was
// already inside, or the model was not good enough)
bool
ModelCache::insert(const ModelPar& par,
                      const ModelInfo& info)
{
    // first check size of cache
    if(isFull())
    {
        // if we are full, then check if this log posterior is better than
        // the worst cached model, which is pointed to by
        MapType::iterator worstModelIter = *(modelIterSet.begin());

        // the comparison
        if((worstModelIter->second) < info)
        {
            // new model is better than worst model cached.
            // so we delete the worst model from the cache.

            // first from the map
            modelMap.erase(worstModelIter);
            // and then from the set
            modelIterSet.erase(modelIterSet.begin());
        }
        else
        {
            // the new model is not better than the worst model cached,
            // so we do not cache it.
            return false;
        }
    }

    // so now we know that we want to insert the model into the cache,
    // either because the cache was not full or because the new model was better
    // than the worst model cached.

    // -> try inserting into the map:
    std::pair<MapType::iterator, bool> ret = modelMap.insert(MapType::value_type(par, info));

    // if we were successful:
    if(ret.second)
    {
        // then also insert the iterator pointing to the map element into the set.
        modelIterSet.insert(ret.first);

        // return success
        return true;
    }
    else
    {
        return false;
        Rf_error("Should not happen: model already contained in model cache!");
    }
}

// search for the model info of a model config in the map,
// and return false if not found
bool
ModelCache::getModelInfo(const ModelPar& par,
                         ModelInfo& info) const
{
    // search for the config in the map
    MapType::const_iterator ret = modelMap.find(par);

    // if found, return the log marg lik
    if(ret != modelMap.end())
    {
        info = ret->second;
        return true;
    }
    else
        return false;
}

// increment the sampling frequency for a model configuration
// (of course, if this config is not cached nothing is done!)
void
ModelCache::incrementFrequency(const ModelPar& par)
{
    // search for the config in the map
    MapType::iterator ret = modelMap.find(par);

    // if found, increment the hits
    if(ret != modelMap.end())
        ret->second.hits++;
}

// compute the log normalising constant from all cached models
long double
ModelCache::getLogNormConstant() const
{
    // use safe summation
    SafeSum vec;

    // traverse the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // and add all unnormalized log posteriors
        vec.add(m->second.logPost);
    }

    // return the log of the sum of the exp'ed saved elements
    return vec.logSumExp();
}

// compute the inclusion probabilities from all cached models,
// taking the log normalising constant and the total number of FPs / UC groups
List
ModelCache::getInclusionProbs(long double logNormConstant, int nCovs) const
{
    // abbreviation
    typedef std::vector<SafeSum> SafeSumVector;

    // allocate vector of safeSum objects for all covariates
    SafeSumVector included(nCovs);
    SafeSumVector smooth(nCovs);

    // now process each model in the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // abbrevs
        const ModelPar& thisPar = m->first;
        const ModelInfo& thisInfo = m->second;

        SafeSumVector::iterator i = included.begin();
        SafeSumVector::iterator s = smooth.begin();
        for (int c = 0; c < nCovs; ++c, ++i, ++s)
        {
            // is this covariate included in the model m?
            if (thisPar.config[c] > 0.0)
            {
                // then add the normalized model probability onto his inclusion stack
                double normProb = exp(thisInfo.logPost - logNormConstant);
                i->add(normProb);

                // is it included smoothly?
                if(thisPar.config[c] > 1.0)
                {
                    // then add the normalized model probability onto his smooth stack
                    s->add(normProb);
                }
            }
        }

    } // end processing all models in the cache

    // so now we can sum up safesum-wise to the return double vectors:

    NumericVector includedProbs;
    for(SafeSumVector::iterator
            i = included.begin();
            i != included.end();
            ++i)
    {
        includedProbs.push_back(i->sum());
    }

    NumericVector smoothProbs;
    for(SafeSumVector::iterator
            s = smooth.begin();
            s != smooth.end();
            ++s)
    {
        smoothProbs.push_back(s->sum());
    }

    return List::create(_["included"] = includedProbs,
                        _["smooth"] = smoothProbs);
}

// convert the best nModels from the cache into an R list
List
ModelCache::getListOfBestModels(int nModels, long double logNormConst) const
{
    // allocate the return R-list
    List ret(std::min(nModels, static_cast<int>(modelIterSet.size())));

    // process the ordered list of best models from the end (because the set is ordered increasingly)
    int i = 0;
    for(SetType::const_reverse_iterator
            s = modelIterSet.rbegin();
            i < ret.length();
            ++s, ++i)
    {
        // allocate two-element list in the i-th slot of the return list
        ret[i] = List::create(_["configuration"] = (**s).first.config,
                              _["information"] = (**s).second.convert2list(logNormConst));
    }

    return ret;
}

// *************************************************************************************
