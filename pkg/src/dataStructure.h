/*
 * dataStructure.h
 *
 *  Created on: 14.03.2011
 *      Author: daniel
 *
 *  30/03/2012: remove hyperparameter a, instead use string for g-prior
 */

#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <rcppExport.h>

#include <types.h>
#include <distributions.h>
#include <gpriors.h>
#include <links.h>
#include <string>




// GaussHermite: encapsulate the Gauss-Hermite integration stuff
class GaussHermite {
public:

    // compute nodes and log weights for given mode and sd of target unnormalized density
    void
    getNodesAndLogWeights(double mode, double var,
                          DoubleVector& nodes, DoubleVector& logWeights) const; // output

    // constructor
    GaussHermite(Rcpp::List rcpp_gaussHermite) :
        tVec(Rcpp::as<DoubleVector>(rcpp_gaussHermite["nodes"])),
        wVec(Rcpp::as<DoubleVector>(rcpp_gaussHermite["weights"]))
        {
        }

private:
    const DoubleVector tVec; // nodes
    const DoubleVector wVec; // weights (not log!)
};

// Computation: capture all computation options

// 19/5/2011: generalise "binaryLogisticCorrection" to "higherOrderCorrection"
class Computation {
public:
    const bool verbose;
    const bool debug;
    const PosInt nGaussHermite;
    const GaussHermite gaussHermite;
    const bool useOpenMP;
    const bool higherOrderCorrection;

    Computation(SEXP R_computation) :
        verbose(Rcpp::as<Rcpp::List>(R_computation)["verbose"]),
        debug(Rcpp::as<Rcpp::List>(R_computation)["debug"]),
        nGaussHermite(Rcpp::as<Rcpp::List>(R_computation)["nGaussHermite"]),
        gaussHermite(Rcpp::as<SEXP>(Rcpp::as<Rcpp::List>(R_computation)["gaussHermite"])),
        useOpenMP(Rcpp::as<Rcpp::List>(R_computation)["useOpenMP"]),
        higherOrderCorrection(Rcpp::as<Rcpp::List>(R_computation)["higherOrderCorrection"])
    {}

private:
};

// McmcOptions: capture all MCMC options
class McmcOptions {
public:
    const PosInt samples;
    const PosInt burnin;
    const PosInt step;
    const PosInt iterations;
    const PosInt nIwlsIterations;

    McmcOptions(SEXP R_mcmc) :
        samples(Rcpp::as<PosInt>(Rcpp::as<Rcpp::List>(R_mcmc)["samples"])),
        burnin(Rcpp::as<PosInt>(Rcpp::as<Rcpp::List>(R_mcmc)["burnin"])),
        step(Rcpp::as<PosInt>(Rcpp::as<Rcpp::List>(R_mcmc)["step"])),
        iterations(Rcpp::as<PosInt>(Rcpp::as<Rcpp::List>(R_mcmc)["iterations"])),
        nIwlsIterations(Rcpp::as<PosInt>(Rcpp::as<Rcpp::List>(R_mcmc)["nIwlsIterations"]))
    {}

private:
};

// only one parameter set for the GLM
struct Parameter
{
    // ctr
    Parameter(const AVector& coefs,
              double z) :
                  coefs(coefs),
                  z(z)
                  {
                  }

    // default ctr
    Parameter(PosInt nCoefs) :
        coefs(nCoefs),
        z(R_NaReal)
        {
        }

    // the sampled regression coefficients
    // (comprising the intercept, the fixed and the spline coefficients)
    AVector coefs;

    // the sampled log covariance factor
    double z;
};

// all info from the iwls results

// 19/5/2011 modified to accomodate offsets
struct IwlsResults
{
    // ctr
    IwlsResults(const AVector& linPred,
                PosInt nCoefs) :
        linPred(linPred),
        coefs(nCoefs),
        qFactor(nCoefs, nCoefs)
        {
        }

    // second ctr
    IwlsResults(PosInt nObs,
                PosInt nCoefs) :
        linPred(nObs),
        coefs(nCoefs),
        qFactor(nCoefs, nCoefs)
        {
        }

    // default ctr
    IwlsResults(){}

    // the linear predictor
    // (note that this includes offsets!)
    AVector linPred;

    // the corresponding coefficient vector,
    AVector coefs;

    // lower-triangular Cholesky factor of the precision matrix
    AMatrix qFactor;

    // the log determinant of the precision matrix
    double logPrecisionDeterminant;
};

// model parameter object: must have a strict weak ordering
struct ModelPar
{
    // for each covariate we have the degree of freedom in here:
    DoubleVector config;

    // and we also have the degree indices in here
    IntVector degIndex;

    // compute degree index vector, given a degree vector
    void
    compDegIndex(const IntVector& degrees);

    // start with an empty (null) model:
    ModelPar(int nCovs) :
        config(nCovs, 0.0),
        degIndex(nCovs, 0)
    {
    }

    // copy integer degrees of freedom configuration
    ModelPar(const Rcpp::IntegerVector& config) :
        config(Rcpp::as<DoubleVector>(config))
    {
    }

    // copy double degrees of freedom configuration
    ModelPar(const Rcpp::NumericVector& config) :
        config(Rcpp::as<DoubleVector>(config))
    {
    }

    // is this the null model?
    bool
    isNullModel() const
    {
        return std::accumulate(config.begin(), config.end(), 0.0) == 0.0;
    }

    // comparison of configurations
    bool
    operator<(const ModelPar& m) const
    {
        return config < m.config;
    }

    // return a textual description of this model configuration
    std::string
    print() const;

};

// ModelData: this wraps the general input like X, y etc.
class ModelData
{
public:

    Rcpp::List modelData;
    Rcpp::List rhoList;
    Rcpp::List zList;
    Rcpp::List lambdasList;

    const Rcpp::NumericMatrix Xfull;
    const double * Xfull_ptr;

    typedef std::vector<arma::mat> MatVector;
    typedef std::vector<MatVector> MatArray;
    MatArray ZtZarray;

    const int nObs;
    const int nCovs;
    const int dimSplineBasis;

    const Rcpp::LogicalVector continuous;

    const std::string gPriorString;

    // the vector of degrees. assumed to be sorted increasingly,
    // starting with 0, 1 and continuing with the "spline" degrees.
    const IntVector degrees;
    // number of different degrees
    const int nDegrees;

    const arma::colvec y;
    const arma::colvec y_centered;
    const double yCenterNormSq;

    // indices of continuous covariates
    IntVector contCovs;
    int nContCovs;

    // ctr: construct from R list
    ModelData(SEXP R_modelData);

    // compute the log marginal likelihood (and as a byproduct the R2)
    // of a specific model
    double
    getLogMargLik(const ModelPar& modPar,
                  double& R2);

private:

};



// GlmModelData: this wraps the general input like X, y etc.
// todo: later on, this could be simplified by making this a child of a general
// ModelData class!

// 19/5/2011:  save also the names of link and distribution
class GlmModelData
{
public:

    Rcpp::List rcpp_modelData;
    Rcpp::List rhoList;
    Rcpp::List zList;
    Rcpp::List lambdasList;
    Rcpp::List familyList;

    const Rcpp::NumericMatrix Xfull;
    const double * Xfull_ptr;

    typedef std::vector<arma::mat> MatVector;
    typedef std::vector<MatVector> MatArray;
    MatArray ZtZarray;

    const int nObs;
    const int nCovs;
    const int dimSplineBasis;

    const Rcpp::LogicalVector continuous;

    const std::string gPriorString;

    // the vector of degrees. assumed to be sorted increasingly,
    // starting with 0, 1 and continuing with the "spline" degrees.
    const IntVector degrees;
    // number of different degrees
    const int nDegrees;

    const arma::colvec y;
    const arma::colvec y_centered;
    const double yCenterNormSq;

    // indices of continuous covariates
    IntVector contCovs;
    int nContCovs;

    // the weight matrix entries
    const AVector weightMatrixEntries;
    const AVector sqrtWeightMatrixEntries;

    // the start linear predictor vector
    const AVector linPredStart;

    // the vector of dispersions (phi / weights)
    const AVector dispersions;

    // the vector of offsets
    const AVector offsets;

    // the log marginal likelihood in the null model
    const double nullModelLogMargLik;

    // the g-prior information
    const GPrior* gPrior;

    // the link information
    const Link* link;

    // the distribution information
    const Distribution* distribution;

    // save also the names of link and distribution
    const std::string familyString;
    const std::string linkString;

    // does this model use the canonical link?
    const bool canonicalLink;

    // constructor
    GlmModelData(SEXP R_modelData);

    // destructor
    ~GlmModelData()
    {
        delete gPrior;
        delete link;
        delete distribution;
    }


private:

};

// status bar
class StatusBar
{
public:
    // ctr
    StatusBar(bool print,
              int nBars,
              int numberGap,
              int totaliters);

    // progress function
    void
    showProgress(int iteration);

    // dtr
    ~StatusBar()
    {
        if(print)
            Rprintf("\n");
    }

private:
    const bool print;
    const int nBars;
    const int numberGap;
    const int totaliters;
    const int barInterval;
};

// implements the model prior
class ModelPrior
{
public:
    ModelPrior(SEXP R_modelPrior,
               int nCovs,
               int nDegrees,
               const Rcpp::LogicalVector& continuous) :
                   nCovs(nCovs),
                   nInclDegrees(nDegrees - 1),
                   nSplineDegrees(nDegrees - 2),
                   continuous(continuous),
                   type(Rcpp::as<Rcpp::StringVector>(R_modelPrior)[0])
    {
    }

    double
    getLogPrior(const ModelPar& mod) const;

private:
    const int nCovs;
    const int nInclDegrees;
    const int nSplineDegrees;
    const Rcpp::LogicalVector continuous;
    std::string type;
};

// data structure: model information for general models.
// keep in mind that contents must be assignable.
struct ModelInfo
{
    // most important:
    double logMargLik;
    double logPrior;

    // just the sum of logMargLik and logPrior.
    // But we need it often, so save it here (the ctr computes it).
    double logPost;

    // only needed for MCMC:
    int hits;

    // simple constructor
    ModelInfo(double logMargLik,
              double logPrior) :
        logMargLik(logMargLik), logPrior(logPrior), logPost(logMargLik
                + logPrior), hits(0)
    {
    }

    // compare two model infos bei their posterior probability:
    // the model info with the lower one is lower.
    bool
    operator<(const ModelInfo& m) const
    {
        return logPost < m.logPost;
    }

    // convert to an R list
    Rcpp::List
    convert2list(long double logNormConst) const;
};


// Caches the best models in a map of a given maximum size, and also stores the
// (unnormalized) log posterior probabilities in an ordered set, pointing to the models in the map.
class ModelCache
{
public:

    // create a new ModelCache with given maximum size.
    ModelCache(int maxSize) :
        nanCounter(0),
        maxSize(maxSize),
        modelMap(),
        modelIterSet()
        {
        }

    // check if max size was reached
    bool
    isFull() const
    {
        return modelMap.size() == maxSize;
    }

    // return size of cache
    int
    size() const
    {
        return modelMap.size();
    }

    // insert model parameter and belonging model info into the cache.
    // returns false if not inserted (e.g. because the par was
    // already inside, or the model was not good enough)
    bool
    insert(const ModelPar& par, const ModelInfo& info);


    // search for the model info of a model config in the map,
    // and return false if not found
    bool
    getModelInfo(const ModelPar& par,
                 ModelInfo& info) const;

    // increment the sampling frequency for a model configuration
    // (of course, if this config is not cached nothing is done!)
    void
    incrementFrequency(const ModelPar& par);

    // compute the log normalising constant from all cached models
    long double
    getLogNormConstant() const;

    // compute the inclusion probabilities from all cached models,
    // taking the log normalising constant and the total number of FPs / UC groups
    Rcpp::List
    getInclusionProbs(long double logNormConstant, int nCovs) const;

    // convert the best nModels from the cache into an R list
    Rcpp::List
    getListOfBestModels(int nModels, long double logNormConst) const;

    // for counting the number of times a NaN model was encountered
    int nanCounter;

private:
    // the map type
    typedef std::map<ModelPar, ModelInfo> MapType;

    // define comparison function for iterators
    struct Compare_map_iterators
    {
        bool
        operator()(const MapType::iterator& first, const MapType::iterator& second) const
        {
            return (first->second.logPost) < (second->second.logPost);
        }
    };

    // the set type of ordered map iterators
    typedef std::set<MapType::iterator, Compare_map_iterators> SetType;

    // and finally the data members
    const MapType::size_type maxSize;
    MapType modelMap;
    SetType modelIterSet;
};


// helper class for stochastic search
struct ModelMcmc
{
    // initialize with a given model config and log marg lik
    ModelMcmc(const ModelPar& modPar, double logMargLik, double logPrior) :
        modPar(modPar), logMargLik(logMargLik), logPrior(logPrior)
    {
    }

    ModelPar modPar;

    double logMargLik;
    double logPrior;
};


#endif /* DATASTRUCTURE_H_ */
