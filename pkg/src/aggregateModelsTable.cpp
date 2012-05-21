/*
 * aggregateModelsTable.cpp
 *
 *  Created on: 04.08.2011
 *      Author: daniel
 */

#include <rcppExport.h>
#include <types.h>

#include <map>
#include <string>
#include <sstream>

using namespace Rcpp;
//using namespace std;


// 25/04/2012: erase last "-" from meta model strings
SEXP cpp_aggregateModelsTable(SEXP R_modelsTable,
                              SEXP R_posterior,
                              SEXP R_cut)
{
    // coerce arguments
    List modelsTable = as<List>(R_modelsTable);
    NumericVector posterior = as<NumericVector>(R_posterior);
    int cut = as<int>(R_cut);

    // the map of meta-models with their posterior probabilities
    typedef std::map<std::string, double> MapType;
    MapType metaModels;

    // ...today: iterator action!!

    // construct vector of iterators for separate columns of the modelsTable
    typedef std::vector<IntegerVector::iterator> IterVector;
    IterVector configIters;
    for(List::iterator
            i_column = modelsTable.begin();
            i_column != modelsTable.end();
            ++i_column)
    {
        configIters.push_back(as<IntegerVector>(*i_column).begin());
    }

    // here the meta configurations of the models will be saved
    std::vector<std::string> metaConfigs;

    // now iterate through all model configs
    for(NumericVector::iterator
            i_post = posterior.begin();
            i_post != posterior.end();
            ++i_post)
    {
        // get the meta-configuration vector for this model
        std::ostringstream stream;
        for(IterVector::iterator
                i_iter = configIters.begin();
                i_iter != configIters.end();
                ++i_iter)
        {
            // push back into config vector
            const int metaDf = **i_iter < cut ? **i_iter : cut;
            stream << metaDf << "-";

            // and go on to next element for next time
            ++(*i_iter);
        }

        // convert to string
        std::string thisMetaConfig = stream.str();

        // erase the last character, which is "-"
        thisMetaConfig.erase(thisMetaConfig.size() - 1);

        // insert into vector
        metaConfigs.push_back(thisMetaConfig);

        // check if it is already contained in the map
        MapType::iterator found = metaModels.find(thisMetaConfig);
        if(found != metaModels.end())
        {
            // if found, add the posterior probability of this sub-model
            found->second += *i_post;
        }
        else
        {
            // insert into map
            metaModels.insert(MapType::value_type(thisMetaConfig, *i_post));
        }
    }

    // return the result
    return List::create(_["metaConfig"] = metaConfigs,
                        _["metaProb"] = metaModels);
}
