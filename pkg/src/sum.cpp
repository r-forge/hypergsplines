/*
 * sum.cpp
 *
 *  Created on: 21.10.2010
 *      Author: daniel
 */

#include <sum.h>
#include <rcppExport.h>

// ***************************************************************************************************//

// safeSum //

void SafeSum::add(const long double &val)
{
        vals.push_back(val);
}

long double SafeSum::sum()
{
        long double ret = modified_deflation(vals);
        return ret;
}

// compute the log of the sum of the exp of the elements using accurate algorithm,
// and avoiding infinite contributions.
long double SafeSum::logSumExp()
{
    // be sure that there is at least 1 value in "vals"
    if(vals.empty())
    {
        return R_NaN;
    }

    // the maximum of the log contributions is:
    long double maxLogContrib = *std::max_element(vals.begin(), vals.end());

    // now compute the constant which is added to all log contributions,
    // in order to avoid infinite contributions and at the same time use
    // the whole number space (i.e. possibly avoid zero contributions)
    long double constant = logl(LDBL_MAX) - 100.0L - maxLogContrib;
    // 100 is for safety.

    // so now the contributions, offset by the constant
    LongDoubleVector expVals;
    for(LongDoubleVector::const_iterator
            l = vals.begin();
            l != vals.end();
            ++l)
    {
        expVals.push_back(expl(*l + constant));
    }

    // the result is the log of the sum, corrected with the constant:
    long double ret = logl(modified_deflation(expVals)) - constant;
    return ret;
}

long double SafeSum::simpleSum()
{
        long double ret = 0.0;
        for(LongDoubleVector::const_iterator
                v = vals.begin();
                v != vals.end();
                ++v)
        {
            ret += *v;
        }
        return ret;
}

// ***************************************************************************************************//

// indexSafeSum //

void IndexSafeSum::add(const Ind& ind)
{
        indices.insert(ind);
}

long double
IndexSafeSum::sumNormalizedExp(const SafeSum& s, long double logNormConst) const
{
    LongDoubleVector tempVec;
    for (IndSet::const_iterator i = indices.begin(); i != indices.end(); i++)
    {
        tempVec.push_back(expl(s.vals.at(* i) - logNormConst));
    }
    return modified_deflation(tempVec);
}

// ***************************************************************************************************//
