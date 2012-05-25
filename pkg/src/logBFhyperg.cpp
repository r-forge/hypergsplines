/*
 * logBFhyperg.cpp
 *
 *  Created on: 30.03.2012
 *      Author: daniel
 */

#include <hyp2f1.h>
#include <logBFhyperg.h>

#include <Rmath.h>

// compute the log BF against the null model under the hyper-g prior
double
logBFhyperg(int n,
            int p,
            double R2)
{
    double ret = - log(p + 2.0)
            + log_hyp2f1((n - 1.0) / 2.0,
                         1.0,
                         (p + 4.0) / 2.0,
                         R2);

    return ret;
}
