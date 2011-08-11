/*
 * getLogMargLikEst.h
 *
 *  Created on: 07.04.2011
 *      Author: daniel
 */

#ifndef GETLOGMARGLIKEST_H_
#define GETLOGMARGLIKEST_H_

#include <types.h>

double
getLogMargLikEst(const DoubleVector& numeratorTerms,
                 const DoubleVector& denominatorTerms,
                 double highDensityPointLogUnPosterior,
                 double& se,
                 PosInt endLag = 40);

#endif /* GETLOGMARGLIKEST_H_ */
