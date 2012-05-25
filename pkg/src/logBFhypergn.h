/*
 * logBFhypergn.h
 *
 *  Created on: 30.03.2012
 *      Author: daniel
 */

#ifndef LOGBFHYPERGN_H_
#define LOGBFHYPERGN_H_

// compute the log BF against the null model under the hyper-g/n prior,
// either with the Appell function or with the Laplace approximation
double
logBFhypergn(int n,
             int p,
             double R2);

#endif /* LOGBFHYPERGN_H_ */
