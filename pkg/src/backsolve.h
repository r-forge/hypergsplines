/*
 * backsolve.h
 *
 *  Created on: 14.10.2010
 *      Author: daniel
 */

#ifndef BACKSOLVE_H_
#define BACKSOLVE_H_

#include <RcppArmadillo.h>

arma::mat
backsolve(const arma::mat& r,
          const arma::mat& x,
          const bool upperTri,
          const bool transpose);


#endif /* BACKSOLVE_H_ */
