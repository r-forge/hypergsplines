/*
 * getRho.h
 *
 *  Created on: 28.07.2011
 *      Author: daniel
 */

#ifndef GETRHO_H_
#define GETRHO_H_

#include <rcppExport.h>

double
getRho(const Rcpp::NumericVector& lambdas,
       double degree);

#endif /* GETRHO_H_ */
