/*
 * hyp2f1.h
 *
 *  Created on: 14.10.2010
 *      Author: daniel
 */

#ifndef HYP2F1_H_
#define HYP2F1_H_

//double
//hyp2f1(double a,
//       double b,
//       double c,
//       double x);

// helper function which tries to compute the exact value first,
// and if that is not finite uses the Laplace approximation.
double
log_hyp2f1(double a,
           double b,
           double c,
           double x);

#endif /* HYP2F1_H_ */
