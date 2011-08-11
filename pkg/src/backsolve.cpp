/*
 * backsolve.cpp
 *
 *  Created on: 14.10.2010
 *      Author: daniel
 */

#include <backsolve.h>
#include <R_ext/BLAS.h>

// modified from R/src/appl/bakslv.c:
int
bakslv(const double * t,
       const int ldt,
       const int n,
       const double * b,
       const int ldb,
       const int nb,
       double * x,
       const bool upperTri,
       const bool transpose)
{

/* bakslv is a subroutine to solve triangular systems of
 * the form
 *                   t * x = b
 * or
 *                   t' * x = b         [ t' := transpose(t) ]

 * where t is a triangular matrix of order n.
 * The subroutine handles the multiple right-hand side case.


 * on entry

 *      t      double (ldt,n'). n' >= n (below)
 *             t[] contains the coefficient matrix of the system
 *             to be solved.  only the elements above or below
 *             the diagonal are referenced.

 *      ldt    int;     ldt is the leading dimension of the array t.
 *      n      int;     n is the order of the system.  n <= min(ldt,ldb)


 *      b      double (ldb,nb').  nb' >= nb (below)
 *             b[] contains the right hand side(s) of the system.

 *      ldb    int;     ldb is the leading dimension of the array b.
 *      nb     int;     the number of right hand sides of the system.

 *      upperTri, transpose (bool);     specifies what kind of system is to be solved.
 *

 * on return

 *      x      double precision(ldb, nb)
 *             contains the solution(s) if info == 0.

 *      info   int
 *             info contains zero if the system is nonsingular.
 *             otherwise info contains the index of
 *             the first zero diagonal element of t.

 */

    const char * transa = transpose ? "T" : "N";
    const char * uplo = upperTri ? "U" : "L";
    const char * side = "L";
    const char * diag = "N";

    const int ione = 1;
    const double one = 1.0;

    // check for zeros on diagonal:
    for(int i = 0; i < n; i++)
    {
        if (t[i * (ldt + 1)] == 0.0)
        {
            return i + 1;
        }
    }

    // copy b to x:
    for(int j = 0; j < nb; j++)
    {
        F77_CALL(dcopy)(&n, &b[j * ldb], &ione, &x[j * ldb], &ione);
    }

    // then call the level-3 BLAS routine:
    if (n > 0 && nb > 0 && ldt > 0 && ldb > 0)
    {
        F77_CALL(dtrsm)(side, uplo, transa, diag, &n, &nb, &one,
                        t, &ldt, x, &ldb);
    }

    // everything went OK
    return 0;
}


// the frontend function:
arma::mat
backsolve(const arma::mat& r,
          const arma::mat& x,
          const bool upperTri,
          const bool transpose)
{
    // check and extract dimensions:

    const int k = r.n_cols;

    if (k <= 0 || x.n_rows < k)
    {
        Rf_error("invalid argument values in 'backsolve'");
    }

    const int nb = x.n_cols;

    // allocate space for the result matrix
    arma::mat ret(k, nb);

    // run the workhorse function
    int info = bakslv(r.memptr(),
                      r.n_rows,
                      k,
                      x.memptr(),
                      k,
                      nb,
                      ret.memptr(),
                      upperTri,
                      transpose);

    // check the error code
    if(info != 0)
    {
        Rf_error("singular matrix in 'backsolve'. First zero in diagonal [%d]",
                 info);
    }

    // return the result matrix
    return ret;
}
