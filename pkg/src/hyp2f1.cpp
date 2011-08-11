/*
 * hyp2f1.cpp
 *
 *  Created on: 13.09.2010
 *      Author: daniel
 */

/*
 *
 *      Gauss hypergeometric function   F
 *                                     2 1
 *
 *
 * SYNOPSIS:
 *
 * double a, b, c, x, y, hyp2f1();
 *
 * y = hyp2f1( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 *  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
 *                           2 1
 *
 *           inf.
 *            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
 *   =  1 +   >   -----------------------------  x   .
 *            -         c(c+1)...(c+k) (k+1)!
 *          k = 0
 *
 *  Cases addressed are
 *      Tests and escapes for negative integer a, b, or c
 *      Linear transformation if c - a or c - b negative integer
 *      Special case c = a or c = b
 *      Linear transformation for  x near +1
 *      Transformation for x < -0.5
 *      Psi function expansion if x > 0.5 and c - a - b integer
 *      Conditionally, a recurrence on c to make c-a-b > 0
 *
 * |x| > 1 is rejected.
 *
 * The parameters a, b, c are considered to be integer
 * valued if they are within 1.0e-14 of the nearest integer
 * (1.0e-13 for IEEE arithmetic).
 *
 * ACCURACY:
 *
 *
 *               Relative error (-1 < x < 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1,7        230000      1.2e-11     5.2e-14
 *
 * Several special cases also tested with a, b, c in
 * the range -7 to 7.
 *
 * ERROR MESSAGES:
 *
 * A "partial loss of precision" message is printed if
 * the internally estimated relative error exceeds 1^-12.
 * A "singularity" message is printed on overflow or
 * in cases not addressed (such as x < -1).
 */



/*
 *
 * modified from:
 *
 Cephes Math Library Release 2.8:  June, 2000
 Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
 *
 * by Daniel Sabanés Bové 2010.
 *
 */

#include <rcppExport.h>
#include <string>

#include <hyp2f1.h>

using std::string;

#define gamma Rf_gammafn
#define lgamma Rf_lgammafn
#define psi Rf_digamma
#define pow R_pow
#define powi R_pow_di

const double EPS = 1.0e-13;
const double EPS2 = 1.0e-10;
const double ETHRESH = 1.0e-12;
const double MACHEP = DBL_EPSILON;
const double MAXNUM = DBL_MAX;

static const string ermsg[7] =
{ "unknown", "domain", "singularity", "overflow", "underflow",
        "total loss of precision", "partial loss of precision" };

static int
mtherr(string name,
       int code)
{

    /* Display string passed by calling program,
     * which is supposed to be the name of the
     * function in which the error occurred:
     */
    Rprintf("\n%s ",
            name.c_str());

    /* Display error message defined
     * by the code argument.
     */
    if ((code <= 0) || (code >= 7))
        code = 0;
    Rprintf("%s error\n",
            ermsg[code].c_str());

    /* Return to calling
     * program
     */
    return (0);
}

static double
hyt2f1(double,
       double,
       double,
       double,
       double *);
static double
hys2f1(double,
       double,
       double,
       double,
       double *);

double
hyp2f1(double a,
       double b,
       double c,
       double x)
{
    double d, d1, d2, e;
    double p, q, r, s, y, ax;
    double ia, ib, ic, id, err;
    int flag, i, aid;

    err = 0.0;
    ax = fabs(x);
    s = 1.0 - x;
    flag = 0;
    ia = round(a); /* nearest integer to a */
    ib = round(b);

    if (a <= 0)
    {
        if (fabs(a - ia) < EPS) /* a is a negative integer */
            flag |= 1;
    }

    if (b <= 0)
    {
        if (fabs(b - ib) < EPS) /* b is a negative integer */
            flag |= 2;
    }

    if (ax < 1.0)
    {
        if (fabs(b - c) < EPS) /* b = c */
        {
            y = pow(s,
                    - a); /* s to the -a power */
            goto hypdon;
        }
        if (fabs(a - c) < EPS) /* a = c */
        {
            y = pow(s,
                    - b); /* s to the -b power */
            goto hypdon;
        }
    }

    if (c <= 0.0)
    {
        ic = round(c); /* nearest integer to c */
        if (fabs(c - ic) < EPS) /* c is a negative integer */
        {
            /* check if termination before explosion */
            if ((flag & 1) && (ia > ic))
                goto hypok;
            if ((flag & 2) && (ib > ic))
                goto hypok;
            goto hypdiv;
        }
    }

    if (flag) /* function is a polynomial */
        goto hypok;

    if (ax > 1.0) /* series diverges    */
        goto hypdiv;

    p = c - a;
    ia = round(p); /* nearest integer to c-a */
    if ((ia <= 0.0) && (fabs(p - ia) < EPS)) /* negative int c - a */
        flag |= 4;

    r = c - b;
    ib = round(r); /* nearest integer to c-b */
    if ((ib <= 0.0) && (fabs(r - ib) < EPS)) /* negative int c - b */
        flag |= 8;

    d = c - a - b;
    id = round(d); /* nearest integer to d */
    q = fabs(d - id);

    /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
     * for reporting a bug here.  */
    if (fabs(ax - 1.0) < EPS) /* |x| == 1.0     */
    {
        if (x > 0.0)
        {
            if (flag & 12) /* negative int c-a or c-b */
            {
                if (d >= 0.0)
                    goto hypf;
                else
                    goto hypdiv;
            }
            if (d <= 0.0)
                goto hypdiv;
            /*          y = gamma(c)*gamma(d)/(gamma(p)*gamma(r)); */
            y = exp(lgamma(c) + lgamma(d) - (lgamma(p) + lgamma(r)));
            goto hypdon;
        }

        if (d <= - 1.0)
            goto hypdiv;

    }

    /* Conditionally make d > 0 by recurrence on c
     * AMS55 #15.2.27
     */
    if (d < 0.0)
    {
        /* Try the power series first */
        y = hyt2f1(a,
                   b,
                   c,
                   x,
                   & err);
        if (err < ETHRESH)
            goto hypdon;
        /* Apply the recurrence if power series fails */
        err = 0.0;
        aid = 2 - id;
        e = c + aid;
        d2 = hyp2f1(a,
                    b,
                    e,
                    x);
        d1 = hyp2f1(a,
                    b,
                    e + 1.0,
                    x);
        q = a + b + 1.0;
        for (i = 0; i < aid; i++)
        {
            r = e - 1.0;
            y = (e * (r - (2.0 * e - q) * x) * d2 + (e - a) * (e - b) * x * d1)
                    / (e * r * s);
            e = r;
            d1 = d2;
            d2 = y;
        }
        goto hypdon;
    }

    if (flag & 12)
        goto hypf; /* negative integer c-a or c-b */

    hypok: y = hyt2f1(a,
                      b,
                      c,
                      x,
                      & err);

    hypdon: if (err > ETHRESH)
    {
        mtherr("hyp2f1",
               6);
        /*      printf( "Estimated err = %.2e\n", err ); */
    }
    return (y);

    /* The transformation for c-a or c-b negative integer
     * AMS55 #15.3.3
     */
    hypf: y = pow(s,
                  d) * hys2f1(c - a,
                              c - b,
                              c,
                              x,
                              & err);
    goto hypdon;

    /* The alarm exit */
    hypdiv: mtherr("hyp2f1",
                   3);
    // Overflow -> Return infinity.
    return (R_PosInf);
}



// R interface function
SEXP
cpp_hyp2f1(SEXP R_interface)
{
    // unpack the arguments
    R_interface = CDR(R_interface);
    SEXP R_a = CAR(R_interface);
    double a = Rf_asReal(R_a);

    R_interface = CDR(R_interface);
    SEXP R_b = CAR(R_interface);
    double b = Rf_asReal(R_b);

    R_interface = CDR(R_interface);
    SEXP R_c = CAR(R_interface);
    double c = Rf_asReal(R_c);

    R_interface = CDR(R_interface);
    SEXP R_x = CAR(R_interface);
    double * x = REAL(R_x);

    // call the C function for all values of R_x:
    SEXP ret;
    Rf_protect(ret = Rf_allocVector(REALSXP, Rf_length(R_x)));

    for(int i = 0; i != Rf_length(R_x); ++i)
    {
        REAL(ret)[i] = hyp2f1(a,
                              b,
                              c,
                              x[i]);
    }

    // unprotect and return the result
    Rf_unprotect(1);
    return ret;
}


double
log_psi_laplace(double a,
                double b,
                double c,
                double x)
{
    // a: the first parameter, will be (n-1)/2
    // b: the second parameter, will be 1
    // c: will be (p + a) / 2
    // x: will be R^2

    // components of quadratic equation in exp(tau) are
    double alpha = 2.0 * (b - c) * (1.0 - x);
    double beta = 2.0 * (2.0 * b - c + x * (a - b));
    double gama = 2.0 * b;

    // exp(hat(tau)) is the mode of g:
    double gMode = (- beta - sqrt(beta * beta - 4.0 * alpha * gama)) / (2.0 * alpha);

    // the log integrand evaluated at the mode is
    double logIntegrandAtTauMode = b * log(gMode) + (a - c) * log1p(gMode) - a * log1p((1 - x) * gMode);

    // the variance is
    double sigma2 = 1.0 /
            (- (a + c) * gMode / pow(1.0 + gMode, 2.0) +
                    a / 2.0 * ((1.0 - x) * gMode) / pow(1.0 + (1.0 - x) * gMode, 2.0));


    // so the Laplace approximation is
    double ret = M_LN_SQRT_2PI + 0.5 * log(sigma2) + logIntegrandAtTauMode;
    return(ret);
}


// log of the Laplace approximation
double
log_hyp2f1_laplace(double a,
                   double b,
                   double c,
                   double x)
{
    return log_psi_laplace(a, b, c, x) - Rf_lbeta(b, c - b);
}

// R interface function to the log of the Laplace approximation
SEXP
cpp_log_hyp2f1_laplace(SEXP R_interface)
{
    // unpack the arguments
    R_interface = CDR(R_interface);
    SEXP R_a = CAR(R_interface);
    double a = Rf_asReal(R_a);

    R_interface = CDR(R_interface);
    SEXP R_b = CAR(R_interface);
    double b = Rf_asReal(R_b);

    R_interface = CDR(R_interface);
    SEXP R_c = CAR(R_interface);
    double c = Rf_asReal(R_c);

    R_interface = CDR(R_interface);
    SEXP R_x = CAR(R_interface);
    double * x = REAL(R_x);

    // call the C function for all values of R_x:
    SEXP ret;
    Rf_protect(ret = Rf_allocVector(REALSXP, Rf_length(R_x)));

    for(int i = 0; i != Rf_length(R_x); ++i)
    {
        REAL(ret)[i] = log_hyp2f1_laplace(a,
                                          b,
                                          c,
                                          x[i]);
    }

    // unprotect and return the result
    Rf_unprotect(1);
    return ret;
}


// helper function which tries to compute the exact value first,
// and if that is not finite uses the Laplace approximation.
double
log_hyp2f1(double a,
           double b,
           double c,
           double x)
{
    double v = hyp2f1(a, b, c, x);

    if(! R_finite(v))
    {
        Rf_warning("Using the Laplace approximation for call hyp2f1(%f, %f, %f, %f) because result was infinite", a, b, c, x);
        return log_hyp2f1_laplace(a, b, c, x);
    }
    else
    {
        return log(v);
    }

}


/* Apply transformations for |x| near 1
 * then call the power series
 */
static double
hyt2f1(double a,
       double b,
       double c,
       double x,
       double * loss)
{
    double p, q, r, s, t, y, d, err, err1;
    double ax, id, d1, d2, e, y1;
    int i, aid;

    err = 0.0;
    s = 1.0 - x;
    if (x < - 0.5)
    {
        if (b > a)
            y = pow(s,
                    - a) * hys2f1(a,
                                  c - b,
                                  c,
                                  - x / s,
                                  & err);

        else
            y = pow(s,
                    - b) * hys2f1(c - a,
                                  b,
                                  c,
                                  - x / s,
                                  & err);

        goto done;
    }

    d = c - a - b;
    id = round(d); /* nearest integer to d */

    if (x > 0.9)
    {
        if (fabs(d - id) > EPS) /* test for integer c-a-b */
        {
            /* Try the power series first */
            y = hys2f1(a,
                       b,
                       c,
                       x,
                       & err);
            if (err < ETHRESH)
                goto done;
            /* If power series fails, then apply AMS55 #15.3.6 */
            q = hys2f1(a,
                       b,
                       1.0 - d,
                       s,
                       & err);
            /*  q *= gamma(d) /(gamma(c-a) * gamma(c-b)); */
            q *= exp(lgamma(d) - (lgamma(c - a) + lgamma(c - b)));
            r = pow(s,
                    d) * hys2f1(c - a,
                                c - b,
                                d + 1.0,
                                s,
                                & err1);
            r *= exp(lgamma(- d) - (lgamma(a) + lgamma(b)));
            y = q + r;

            q = fabs(q); /* estimate cancellation error */
            r = fabs(r);
            if (q > r)
                r = q;
            err += err1 + (MACHEP * r) / y;

            y *= gamma(c);
            goto done;
        }
        else
        {
            /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12 */
            if (id >= 0.0)
            {
                e = d;
                d1 = d;
                d2 = 0.0;
                aid = id;
            }
            else
            {
                e = - d;
                d1 = 0.0;
                d2 = d;
                aid = - id;
            }

            ax = log(s);

            /* sum for t = 0 */
            y = psi(1.0) + psi(1.0 + e) - psi(a + d1) - psi(b + d1) - ax;
            y /= gamma(e + 1.0);

            p = (a + d1) * (b + d1) * s / gamma(e + 2.0); /* Poch for t=1 */
            t = 1.0;
            do
            {
                r = psi(1.0 + t) + psi(1.0 + t + e) - psi(a + t + d1) - psi(b
                        + t + d1) - ax;
                q = p * r;
                y += q;
                p *= s * (a + t + d1) / (t + 1.0);
                p *= (b + t + d1) / (t + 1.0 + e);
                t += 1.0;
            }
            while (fabs(q / y) > EPS);

            if (id == 0.0)
            {
                y *= exp(lgamma(c) - (lgamma(a) + lgamma(b)));
                goto psidon;
            }

            y1 = 1.0;

            if (aid == 1)
                goto nosum;

            t = 0.0;
            p = 1.0;
            for (i = 1; i < aid; i++)
            {
                r = 1.0 - e + t;
                p *= s * (a + t + d2) * (b + t + d2) / r;
                t += 1.0;
                p /= t;
                y1 += p;
            }
            nosum: p = gamma(c);
            y1 *= exp(lgamma(e) + log(p) - (lgamma(a + d1) + lgamma(b + d1)));

            y *= exp(log(p) - (lgamma(a + d2) + lgamma(b + d2)));
            if ((aid & 1) != 0)
                y = - y;

            q = pow(s,
                    id); /* s to the id power */
            if (id > 0.0)
                y *= q;
            else
                y1 *= q;

            y += y1;
            psidon: goto done;
        }

    }

    /* Use defining power series if no special cases */
    y = hys2f1(a,
               b,
               c,
               x,
               & err);

    done: * loss = err;
    return (y);
}

/* Defining power series expansion of Gauss hypergeometric function */
static double
hys2f1(double a,
       double b,
       double c,
       double x,
       double * loss)
{
    double f, g, h, k, m, s, u, umax;
    int i;

    i = 0;
    umax = 0.0;
    f = a;
    g = b;
    h = c;
    s = 1.0;
    u = 1.0;
    k = 0.0;
    do
    {
        if (fabs(h) < EPS)
        {
            * loss = 1.0;
            return (MAXNUM);
        }
        m = k + 1.0;
        u = u * ((f + k) * (g + k) * x / ((h + k) * m));
        s += u;
        k = fabs(u); /* remember largest term summed */
        if (k > umax)
            umax = k;
        k = m;
        if (++i > 10000) /* should never happen */
        {
            * loss = 1.0;
            return (s);
        }
    }
    while (fabs(u / s) > MACHEP);

    /* return estimated relative error */
    * loss = (MACHEP * umax) / fabs(s) + (MACHEP * i);

    return (s);
}
