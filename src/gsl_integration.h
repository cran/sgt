/*
This file is part of SGT R Package.

    This package is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this package. If not, see <http://www.gnu.org/licenses/>.

    Author: Carter Davis
    Contact: carterdavis@byu.edu

    Note: Parts of the GNU Scientific Library, version 1.16, were copied and
    modified for this package.
 */

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <R.h>
#include <Rmath.h>

/// Setting Up important Macros - especially to change GSL code to output errors to R

/// These are commented out, only used to test within C only
//#define GSL_ERROR(message, code) printf("GSL ERROR %d: %s", (int)code, message)
//#define GSL_FINITE(x) isfinite(x)

#define GSL_ERROR(message, code) error("GSL ERROR %d: %s", (int)code, message)
#define GSL_FINITE(x) R_FINITE(x)
#define GSL_ERROR_VAL(message, code, value) \
       do { \
       GSL_ERROR(message, code); \
       return value; \
       } while (0)
#define RETURN_IF_NULL(x) if (x == NULL) { return; }
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MAX_DBL(a,b) GSL_MAX(a, b)
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))
#define GSL_MIN_DBL(a,b) GSL_MIN(a,b)
#define GSL_COERCE_DBL(x) (double)x

/// PULLED FROM gsl_machine.h

#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/// PULLED FROM gsl_errno.h

enum {
  GSL_SUCCESS  = 0,
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

/// PULLED FROM gsl_math.h

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))

/// PULLED FROM gsl_integration.h

typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;

gsl_integration_workspace *
  gsl_integration_workspace_alloc (const size_t n);

void
  gsl_integration_workspace_free (gsl_integration_workspace * w);


/* Workspace for QAWS integrator */

typedef struct
{
  double alpha;
  double beta;
  int mu;
  int nu;
  double ri[25];
  double rj[25];
  double rg[25];
  double rh[25];
}
gsl_integration_qaws_table;

gsl_integration_qaws_table *
gsl_integration_qaws_table_alloc (double alpha, double beta, int mu, int nu);

int
gsl_integration_qaws_table_set (gsl_integration_qaws_table * t,
                                double alpha, double beta, int mu, int nu);

void
gsl_integration_qaws_table_free (gsl_integration_qaws_table * t);

/* Workspace for QAWO integrator */

enum gsl_integration_qawo_enum { GSL_INTEG_COSINE, GSL_INTEG_SINE };

typedef struct
{
  size_t n;
  double omega;
  double L;
  double par;
  enum gsl_integration_qawo_enum sine;
  double *chebmo;
}
gsl_integration_qawo_table;

gsl_integration_qawo_table *
gsl_integration_qawo_table_alloc (double omega, double L,
                                  enum gsl_integration_qawo_enum sine,
                                  size_t n);

int
gsl_integration_qawo_table_set (gsl_integration_qawo_table * t,
                                double omega, double L,
                                enum gsl_integration_qawo_enum sine);

int
gsl_integration_qawo_table_set_length (gsl_integration_qawo_table * t,
                                       double L);

void
gsl_integration_qawo_table_free (gsl_integration_qawo_table * t);


/* Definition of an integration rule */

typedef void gsl_integration_rule (const gsl_function * f,
                                   double a, double b,
                                   double *result, double *abserr,
                                   double *defabs, double *resabs);

void gsl_integration_qk15 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk21 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk31 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk41 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk51 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk61 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qcheb (gsl_function * f, double a, double b,
                            double *cheb12, double *cheb24);

/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
  {
    GSL_INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS21 = 2,      /* 21 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS31 = 3,      /* 31 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS41 = 4,      /* 41 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS51 = 5,      /* 51 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS61 = 6       /* 61 point Gauss-Kronrod rule */
  };

void
gsl_integration_qk (const int n, const double xgk[],
                    const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function *f, double a, double b,
                    double * result, double * abserr,
                    double * resabs, double * resasc);


int gsl_integration_qng (const gsl_function * f,
                         double a, double b,
                         double epsabs, double epsrel,
                         double *result, double *abserr,
                         size_t * neval);

int gsl_integration_qag (const gsl_function * f,
                         double a, double b,
                         double epsabs, double epsrel, size_t limit,
                         int key,
                         gsl_integration_workspace * workspace,
                         double *result, double *abserr);

int gsl_integration_qagi (gsl_function * f,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qagiu (gsl_function * f,
                           double a,
                           double epsabs, double epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           double *result, double *abserr);

int gsl_integration_qagil (gsl_function * f,
                           double b,
                           double epsabs, double epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           double *result, double *abserr);


int gsl_integration_qags (const gsl_function * f,
                          double a, double b,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qagp (const gsl_function * f,
                          double *pts, size_t npts,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qawc (gsl_function *f,
                          const double a, const double b, const double c,
                          const double epsabs, const double epsrel, const size_t limit,
                          gsl_integration_workspace * workspace,
                          double * result, double * abserr);

int gsl_integration_qaws (gsl_function * f,
                          const double a, const double b,
                          gsl_integration_qaws_table * t,
                          const double epsabs, const double epsrel,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qawo (gsl_function * f,
                          const double a,
                          const double epsabs, const double epsrel,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          gsl_integration_qawo_table * wf,
                          double *result, double *abserr);

int gsl_integration_qawf (gsl_function * f,
                          const double a,
                          const double epsabs,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          gsl_integration_workspace * cycle_workspace,
                          gsl_integration_qawo_table * wf,
                          double *result, double *abserr);

/* Workspace for fixed-order Gauss-Legendre integration */

typedef struct
  {
    size_t n;         /* number of points */
    double *x;        /* Gauss abscissae/points */
    double *w;        /* Gauss weights for each abscissae */
    int precomputed;  /* high precision abscissae/weights precomputed? */
  }
gsl_integration_glfixed_table;


gsl_integration_glfixed_table *
  gsl_integration_glfixed_table_alloc (size_t n);

void
  gsl_integration_glfixed_table_free (gsl_integration_glfixed_table * t);

/* Routine for fixed-order Gauss-Legendre integration */

double
  gsl_integration_glfixed (const gsl_function *f,
                           double a,
                           double b,
                           const gsl_integration_glfixed_table * t);

/* Routine to retrieve the i-th Gauss-Legendre point and weight from t */

int
  gsl_integration_glfixed_point (double a,
                                 double b,
                                 size_t i,
                                 double *xi,
                                 double *wi,
                                 const gsl_integration_glfixed_table * t);


/* Cquad integration - Pedro Gonnet */

/* Data of a single interval */
typedef struct
{
  double a, b;
  double c[64];
  double fx[33];
  double igral, err;
  int depth, rdepth, ndiv;
} gsl_integration_cquad_ival;


/* The workspace is just a collection of intervals */
typedef struct
{
  size_t size;
  gsl_integration_cquad_ival *ivals;
  size_t *heap;
} gsl_integration_cquad_workspace;

gsl_integration_cquad_workspace *
gsl_integration_cquad_workspace_alloc (const size_t n);

void
gsl_integration_cquad_workspace_free (gsl_integration_cquad_workspace * w);

int
gsl_integration_cquad (const gsl_function * f, double a, double b,
		       double epsabs, double epsrel,
		       gsl_integration_cquad_workspace * ws,
		       double *result, double *abserr, size_t * nevals);

