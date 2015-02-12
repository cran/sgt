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

#include "gsl_integration.h"

/// PULLED FROM gsl_types.h

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

/// PULLED FROM gsl_roots.h

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double * root, double x_lower, double x_upper);
    int (*iterate) (void *state, gsl_function * f, double * root, double * x_lower, double * x_upper);
  }
gsl_root_fsolver_type;

typedef struct
  {
    const gsl_root_fsolver_type * type;
    gsl_function * function ;
    double root ;
    double x_lower;
    double x_upper;
    void *state;
  }
gsl_root_fsolver;

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function_fdf * f, double * root);
    int (*iterate) (void *state, gsl_function_fdf * f, double * root);
  }
gsl_root_fdfsolver_type;

typedef struct
  {
    const gsl_root_fdfsolver_type * type;
    gsl_function_fdf * fdf ;
    double root ;
    void *state;
  }
gsl_root_fdfsolver;

gsl_root_fsolver *
gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T);
void gsl_root_fsolver_free (gsl_root_fsolver * s);

int gsl_root_fsolver_set (gsl_root_fsolver * s,
                          gsl_function * f,
                          double x_lower, double x_upper);

int gsl_root_fsolver_iterate (gsl_root_fsolver * s);

const char * gsl_root_fsolver_name (const gsl_root_fsolver * s);
double gsl_root_fsolver_root (const gsl_root_fsolver * s);
double gsl_root_fsolver_x_lower (const gsl_root_fsolver * s);
double gsl_root_fsolver_x_upper (const gsl_root_fsolver * s);


gsl_root_fdfsolver *
gsl_root_fdfsolver_alloc (const gsl_root_fdfsolver_type * T);

int
gsl_root_fdfsolver_set (gsl_root_fdfsolver * s,
                         gsl_function_fdf * fdf, double root);

int
gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s);

void
gsl_root_fdfsolver_free (gsl_root_fdfsolver * s);

const char * gsl_root_fdfsolver_name (const gsl_root_fdfsolver * s);
double gsl_root_fdfsolver_root (const gsl_root_fdfsolver * s);

int
gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel);

int
gsl_root_test_residual (double f, double epsabs);

int
gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel);

GSL_VAR const gsl_root_fsolver_type  * gsl_root_fsolver_bisection;
GSL_VAR const gsl_root_fsolver_type  * gsl_root_fsolver_brent;
GSL_VAR const gsl_root_fsolver_type  * gsl_root_fsolver_falsepos;
GSL_VAR const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_newton;
GSL_VAR const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_secant;
GSL_VAR const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson;
