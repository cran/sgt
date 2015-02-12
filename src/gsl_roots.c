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

#include "gsl_roots.h"

/// PULLED FROM convergence.c

int
gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel)
{
  const double abs_lower = fabs(x_lower) ;
  const double abs_upper = fabs(x_upper) ;

  double min_abs, tolerance;

  if (epsrel < 0.0)
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);

  if (epsabs < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  if (x_lower > x_upper)
    GSL_ERROR ("lower bound larger than upper bound", GSL_EINVAL);

  if ((x_lower > 0.0 && x_upper > 0.0) || (x_lower < 0.0 && x_upper < 0.0))
    {
      min_abs = GSL_MIN_DBL(abs_lower, abs_upper) ;
    }
  else
    {
      min_abs = 0;
    }

  tolerance = epsabs + epsrel * min_abs  ;

  if (fabs(x_upper - x_lower) < tolerance)
    return GSL_SUCCESS;

  return GSL_CONTINUE ;
}

int
gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel)
{
  const double tolerance = epsabs + epsrel * fabs(x1)  ;

  if (epsrel < 0.0)
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);

  if (epsabs < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  if (fabs(x1 - x0) < tolerance || x1 == x0)
    return GSL_SUCCESS;

  return GSL_CONTINUE ;
}

int
gsl_root_test_residual (double f, double epsabs)
{
  if (epsabs < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  if (fabs(f) < epsabs)
    return GSL_SUCCESS;

  return GSL_CONTINUE ;
}

/// PULLED FROM steffenson.c

typedef struct
  {
    double f, df;
    double x;
    double x_1;
    double x_2;
    int count;
  }
steffenson_state_t;

static int steffenson_init (void * vstate, gsl_function_fdf * fdf, double * root);
static int steffenson_iterate (void * vstate, gsl_function_fdf * fdf, double * root);

static int
steffenson_init (void * vstate, gsl_function_fdf * fdf, double * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;

  const double x = *root ;

  state->f = GSL_FN_FDF_EVAL_F (fdf, x);
  state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;

  state->x = x;
  state->x_1 = 0.0;
  state->x_2 = 0.0;

  state->count = 1;

  return GSL_SUCCESS;

}

static int
steffenson_iterate (void * vstate, gsl_function_fdf * fdf, double * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;

  double x_new, f_new, df_new;

  double x_1 = state->x_1 ;
  double x = state->x ;

  if (state->df == 0.0)
    {
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    }

  x_new = x - (state->f / state->df);

  GSL_FN_FDF_EVAL_F_DF(fdf, x_new, &f_new, &df_new);

  state->x_2 = x_1 ;
  state->x_1 = x ;
  state->x = x_new;

  state->f = f_new ;
  state->df = df_new ;

  /// EDITED BY CARTER: PUT IN MACRO GSL_FINITE INSTEAD

  if (!GSL_FINITE (f_new))
    {
      GSL_ERROR ("function value is not finite", GSL_EBADFUNC);
    }

  if (state->count < 3)
    {
      *root = x_new ;
      state->count++ ;
    }
  else
    {
      double u = (x - x_1) ;
      double v = (x_new - 2 * x + x_1);

      if (v == 0)
        *root = x_new;  /* avoid division by zero */
      else
        *root = x_1 - u * u / v ;  /* accelerated value */
    }

  /// EDITED BY CARTER: PUT IN MACRO GSL_FINITE INSTEAD

  if (!GSL_FINITE (df_new))
    {
      GSL_ERROR ("derivative value is not finite", GSL_EBADFUNC);
    }

  return GSL_SUCCESS;
}


static const gsl_root_fdfsolver_type steffenson_type =
{"steffenson",                          /* name */
 sizeof (steffenson_state_t),
 &steffenson_init,
 &steffenson_iterate};

const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson = &steffenson_type;

/// PULLED FROM fdfsolver.c

gsl_root_fdfsolver *
gsl_root_fdfsolver_alloc (const gsl_root_fdfsolver_type * T)
{

  gsl_root_fdfsolver * s = (gsl_root_fdfsolver *) malloc (sizeof (gsl_root_fdfsolver));

  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for root solver struct",
                        GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for root solver state",
                        GSL_ENOMEM, 0);
    };

  s->type = T ;
  s->fdf = NULL;

  return s;
}

int
gsl_root_fdfsolver_set (gsl_root_fdfsolver * s, gsl_function_fdf * f, double root)
{
  s->fdf = f;
  s->root = root;

  return (s->type->set) (s->state, s->fdf, &(s->root));
}

int
gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s)
{
  return (s->type->iterate) (s->state, s->fdf, &(s->root));
}

void
gsl_root_fdfsolver_free (gsl_root_fdfsolver * s)
{
  RETURN_IF_NULL (s);
  free (s->state);
  free (s);
}

const char *
gsl_root_fdfsolver_name (const gsl_root_fdfsolver * s)
{
  return s->type->name;
}

double
gsl_root_fdfsolver_root (const gsl_root_fdfsolver * s)
{
  return s->root;
}
