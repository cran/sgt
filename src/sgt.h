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

struct sgt_parm {
	double mu;
	double sigma;
	double lambda;
	double p;
	double q;
};

struct qsgt_parm {
	double mu;
	double sigma;
	double lambda;
	double p;
	double q;
	double y;
};

typedef double (*distfunc)(double, void *);

inline double sgn(double x);

inline double beta(double x, double y);

inline double lbeta(double x, double y);

double slicesampler(int burnin, double start, double *unif, void *params, distfunc leftinv, distfunc rightinv, distfunc pdf);

inline int sgt_parm_check(struct sgt_parm *params);

inline int qsgt_parm_check(struct qsgt_parm *params);

double dsgt_univar(double x, void * params);

double dsgt_univar_nmc(double x, void * params);

double dlsgt_univar(double x, void * params);

double dlsgt_univar_nmc(double x, void * params);

void dsgt(double *out, int *outn, int *fin, double *x, int *xn, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent, int *takelog);

double dsgt_leftinv(double y, void * params);

double dsgt_rightinv(double y, void * params);

double dsgt_leftinv_nmc(double y, void * params);

double dsgt_rightinv_nmc(double y, void * params);

void psgt(double *out, int *outn, int *fin, double *x, int *xn, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent, int *lowertail, int *takelog);

double qsgt_df(double x, void * params);

double qsgt_df_nmc(double x, void * params);

double qsgt_f(double x, void * params);

double qsgt_f_nmc(double x, void * params);

void qsgt_fdf(double x, void *params, double *y, double *dy);

void qsgt_fdf_nmc(double x, void *params, double *y, double *dy);

void qsgt(double *out, int *outn, int *fin, double *x, int *xn, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent, int *lowertail, int *takelog);

void rsgt(double *out, int *outn, int *fin, int *burnin, double *unif, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent);
