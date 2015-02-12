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

#include "sgt.h"

inline double sgn(double x)
{
	return ((x) >= 0.0 ? 1.0 : -1.0);
}

inline double beta(double x, double y)
{
	return exp(lgamma(x) + lgamma(y) - lgamma(x + y));
}

inline double lbeta(double x, double y)
{
	return (lgamma(x) + lgamma(y) - lgamma(x + y));
}

double slicesampler(int burnin, double start, double *unif, void *params, distfunc leftinv, distfunc rightinv, distfunc pdf)
{
	int i;
	double xval, yval, left, right;
	xval = start;
	for (i = 0; i < burnin; i++) {
		yval = pdf(xval, params);
		yval = (unif[2*i]) * yval;
		left = leftinv(yval, params);
		right = rightinv(yval, params);
		xval = (unif[2*i+1])*(right-left)+left;
	}
	return xval;
}

inline int sgt_parm_check(struct sgt_parm *params)
{
	if (params->sigma > 0 && params->p > 0 && params->q > 0 && params->lambda > -1 && params->lambda < 1 && (params->p)*(params->q) > 1)
		return 1;
	else
		return 0;
}

inline int qsgt_parm_check(struct qsgt_parm *params)
{
	if (params->sigma > 0 && params->p > 0 && params->q > 0 && params->lambda > -1 && params->lambda < 1 && (params->p)*(params->q) > 1 && (params->y) > 1e-12 && (params->y) < (1 - 1e-12))
		return 1;
	else
		return 0;
}

double dsgt_univar(double x, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	double mean = (2*(parm->sigma)*(parm->lambda)*pow((parm->q),(1/(parm->p)))*beta(2/(parm->p),(parm->q)-1/(parm->p)))/beta(1/(parm->p),(parm->q));
	return ((parm->p)/(2*(parm->sigma)*pow((parm->q),(1/(parm->p)))*beta(1/(parm->p),(parm->q))*pow((1+pow(fabs(x-(parm->mu)+mean),(parm->p))/((parm->q)*pow((parm->sigma),(parm->p))*pow((1+(parm->lambda)*sgn(x-(parm->mu)+mean)),(parm->p)))),((parm->q)+1/(parm->p)))));
}

double dsgt_univar_nmc(double x, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	return ((parm->p)/(2*(parm->sigma)*pow((parm->q),(1/(parm->p)))*beta(1/(parm->p),(parm->q))*pow((1+pow(fabs(x-(parm->mu)),(parm->p))/((parm->q)*pow((parm->sigma),(parm->p))*pow((1+(parm->lambda)*sgn(x-(parm->mu))),(parm->p)))),((parm->q)+1/(parm->p)))));
}

double dlsgt_univar(double x, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	double mean = (2*(parm->sigma)*(parm->lambda)*pow((parm->q),(1/(parm->p)))*beta(2/(parm->p),(parm->q)-1/(parm->p)))/beta(1/(parm->p),(parm->q));
	return (log((parm->p))-log(2)-log((parm->sigma))-log((parm->q))/(parm->p)-lbeta(1/(parm->p),(parm->q))-(1/(parm->p)+(parm->q))*log(1+pow(fabs(x-(parm->mu)+mean),(parm->p))/((parm->q)*pow((parm->sigma),(parm->p))*pow((1+(parm->lambda)*sgn(x-(parm->mu)+mean)),(parm->p)))));
}

double dlsgt_univar_nmc(double x, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	return (log((parm->p))-log(2)-log((parm->sigma))-log((parm->q))/(parm->p)-lbeta(1/(parm->p),(parm->q))-(1/(parm->p)+(parm->q))*log(1+pow(fabs(x-(parm->mu)),(parm->p))/((parm->q)*pow((parm->sigma),(parm->p))*pow((1+(parm->lambda)*sgn(x-(parm->mu))),(parm->p)))));
}

void dsgt(double *out, int *outn, int *fin, double *x, int *xn, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent, int *takelog)
{
	int i;
	struct sgt_parm params;
	distfunc pdf;
	if (*meancent && !(*takelog))
		pdf = &dsgt_univar;
	else if (*meancent && *takelog)
		pdf = &dlsgt_univar;
	else if (!(*meancent) && !(*takelog))
		pdf = &dsgt_univar_nmc;
	else
		pdf = &dlsgt_univar_nmc;
	for (i = 0; i < *outn; i++) {
		params.mu = mu[i % *mun];
		params.sigma = sigma[i % *sigman];
		params.lambda = lambda[i % *lambdan];
		params.p = p[i % *pn];
		params.q = q[i % *qn];
		if (fin[i] && sgt_parm_check(&params))
			out[i] = pdf(x[i % *xn], &params);
		else
			out[i] = NAN;
	}
}

double dsgt_leftinv(double y, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	double mean = (2*(parm->sigma)*(parm->lambda)*pow((parm->q),(1/(parm->p)))*beta(2/(parm->p),(parm->q)-1/(parm->p)))/beta(1/(parm->p),(parm->q));
	return (parm->mu)-mean-((pow((parm->q),1/(parm->p))*(parm->sigma)*(1-(parm->lambda)))*pow(pow((parm->p)/(2*y*(parm->sigma)*pow((parm->q),1/(parm->p))*beta(1/(parm->p),(parm->q))),(parm->p)/(1+(parm->q)*(parm->p)))-1,1/(parm->p)));
}

double dsgt_rightinv(double y, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	double mean = (2*(parm->sigma)*(parm->lambda)*pow((parm->q),(1/(parm->p)))*beta(2/(parm->p),(parm->q)-1/(parm->p)))/beta(1/(parm->p),(parm->q));
	return (parm->mu)-mean+((pow((parm->q),1/(parm->p))*(parm->sigma)*(1+(parm->lambda)))*pow(pow((parm->p)/(2*y*(parm->sigma)*pow((parm->q),1/(parm->p))*beta(1/(parm->p),(parm->q))),(parm->p)/(1+(parm->q)*(parm->p)))-1,1/(parm->p)));
}

double dsgt_leftinv_nmc(double y, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	return (parm->mu)-((pow((parm->q),1/(parm->p))*(parm->sigma)*(1-(parm->lambda)))*pow(pow((parm->p)/(2*y*(parm->sigma)*pow((parm->q),1/(parm->p))*beta(1/(parm->p),(parm->q))),(parm->p)/(1+(parm->q)*(parm->p)))-1,1/(parm->p)));
}

double dsgt_rightinv_nmc(double y, void * params)
{
	struct sgt_parm *parm = (struct sgt_parm *)params;
	return (parm->mu)+((pow((parm->q),1/(parm->p))*(parm->sigma)*(1+(parm->lambda)))*pow(pow((parm->p)/(2*y*(parm->sigma)*pow((parm->q),1/(parm->p))*beta(1/(parm->p),(parm->q))),(parm->p)/(1+(parm->q)*(parm->p)))-1,1/(parm->p)));
}

void psgt(double *out, int *outn, int *fin, double *x, int *xn, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent, int *lowertail, int *takelog)
{
	int i;
	size_t ws_size = 500*sizeof(double);
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (ws_size);
	double error, result1, result2;
	gsl_function F;
	if (*meancent)
		F.function = &dsgt_univar;
	else
		F.function = &dsgt_univar_nmc;
	struct sgt_parm params;
	for (i = 0; i < *outn; i++) {
		params.mu = mu[i % *mun];
		params.sigma = sigma[i % *sigman];
		params.lambda = lambda[i % *lambdan];
		params.p = p[i % *pn];
		params.q = q[i % *qn];
		if (fin[i] && sgt_parm_check(&params)) {
			F.params = &params;
			if (x[i % *xn] > 2 * sigma[i % *sigman]) {
				gsl_integration_qagil(&F, params.mu, 0, 1e-10, ws_size, w, &result1, &error);
				gsl_integration_qag (&F, params.mu, x[i % *xn], 0, 1e-10, ws_size, 6, w, &result2, &error);
				out[i] = result1 + result2;
			}
			else
				gsl_integration_qagil(&F, x[i % *xn], 0, 1e-10, ws_size, w, out+i, &error);
			if (!(*lowertail))
				out[i] = 1 - out[i];
			if (*takelog)
				out[i] = log(out[i]);
		}
		else {
			out[i] = NAN;
		}
	}
	gsl_integration_workspace_free (w);
}

double qsgt_df(double x, void * params)
{
	struct qsgt_parm *parm = (struct qsgt_parm *)params;
	double mean = (2*(parm->sigma)*(parm->lambda)*pow((parm->q),(1/(parm->p)))*beta(2/(parm->p),(parm->q)-1/(parm->p)))/beta(1/(parm->p),(parm->q));
	return ((parm->p)/(2*(parm->sigma)*pow((parm->q),(1/(parm->p)))*beta(1/(parm->p),(parm->q))*pow((1+pow(fabs(x-(parm->mu)+mean),(parm->p))/((parm->q)*pow((parm->sigma),(parm->p))*pow((1+(parm->lambda)*sgn(x-(parm->mu)+mean)),(parm->p)))),((parm->q)+1/(parm->p)))));
}

double qsgt_df_nmc(double x, void * params)
{
	struct qsgt_parm *parm = (struct qsgt_parm *)params;
	return ((parm->p)/(2*(parm->sigma)*pow((parm->q),(1/(parm->p)))*beta(1/(parm->p),(parm->q))*pow((1+pow(fabs(x-(parm->mu)),(parm->p))/((parm->q)*pow((parm->sigma),(parm->p))*pow((1+(parm->lambda)*sgn(x-(parm->mu))),(parm->p)))),((parm->q)+1/(parm->p)))));
}

double qsgt_f(double x, void * params)
{
	struct qsgt_parm *parm = (struct qsgt_parm *)params;
	size_t ws_size = 500*sizeof(double);
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (ws_size);
	double result, error, result1, result2;
	gsl_function F;
	F.function = &qsgt_df;
	F.params = params;
	if (x > 2 * (parm->sigma)) {
		gsl_integration_qagil(&F, parm->mu, 0, 1e-10, ws_size, w, &result1, &error);
		gsl_integration_qag (&F, parm->mu, x, 0, 1e-10, ws_size, 6, w, &result2, &error);
		result = result1 + result2;
	}
	else
		gsl_integration_qagil(&F, x, 0, 1e-10, ws_size, w, &result, &error);
	gsl_integration_workspace_free (w);
	return result - (parm->y);
}

double qsgt_f_nmc(double x, void * params)
{
	struct qsgt_parm *parm = (struct qsgt_parm *)params;
	size_t ws_size = 500*sizeof(double);
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (ws_size);
	double result, error, result1, result2;
	gsl_function F;
	F.function = &qsgt_df_nmc;
	F.params = params;
	if (x > 2 * (parm->sigma)) {
		gsl_integration_qagil(&F, parm->mu, 0, 1e-10, ws_size, w, &result1, &error);
		gsl_integration_qag (&F, parm->mu, x, 0, 1e-10, ws_size, 6, w, &result2, &error);
		result = result1 + result2;
	}
	else
		gsl_integration_qagil(&F, x, 0, 1e-10, ws_size, w, &result, &error);
	gsl_integration_workspace_free (w);
	return result - (parm->y);
}

void qsgt_fdf(double x, void *params, double *y, double *dy)
{
	*y = qsgt_f(x, params);
	*dy = qsgt_df(x, params);
}

void qsgt_fdf_nmc(double x, void *params, double *y, double *dy)
{
	*y = qsgt_f_nmc(x, params);
	*dy = qsgt_df_nmc(x, params);
}

void qsgt(double *out, int *outn, int *fin, double *x, int *xn, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent, int *lowertail, int *takelog)
{
	int i;
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
	gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (T);
	double xval, x0, r;
	gsl_function_fdf FDF;
	if (*meancent) {
		FDF.f = &qsgt_f;
		FDF.df = &qsgt_df;
		FDF.fdf = &qsgt_fdf;
	}
	else {
		FDF.f = &qsgt_f_nmc;
		FDF.df = &qsgt_df_nmc;
		FDF.fdf = &qsgt_fdf_nmc;
	}
	struct qsgt_parm params;
	for (i = 0; i < *outn; i++) {
		xval = x[i % *xn];
		if (*takelog)
			xval = exp(xval);
		if (!(*lowertail))
			xval = 1 - xval;
		params.mu = mu[i % *mun];
		params.sigma = sigma[i % *sigman];
		params.lambda = lambda[i % *lambdan];
		params.p = p[i % *pn];
		params.q = q[i % *qn];
		params.y = xval;
		if (fin[i] && qsgt_parm_check(&params)) {
			r = params.mu;
			FDF.params = &params;
			gsl_root_fdfsolver_set (s, &FDF, r);
			iter = 0;
			do {
				iter++;
				status = gsl_root_fdfsolver_iterate (s);
				x0 = r;
				r = gsl_root_fdfsolver_root (s);
				status = gsl_root_test_delta (r, x0, 0, 1e-9);
			} while (status == GSL_CONTINUE && iter < max_iter);
			out[i] = r;
		}
		else {
			out[i] = NAN;
		}
	}
	gsl_root_fdfsolver_free (s);
}

void rsgt(double *out, int *outn, int *fin, int *burnin, double *unif, double *mu, int *mun, double *sigma, int *sigman, double *lambda, int *lambdan, double *p, int *pn, double *q, int *qn, int *meancent)
{
	int i;
	struct sgt_parm params;
	distfunc leftinv;
	distfunc rightinv;
	distfunc pdf;
	if (*meancent) {
		leftinv = &dsgt_leftinv;
		rightinv = &dsgt_rightinv;
		pdf = &dsgt_univar;
	}
	else {
		leftinv = &dsgt_leftinv_nmc;
		rightinv = &dsgt_rightinv_nmc;
		pdf = &dsgt_univar_nmc;
	}
	for (i = 0; i < *outn; i++) {
		params.mu = mu[i % *mun];
		params.sigma = sigma[i % *sigman];
		params.lambda = lambda[i % *lambdan];
		params.p = p[i % *pn];
		params.q = q[i % *qn];
		if (fin[i] && sgt_parm_check(&params))
			out[i] = slicesampler(*burnin, mu[i % *mun], unif+(2*(*burnin)*i), &params, leftinv, rightinv, pdf);
		else
			out[i] = NAN;
	}
}
