/*
 * Adapted from ncvreg v3.13.0 by Patrick Breheny
 * Original source: https://github.com/pbreheny/ncvreg
 * Licensed under GPL-3.
 *
 * Modifications:
 *   - Renamed R_init_ncvreg -> R_init_PAWPH2
 *   - Removed registration of C routines not used by PAWPH2
 *     (cdfit_binomial, cdfit_gaussian, cdfit_poisson, rawfit_gaussian,
 *      mfdr_binomial, mfdr_cox, mfdr_gaussian)
 *   - Retained: cdfit_cox_dh, maxprod, standardize, and shared helper
 *     functions (crossprod, wcrossprod, wsqsum, sqsum, sum, MCP, SCAD, lasso)
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>

extern SEXP cdfit_cox_dh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP maxprod(SEXP, SEXP, SEXP, SEXP);
extern SEXP standardize(SEXP);

// List accessor
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

// Cross product of y with jth column of X
double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

// Weighted cross product of y with jth column of X
double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(double *X, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

// Sum of squares of jth column of X
double sqsum(double *X, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += pow(X[nn+i], 2);
  return(val);
}

double sum(double *x, int n) {
  double val=0;
  for (int i=0;i<n;i++) val += x[i];
  return(val);
}

double MCP(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(v*(1+l2-1/gamma)));
  else return(z/(v*(1+l2)));
}

double SCAD(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(v*(1+l2)));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)));
  else return(z/(v*(1+l2)));
}

double lasso(double z, double l1, double l2, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else return(s*(fabs(z)-l1)/(v*(1+l2)));
}

static const R_CallMethodDef CallEntries[] = {
  {"cdfit_cox_dh", (DL_FUNC) &cdfit_cox_dh, 13},
  {"maxprod",      (DL_FUNC) &maxprod,        4},
  {"standardize",  (DL_FUNC) &standardize,    1},
  {NULL, NULL, 0}
};

void R_init_PAWPH2(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
