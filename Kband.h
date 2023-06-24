#ifndef __KBAND__
#define __KBAND__

#include <stddef.h>
#include "mltaln.h"
#if _DEBUG
#include <assert.h>
#endif

typedef short int bool;

long long __min__(long long x, long long y);
long long __max__(long long x, long long y);
long long __abs__(long long x);
long long right_band(int n, int m, int i, int band);
long long left_band(int n, int m, int i, int band);
bool insideband(int x, int y, int band, int n, int m);
double protein_score(char a, char b, double **mtx);
size_t query_place(int x, int y, int band, int n, int m);
int matrix_set(int x, int y, int band, double val, double *mtx, int n, int m);
double matrix_query(int x, int y, int band, double *mtx, int n, int m);
int matrix_update(int x, int y, int band, double val, double *mtx, int n, int m);
int matrix_set_INT(int x, int y, int band, int val, int *mtx, int n, int m);
int matrix_query_INT(int x, int y, int band, int *mtx, int n, int m);

#define inf 1e100
#define eps 1e-5

double simple_banded_score(int len1, int len2, int band, double *val_vec);
double memsave_banded_score(int len1, int len2, int band, double *last_len1, double *last_len2);

double simple_Salignmm_score(int len1, int len2, int band, double *val_vec, double fpenalty_ex, double* og1, double* og2, double* fg1, double* fg2, double* gf1, double* gf2);
double memsave_Salignmm_score(int len1, int len2, int band, double *last_len1, double *last_len2, double fpenalty_ex, double* og1, double* og2, double* fg1, double* fg2, double* gf1, double* gf2);

#endif
