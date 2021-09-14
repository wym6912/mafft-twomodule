#include "Kband.h"
#include <stdio.h>

#define DEBUG 0

extern int nalphabets, amino[];
#define INF (1e+100)

bool insideband(int i, int j, int band)
{
    /* Judge the coordinate in the band */
	if(i < 0 || j < 0) return 0;
    return -band <= i - j && i - j <= band;
}

int mapping(int x, int y, int n)
{
    /* 
	   Mapping the coordinate -> number
       if the coordinate is not in the matrix, return -1
       -1 means out of range, please **process** it
    */
   if (! (0 <= x && x < n && 0 <= y && y < n)) return -1;
   return x * n + y;
}

double protein_score(char a, char b, double **mtx)
{
	/* calcuate the score of c1 and c2, using BLOSUM62 */
	return mtx[a][b];
}

int __abs__(int x)
{
    return x >= 0 ? x : -x;
}

int __max__(int x, int y)
{
    return x > y ? x : y;
}

void init_allline(int n, int band, int *p)
{
	/* Initialize the matrix mapping: matrix coordinate to sequence
	 * See https://gitee.com/wym6912/grad_learn/tree/master/kband_mapping
	 */
	int i;
	for(i = 0; i < n; ++ i)
		p[i] = (band << 1 | 1) - __max__((band - i), 0) - __max__((band + i - n + 1), 0);
	for(i = 1; i < n; ++ i) p[i] += p[i - 1];
}

int query_place(int x, int y, int band, int *s)
{
	/* Query the place in the array 
	 * Now the complexity of the segment is O(1)
	 */
	if(! insideband(x, y, band)) return -1;
	else 
	{
		if(x == 0) return __abs__(y - x);
		else return s[x - 1] + band + (y - x) - __max__((band - x), 0);
	}
	return -1;
}

int matrix_set(int x, int y, int band, double val, double *mtx, int *place_seq)
{
	/* insert node by matrix */
	if(! insideband(x, y, band)) return -1;
	int q = query_place(x, y, band, place_seq);
#if 0//DEBUG
	fprintf(stderr, "(%d, %d) = %d\n", x, y, q);
	fprintf(stderr, "%lf %lf\n", val, mtx[q]);
#endif
	if(q == -1) return -1;
	mtx[q] = val;
	return 0;
}

int matrix_update(int x, int y, int band, double val, double *mtx, int *place_seq)
{
	/* insert node by matrix
	   Add 
	 */
	if(! insideband(x, y, band)) return -1;
	int q = query_place(x, y, band, place_seq);
#if DEBUG
	fprintf(stderr, "(%d, %d) = %d; ", x, y, q);
	fprintf(stderr, "%lf %lf\n", val, mtx[q]);
#endif
	if(q == -1) return -1;
	mtx[q] += val;
#if DEBUG
	fprintf(stderr, "After update, mtx(%d, %d) = %lf\n", x, y, mtx[q]);
#endif
	return 0;
}

double matrix_query(int x, int y, int band, double *mtx, int *place_seq)
{
	/* query node by matrix */
	int key = query_place(x, y, band, place_seq);
	if(key == -1) return -INF;
	return mtx[key];
}

int matrix_set_INT(int x, int y, int band, int val, int *mtx, int *place_seq)
{
	/* insert node by matrix
	 * type: int
	 */
	if(! insideband(x, y, band)) return -1;
	int q = query_place(x, y, band, place_seq);
#if DEBUG
	fprintf(stderr, "(%d, %d) = %d; ", x, y, q);
	fprintf(stderr, "%d %d\n", val, mtx[q]);
#endif
	if(q == -1) return -1;
	mtx[q] = val;
#if DEBUG
	fprintf(stderr, "After add, mtx(%d, %d) = %d\n", x, y, mtx[q]);
#endif
	return 0;
}

int matrix_query_INT(int x, int y, int band, int *mtx, int *place_seq)
{
	/* query node by matrix 
	   type: int
	 */
	int key = query_place(x, y, band, place_seq);
	if(key == -1) return -__INT_MAX__;
	return mtx[key];
}
