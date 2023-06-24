#include "mltaln.h"
#include "dp.h"
#include "mtxutl.h"
#include "Kband.h"

# define FFTNS
# define NOTPRINT
# define NOTPRINTLEN
# define DEBUG 0

# ifndef FFTNS //constants
const int nalphabets = 26;
const int nscoredalphabets = 20;
const int max_seq = 100000 + 7;
#endif

void G__align11_traceback(char *s1, char *s2, char *r1, char *r2, int len1, int len2, int band, int *move_vec, double *val_vec)
{
	/* Traceback */
	char *tmpr1, *tmpr2;
	int nowi, nowj, val, tari, tarj, h, i, first1, first2, j;
	double first_score = matrix_query(len1, len2, band, val_vec, len1 + 1, len2 + 1), this_score;
    // determine the start place
	first1 = len1, first2 = len2;
    for(i = len1 - 1; ~i; -- i)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            if((this_score = matrix_query(i, len2, band, val_vec, len1 + 1, len2 + 1) + penalty * (len1 - i)) >= first_score)
            {
                first_score = this_score;
                first1 = i;
                first2 = len2;
            }
        }
        else break; // must outside the band
    for(j = len2 - 1; ~j; -- j)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            if((this_score = matrix_query(len1, j, band, val_vec, len1 + 1, len2 + 1) + penalty * (len2 - j)) >= first_score)
            {
                first_score = this_score;
                first1 = len1;
                first2 = j;
            }
        }
        else break; // must outside the band
    nowi = len1, nowj = len2;
	tmpr1 = r1 - 1, tmpr2 = r2 - 1;
    // first cycle
    h = nowi - first1;
    while(h -- > 0)
    {
        *++ tmpr1 = s1[first1 + h];
		*++ tmpr2 = *newgapstr;
    }
    h = nowj - first2;
    while(h -- > 0)
    {
        *++ tmpr1 = *newgapstr;
        *++ tmpr2 = s2[first2 + h];
    }
    // there is no need to add k in this cycle
    nowi = first1, nowj = first2;
    // center cycles
	for(; ;)
	{
		val = matrix_query_INT(nowi, nowj, band, move_vec, len1 + 1, len2 + 1);
		// calcuate the gaps
		if(val < 0) // column gap
		{
			tari = nowi - 1; tarj = nowj + val;
		}
		else if(val > 0) // row gap
		{
			tari = nowi - val; tarj = nowj - 1;
		}
		else // no gap
		{
			tari = nowi - 1; tarj = nowj - 1;
		}
#if _DEBUG
		reporterr("(nowi, nowj) = (%d, %d) -> (tari, tarj) = (%d, %d), val = %d\n", nowi, nowj, tari, tarj, val);
#endif
        *++ tmpr1 = s1[nowi - 1];
		*++ tmpr2 = s2[nowj - 1];
		// insert it!
		h = nowi - tari - 1;
        while(-- h >= 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = s1[tari + h];
			*++ tmpr2 = *newgapstr;
		}
		h = nowj - tarj - 1;
		while(-- h >= 0)
		{
			*++ tmpr1 = *newgapstr;
			*++ tmpr2 = s2[tarj + h];
		}
		nowi = tari, nowj = tarj;
		// reporterr("r1 = %s, r2 = %s\n", r1, r2);
		if (nowi <= 0 || nowj <= 0) break;
	}
    // last cycle
	while (--nowi >= 0)
	{
		*++ tmpr1 = s1[nowi];
		*++ tmpr2 = *newgapstr;
	}
	while (--nowj >= 0)
	{
		*++ tmpr1 = *newgapstr;
		*++ tmpr2 = s2[nowj];
	}
	*++ tmpr1 = 0;
	*++ tmpr2 = 0;
	//if(strlen(r1) != strlen(r2)) reporterr("%s %s\n", r1, r2);

}

double G__align11_Kband(char *s1, char *s2, int len1, int len2, double penalty, int band, double **ad, char *r1, char *r2, int headgp, int tailgp, double score_before)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
	// fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	//printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
	int i, j, k, ll, rr;
	static TLS int len, minlen;
	/*
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval;
	static TLS double score = 0, tmpval, max_data_i = -inf, tmpval_row, tmpval_col;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec, *max_data_j, *mdjp;
	static TLS unsigned long long Kband_len, delta_len;
	minlen = MIN(len1, len2);
	delta_len = (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0);
	Kband_len = (unsigned long long)(len1 + 1) * (len2 + 1) - delta_len;

#if _DEBUG
    reporterr("Info: Galign will alloc %llu = %llu - %llu(%llu) characters.\n",
               Kband_len + 10, (unsigned long long)(len1 + 1) * (len2 + 1), delta_len,
               (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0));
#endif

	move_vec = AllocateIntVecLarge(Kband_len + 10);
	val_vec = AllocateDoubleVecLarge(Kband_len + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);

	/* Initial of the matrix */
	matrix_set(0, 0, band, 0, val_vec, len1 + 1, len2 + 1);

	if(headgp == 1)
	{
		for(i = 1; i <= len1; ++ i) matrix_set(i, 0, band, penalty * i, val_vec, len1 + 1, len2 + 1);
		for(i = 1; i <= len2; ++ i) matrix_set(0, i, band, penalty * i, val_vec, len1 + 1, len2 + 1);
	}

	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = matrix_query(0, i, band, val_vec, len1 + 1, len2 + 1);
		*mpjp = 0;
	}
	*max_data_j = 0;
	*max_place_j = 0;
	for(i = 1; i <= len1; ++ i)
	{
		max_data_i = matrix_query(i - 1, MAX(i - 1 - band, 0), band, val_vec, len1 + 1, len2 + 1);
        max_place_i = MAX(i - 1 - band, 0);
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = matrix_query(i - 1, 0, band, val_vec, len1 + 1, len2 + 1) + penalty;
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}
		for(j = left_band(len1, len2, i, band); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			if(1 <= j && j <= len2)
			{
				tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1);
				moveval = 0;
				/* Calcuate now the strategy of moving */
                if((tmpval_row = max_data_i) > tmpval)
                {
                    tmpval = tmpval_row;
                    moveval = -(j - max_place_i);
#if _DEBUG
                    reporterr("Set by i, j - 1, j = %d, max_place_i = %d, val = %f\n", j, max_place_i, tmpval);
#endif
                }
                if((tmpval_col = *(max_data_j + j - 1)) > tmpval)
                {
                    tmpval = tmpval_col;
                    moveval = +(i - *(max_place_j + j - 1));
#if _DEBUG
                    reporterr("Set by i - 1, j, i = %d, *(max_place_j + j) = %d, val = %f\n", i, *(max_place_j + j), tmpval);
#endif
				}
                matrix_set(i, j, band, tmpval + protein_score(s1[i - 1], s2[j - 1], ad), val_vec, len1 + 1, len2 + 1);
				matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if _DEBUG
                reporterr("\n(%d, %d) = %f (direction = %d)\nmax place i = %d, max data i = %f\nmax place j = ",
                        i, j, matrix_query(i, j, band, val_vec, len1 + 1, len2 + 1),
                        matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
                for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
                reporterr("\nmax data j  = ");
                for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
                reporterr("\n(%d, %d) = (%c, %c) = %f\n", i, j, s1[i - 1], s2[j - 1], protein_score(s1[i - 1], s2[j - 1], ad));
                reporterr("\n\n");
#endif
                /* UPDATED in 2023/06/01: add penalty calc in max_data_i */
                max_data_i += penalty;

				tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1) + penalty;
				/* Calcuate the next */
                if(tmpval > max_data_i)
                {
                    max_data_i = tmpval;
                    max_place_i = j - 1;
#if _DEBUG
                    reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
#endif
                }

                if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
                {
                    tmpval = matrix_query(i - 1, j, band, val_vec, len1 + 1, len2 + 1);
#if _DEBUG
                    reporterr("determine: (%d, %d) ? %f %f\n", i - 1, j, tmpval, *(max_data_j + j));
#endif
                    if(tmpval > *(max_data_j + j))
                    {
                        *(max_data_j + j) = tmpval;
                        *(max_place_j + j) = i - 1;
#if _DEBUG
                        reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
#endif
                    }
                }
			}
		}

		/* assure that the place appears in max_place_j is in the band */
		for(k = 0, mdjp = max_data_j, mpjp = max_place_j; k < i - band; ++ k, ++ mdjp, ++ mpjp) // clean the outside band date (only need to it?)
		{
			*mdjp = -inf + penalty;
			*mpjp = -__INT_MAX__;
		}
        /* UPDATED in 2023/06/01: add penalty calc in max_data_j */
        for(k = left_band(len1 + 1, len2 + 1, i, band), mdjp = max_data_j + left_band(len1 + 1, len2 + 1, i, band);
            k <= right_band(len1 + 1, len2 + 1, i, band);
            ++ k, ++ mdjp) *mdjp += penalty;
	}
	for(i = 0; i <= len1; ++ i) matrix_set_INT(i, 0, band, i + 1, move_vec, len1 + 1, len2 + 1);
	for(i = 0; i <= len2; ++ i) matrix_set_INT(0, i, band, -(i + 1), move_vec, len1 + 1, len2 + 1);

#if _DEBUG
	reporterr("\n     ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%6d ", j - 1)); reporterr("\n%-6d", 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%-6d", i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%6d", matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("??????");
			}
		}
	reporterr("\n       ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%7d%c ", j - 1, s2[j - 1] ? s2[j - 1] : '0')); reporterr("\n%c%-6d", s1[0], 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%c%-6d", s1[i] ? s1[i] : '0', i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%8.1f", matrix_query(i, j, band, val_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("????????");
			}
		}
	reporterr("\n");
#endif
	score = simple_banded_score(len1, len2, band, val_vec);
	if(score == score_before || score_before == inf)
		G__align11_traceback(s1, s2, r1, r2, len1, len2, band, move_vec, val_vec);
	FreeIntVec(move_vec);
	FreeDoubleVec(val_vec);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
	return score;
}

void MSG__align11_traceback(char *s1, char *s2, char *r1, char *r2, int len1, int len2, int band, int *move_vec, double *last_len1, double *last_len2)
{
	/* Traceback */
	char *tmpr1, *tmpr2;
	int nowi, nowj, val, tari, tarj, h, i, first1, first2, j;
	double first_score = last_len2[len2], this_score;
    // determine the start place
	first1 = len1, first2 = len2;
    for(i = len1 - 1; ~i; -- i)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            if((this_score = last_len1[i] + penalty * (len1 - i)) >= first_score)
            {
                first_score = this_score;
                first1 = i;
                first2 = len2;
            }
        }
        else break; // must outside the band
    for(j = len2 - 1; ~j; -- j)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            if((this_score = last_len2[j] + penalty * (len2 - j)) >= first_score)
            {
                first_score = this_score;
                first1 = len1;
                first2 = j;
            }
        }
        else break; // must outside the band
    nowi = len1, nowj = len2;
	tmpr1 = r1 - 1, tmpr2 = r2 - 1;
    // first cycle
    h = nowi - first1;
    while(h -- > 0)
    {
        *++ tmpr1 = s1[first1 + h];
		*++ tmpr2 = *newgapstr;
    }
    h = nowj - first2;
    while(h -- > 0)
    {
        *++ tmpr1 = *newgapstr;
        *++ tmpr2 = s2[first2 + h];
    }
    // there is no need to add k in this cycle
    nowi = first1, nowj = first2;
    // center cycles
	for(; ;)
	{
		val = matrix_query_INT(nowi, nowj, band, move_vec, len1 + 1, len2 + 1);
		// calcuate the gaps
		if(val < 0) // column gap
		{
			tari = nowi - 1; tarj = nowj + val;
		}
		else if(val > 0) // row gap
		{
			tari = nowi - val; tarj = nowj - 1;
		}
		else // no gap
		{
			tari = nowi - 1; tarj = nowj - 1;
		}
#if _DEBUG
		reporterr("(nowi, nowj) = (%d, %d) -> (tari, tarj) = (%d, %d), val = %d\n", nowi, nowj, tari, tarj, val);
#endif
        *++ tmpr1 = s1[nowi - 1];
		*++ tmpr2 = s2[nowj - 1];
		// insert it!
		h = nowi - tari - 1;
        while(-- h >= 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = s1[tari + h];
			*++ tmpr2 = *newgapstr;
		}
		h = nowj - tarj - 1;
		while(-- h >= 0)
		{
			*++ tmpr1 = *newgapstr;
			*++ tmpr2 = s2[tarj + h];
		}
		nowi = tari, nowj = tarj;
		// reporterr("r1 = %s, r2 = %s\n", r1, r2);
		if (nowi <= 0 || nowj <= 0) break;
	}
    // last cycle
	while (--nowi >= 0)
	{
		*++ tmpr1 = s1[nowi];
		*++ tmpr2 = *newgapstr;
	}
	while (--nowj >= 0)
	{
		*++ tmpr1 = *newgapstr;
		*++ tmpr2 = s2[nowj];
	}
	*++ tmpr1 = 0;
	*++ tmpr2 = 0;
	//if(strlen(r1) != strlen(r2)) reporterr("%s %s\n", r1, r2);

}

double MSGalign11_Kband(char *s1, char *s2, int len1, int len2, double penalty, int band, double **ad, char *r1, char *r2, int headgp, int tailgp, double score_before)
{
	/*
		Every iteration of Kband with memsave mode.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
	// fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	//printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
	int i, j, k, ll, rr;
	static TLS int len, minlen;
	/*
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval;
	static TLS double score = 0, tmpval, max_data_i = -inf, tmpval_row, tmpval_col;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec_before, *val_vec_now, *first_len1, *last_len1, *max_data_j, *mdjp, *val_tmp;
	static TLS unsigned long long Kband_len, delta_len;
	minlen = MIN(len1, len2);
	delta_len = (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0);
	Kband_len = (unsigned long long)(len1 + 1) * (len2 + 1) - delta_len;

#if _DEBUG
    reporterr("Info: Galign will alloc %llu = %llu - %llu(%llu) characters.\n",
               Kband_len + 10, (unsigned long long)(len1 + 1) * (len2 + 1), delta_len,
               (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0));
#endif

	move_vec = AllocateIntVecLarge(Kband_len + 10);
	val_vec_before = AllocateDoubleVec(len2 + 10);
	val_vec_now = AllocateDoubleVec(len2 + 10);
	first_len1 = AllocateDoubleVec(len1 + 10);
	last_len1 = AllocateDoubleVec(len1 + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);

	/* Initial of the matrix */
	val_vec_before[0] = 0;

	if(headgp == 1)
	{
		for(i = 1; i <= len1; ++ i) first_len1[i]     = penalty * i;
		for(i = 1; i <= len2; ++ i) val_vec_before[i] = penalty * i;
	}

	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = val_vec_before[i];
		*mpjp = 0;
	}
	*max_data_j = 0;
	*max_place_j = 0;
	if (insideband(0, len1, band, len1 + 1, len2 + 1)) last_len1[0] = val_vec_before[len2];
	for(i = 1; i <= len1; ++ i)
	{
		max_data_i = val_vec_before[MAX(i - 1 - band, 0)];
        max_place_i = MAX(i - 1 - band, 0);
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = val_vec_before[0] + penalty;
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}
		val_vec_now[0] = first_len1[i];
		for(j = left_band(len1, len2, i, band); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			if(1 <= j && j <= len2)
			{
				tmpval = val_vec_before[j - 1];
				moveval = 0;
				/* Calcuate now the strategy of moving */
                if((tmpval_row = max_data_i) > tmpval)
                {
                    tmpval = tmpval_row;
                    moveval = -(j - max_place_i);
#if _DEBUG
                    reporterr("Set by i, j - 1, j = %d, max_place_i = %d, val = %f\n", j, max_place_i, tmpval);
#endif
                }
                if((tmpval_col = *(max_data_j + j - 1)) > tmpval)
                {
                    tmpval = tmpval_col;
                    moveval = +(i - *(max_place_j + j - 1));
#if _DEBUG
                    reporterr("Set by i - 1, j, i = %d, *(max_place_j + j) = %d, val = %f\n", i, *(max_place_j + j), tmpval);
#endif
				}
                val_vec_now[j] = tmpval + protein_score(s1[i - 1], s2[j - 1], ad);
				matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if _DEBUG
                reporterr("\n(%d, %d) = %f (direction = %d)\nmax place i = %d, max data i = %f\nmax place j = ",
                        i, j, val_vec_now[j], matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
                for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
                reporterr("\nmax data j  = ");
                for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
                reporterr("\n(%d, %d) = (%c, %c) = %f\n", i, j, s1[i - 1], s2[j - 1], protein_score(s1[i - 1], s2[j - 1], ad));
                reporterr("\n\n");
#endif
                /* UPDATED in 2023/06/01: add penalty calc in max_data_i */
                max_data_i += penalty;

				tmpval = val_vec_before[j - 1] + penalty;
				/* Calcuate the next */
                if(tmpval > max_data_i)
                {
                    max_data_i = tmpval;
                    max_place_i = j - 1;
#if _DEBUG
                    reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
#endif
                }

                if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
                {
                    tmpval = val_vec_before[j];
#if _DEBUG
                    reporterr("determine: (%d, %d) ? %f %f\n", i - 1, j, tmpval, *(max_data_j + j));
#endif
                    if(tmpval > *(max_data_j + j))
                    {
                        *(max_data_j + j) = tmpval;
                        *(max_place_j + j) = i - 1;
#if _DEBUG
                        reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
#endif
                    }
                }
			}
		}
		last_len1[i] = val_vec_now[len2];
        for(k = left_band(len1 + 1, len2 + 1, i, band), mdjp = max_data_j + left_band(len1 + 1, len2 + 1, i, band);
            k <= right_band(len1 + 1, len2 + 1, i, band);
            ++ k, ++ mdjp) *mdjp += penalty;
#if _DEBUG
		// print val_before and val_now
		reporterr("val_vec before: \n");
		for(k = 0; k <= len2; ++ k)
		{
			if(insideband(i - 1, k, band, len1 + 1, len2 + 1)) reporterr("%6.0f ", val_vec_before[k]);
			else reporterr("???????? ");
		}
		reporterr("\nval_vec now: \n");
		for(k = 0; k <= len2; ++ k)
		{
			if(insideband(i, k, band, len1 + 1, len2 + 1)) reporterr("%6.0f ", val_vec_now[k]);
			else reporterr("???????? ");
		}
		reporterr("\n");
#endif
		// swap the val_before and val_now
		val_tmp = val_vec_before;
		val_vec_before = val_vec_now;
		val_vec_now = val_tmp;
	}
	for(i = 0; i <= len1; ++ i) matrix_set_INT(i, 0, band, i + 1, move_vec, len1 + 1, len2 + 1);
	for(i = 0; i <= len2; ++ i) matrix_set_INT(0, i, band, -(i + 1), move_vec, len1 + 1, len2 + 1);

#if _DEBUG
	reporterr("\n     ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%6d ", j - 1)); reporterr("\n%-6d", 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%-6d", i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%6d", matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("??????");
			}
		}
	reporterr("\nlast_len1: \n");
	for (i = 0; i <= len1; ++i)
	{
		if (insideband(i, len2, band, len1 + 1, len2 + 1)) reporterr("%6.0f ", last_len1[i]);
		else reporterr("???????? ");
	}
	reporterr("\nlast_len2: \n");
	for (j = 0; j <= len2; ++j)
	{
		if (insideband(len1, j, band, len1 + 1, len2 + 1)) reporterr("%6.0f ", val_vec_before[j]);
		else reporterr("???????? ");
	}
#endif
	score = memsave_banded_score(len1, len2, band, last_len1, val_vec_before);
	if(score == score_before || score_before == inf)
		MSG__align11_traceback(s1, s2, r1, r2, len1, len2, band, move_vec, last_len1, val_vec_before);
	FreeIntVec(move_vec);
	FreeDoubleVec(val_vec_before);
	FreeDoubleVec(val_vec_now);
	FreeDoubleVec(first_len1);
	FreeDoubleVec(last_len1);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
	return score;
}

double pairalign(char *s1, char *s2, int len1, int len2, double penalty, int band, double **ad, char *r1, char *r2, int headgp, int tailgp, double score_before)
{
#if 1
	if (nevermemsave) return G__align11_Kband(s1, s2, len1, len2, penalty, band, ad, r1, r2, headgp, tailgp, score_before);
	else
#endif
	return MSGalign11_Kband(s1, s2, len1, len2, penalty, band, ad, r1, r2, headgp, tailgp, score_before);
	ErrorExit("Error: can not determine algorithm. Program will exit.\n");
}

double G__align11(double** n_dynamicmtx, char** seq1, char** seq2, int alloclen, int headgp, int tailgp)
{
	//puts("In G__align11");
	/*
		Kband Algorithm: align **two** sequences
		Return: Alignment Score. The sequence seq1 & seq2 must be aligned.
	*/
	// Part -1: loop varibles
	int i, j;
# ifndef NOTPRINT
	// Part 0: Print some arguments
	printf("nalphabets = %d\n", nalphabets);
	for (i = 0; i < nalphabets; ++i, puts(""))
		for (j = 0; j < nalphabets; ++j, putchar(' '))
			printf("%f", n_dynamicmtx[i][j]);
	puts(seq1[0]);
	puts(seq2[0]);
	printf("alloclen = %d, headgp = %d, tailgp = %d\n", alloclen, headgp, tailgp);
	// Protein alphabets
	for (i = 0; i < nscoredalphabets; ++i, putchar(' ')) printf("%c", amino[i]); puts("");
	// Gap char
	printf("Newgapstr = %s\n", newgapstr);
	//Gap penalty
	printf("Gap penalty = %d\n", penalty);
# endif

	/* Part 1: Defining the varibles of KBand */
	static TLS double** amino_dynamicmtx = NULL;
	static TLS int length1, length2, mxlength, band, needrerun;
	static TLS char* s1, * s2, * tmp, * res1, * res2, * tmp2;
	static TLS double oldval, val, gappenalty;

	/* Part 2: Allocing the varibles of KBand */
	amino_dynamicmtx = AllocateDoubleMtx(0x100, 0x100);

	/* Part 3: Give initial varibles of the varibles of Part 1 */
	//amino__dynamicmtx: weight of protein sequence
	for (i = 0; i < nalphabets; ++i)
		for (j = 0; j < nalphabets; ++j)
			amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];
	//lenght1, length2: length of sequence
	length1 = strlen(seq1[0]);
	length2 = strlen(seq2[0]);

	if (length1 == 0 || length2 == 0)
	{
		if (length1 == 0 && length2 == 0) return 0.0;
		else if (length1 == 0)
		{
			for (j = 0; j < length2; ++j) seq1[0][j] = *newgapstr;
			seq1[0][length2] = 0;
			return 0.0;
		}
		else //len2 == 0
		{
			for (j = 0; j < length1; ++j) seq2[0][j] = *newgapstr;
			seq2[0][length1] = 0;
			return 0.0;
		}
	}

	mxlength = MAX(length1, length2);
	oldval = -inf;
	s1 = seq1[0];
	s2 = seq2[0];
	band = alignband;
	gappenalty = (double)penalty;

	/* Part 4: Exception */
	// No data
	if (length1 == 0 || length2 == 0)
	{
		FreeDoubleMtx(amino_dynamicmtx);
		return 0.0;
	}

	/* Part 5: KBand algorithm and alloc the matrix */
	res1 = AllocateCharVec((MAX(length1, length2) << 1) + 10);
	res2 = AllocateCharVec((MAX(length1, length2) << 1) + 10);


	if (alignband != NOTSPECIFIED)
	{
		if (length1 < length2)
			val = pairalign(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp, inf);
		else
			val = pairalign(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp, inf);
	}
	else
	{
		band = 10;
		if (length1 > length2)
		{
			// reporterr("Swapped\n");
			if (band > mxlength) val = pairalign(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp, inf);
			while (band <= mxlength)
			{
				// According to the paper, band must be larger than abs(length1 - length2)
				val = pairalign(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp, val);
#if DEBUG
				reporterr("Score = %f, band = %d\n", val, band);
#endif
				if (val == -inf)
				{
					band <<= 1;
					fprintf(stderr, "Warning: can not found DP path in Galign11 Kband. Program will retry.\n");
				}
				else if (val < oldval) { band >>= 1; needrerun = 1; break; }
				else if (val == oldval) { needrerun = 0; break; }
				else
				{
					oldval = val;
					band <<= 1;
					if (band > mxlength) { needrerun = 1; break; }
				}
			}
			if(needrerun) pairalign(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp, inf);
		}
		else
		{
			if (band > mxlength) val = pairalign(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp, inf);
			while (band <= mxlength)
			{
				// According to the paper, band must be larger than abs(length1 - length2)
				val = pairalign(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp, val);
#if DEBUG
				reporterr("Score = %f, band = %d\n", val, band);
#endif
				if (val == -inf)
				{
					band <<= 1;
					fprintf(stderr, "Warning: can not found DP path in Galign11 Kband. Program will retry.\n");
				}
				else if (val < oldval) { band >>= 1; needrerun = 1; break; }
				else if (val == oldval) { needrerun = 0; break; }
				else
				{
					oldval = val;
					band <<= 1;
					if (band > mxlength) { needrerun = 1; break; }
				}
			}
			if(needrerun) pairalign(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp, inf);
		}
	}


#if DEBUG
	printf("After alignment, s1 = %s, s2 = %s (Line 356 in Galign11.c)\n", res1, res2);
#endif

#if 1 // reverse tmp1, tmp2
	for (tmp2 = res1; *tmp2; ++tmp2);
	for (tmp = s1; tmp2 != res1; ++tmp)
	{
		*tmp = *--tmp2;
		//printf("%s\n", s1);
	}
	*tmp = 0;
	for (tmp2 = res2; *tmp2; ++tmp2);
	for (tmp = s2; tmp2 != res2; ++tmp)
	{
		//printf("%d %s\n", *tmp, res2);
		*tmp = *--tmp2;
		//printf("%s\n", s2);
	}
	*tmp = 0;
#else // copy directly
	strncpy(s1, res1, strlen(res1));
	strncpy(s2, res2, strlen(res2));
#endif

#if _DEBUG
	if (strlen(s1) != strlen(res1) || strlen(s2) != strlen(res2))
	{
		reporterr("In G__align11, aligned length = (%d, %d)\n", strlen(s1), strlen(s2));
		if(strlen(s1) != strlen(s2)) reporterr("%s\n%s\n", s1, s2);
		reporterr("(%d, %d)\n", strlen(res1), strlen(res2));
	}
#endif

#if DEBUG
	printf("After process, s1 = %s, s2 = %s (Line 388 in Galign11.c)\n", s1, s2);
#endif

	/* Part 6: Free All varibles and return */
	FreeDoubleMtx(amino_dynamicmtx);
	free(res1);
	free(res2);
	return val;
}
