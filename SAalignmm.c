#include "mltaln.h"
#include "dp.h"
#include "Kband.h"

#define DEBUG 0

double Aalign_traceback(char* r1, char* r2, int len1, int len2, int band, int *move_vec, double *val_vec)
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
        *++ tmpr1 = 'o';
		*++ tmpr2 = *newgapstr;
    }
    h = nowj - first2;
    while(h -- > 0)
    {
        *++ tmpr1 = *newgapstr;
        *++ tmpr2 = 'o';
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
        *++ tmpr1 = 'o';
		*++ tmpr2 = 'o';
		// insert it!
		h = nowi - tari - 1;
        while(-- h >= 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = 'o';
			*++ tmpr2 = *newgapstr;
		}
		h = nowj - tarj - 1;
		while(-- h >= 0)
		{
			*++ tmpr1 = *newgapstr;
			*++ tmpr2 = 'o';
		}
		nowi = tari, nowj = tarj;
		// reporterr("r1 = %s, r2 = %s\n", r1, r2);
		if (nowi <= 0 || nowj <= 0) break;
	}
    // last cycle
	while (--nowi >= 0)
	{
		*++ tmpr1 = 'o';
		*++ tmpr2 = *newgapstr;
	}
	while (--nowj >= 0)
	{
		*++ tmpr1 = *newgapstr;
		*++ tmpr2 = 'o';
	}
	*++ tmpr1 = 0;
	*++ tmpr2 = 0;
	//if(strlen(r1) != strlen(r2)) reporterr("%s %s\n", r1, r2);
}

double Aalign_Kband(int icyc, int jcyc, int len1, int len2, double **cpmx1, double **cpmx2, char *r1, char *r2, int band, double penalty, double score_before)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
	//fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	//printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
	int i, j, h, k;
	static TLS int len, swapped, minlen;
	static TLS char *s;
	swapped = 0;
	if(len1 > len2)
	{
		swapped = 1;
		len = len2; len2 = len1; len1 = len;
		s = r1; r1 = r2; r2 = s;
		minlen = len2;
	}
	else minlen = len1;
	/*
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval;
	static TLS double score = 0, tmpval, max_data_i = -inf, tmpval_row, tmpval_col;
	static TLS char *tmpr1, *tmpr2;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec, *max_data_j, *mdjp, *homo_this_line;

	static TLS unsigned long long Kband_len, delta_len;
	delta_len = (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0);
	Kband_len = (unsigned long long)(len1 + 1) * (len2 + 1) - delta_len;

	move_vec = AllocateIntVecLarge(Kband_len + 10);
	val_vec = AllocateDoubleVecLarge(Kband_len + 10);
	homo_this_line = AllocateDoubleVec(len2 + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);


	//initial varibles
	memset(r1, '?', sizeof(char) * (len1 + len2));
	memset(r2, '?', sizeof(char) * (len1 + len2));

	/* Initial of the matrix */
	matrix_set(0, 0, band, 0, val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len1; ++ i) matrix_set(i, 0, band, penalty * i, val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len2; ++ i) matrix_set(0, i, band, penalty * i, val_vec, len1 + 1, len2 + 1);
	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = matrix_query(0, i, band, val_vec, len1 + 1, len2 + 1);
		*mpjp = 0;
	}
	for(i = 1; i <= len1; ++ i)
	{
		reporterr("%f\n\n", matrix_query(i, 0, band, val_vec, len1 + 1, len2 + 1));
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
		// pre-load the previous line homo matrix
		copy_protein_score(cpmx1, cpmx2, len1, len2, band, homo_this_line, i - 1);
#if _DEBUG
		for (j = 0; j <= len2; ++j)
		{
			if (insideband(i - 1, j, band, len1 + 1, len2 + 1)) reporterr("%f ", homo_this_line[j]);
			else reporterr("???? ");
		}
		reporterr("\n");
#endif
		for(j = MAX(left_band(len1 + 1, len2 + 1, i, band), 1); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1);
			moveval = 0;
			/* Calcuate now the strategy of moving */
			if(insideband(i, j - 1, band, len1 + 1, len2 + 1))
			{
				if((tmpval_row = max_data_i) > tmpval)
				{
					tmpval = tmpval_row;
					moveval = -(j - max_place_i);
#if _DEBUG
					reporterr("Set by i, j - 1, j = %d, max_place_i = %d, val = %f\n", j, max_place_i, tmpval);
#endif
				}

			}
			if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
			{
				if((tmpval_col = *(max_data_j + j - 1)) > tmpval)
				{
					tmpval = tmpval_col;
					moveval = +(i - *(max_place_j + j - 1));
#if _DEBUG
					reporterr("Set by i - 1, j, i = %d, *(max_place_j + j) = %d, val = %f\n", i, *(max_place_j + j), tmpval);
#endif
				}
			}
			matrix_set(i, j, band, tmpval + homo_this_line[j - 1], val_vec, len1 + 1, len2 + 1); // protein score = homo[i - 1][j - 1]
			matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if _DEBUG
			reporterr("\n(%d, %d) = %d\nmax place i = %d, max data i = %f\nmax place j = ", i, j, matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
			for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
			reporterr("\nmax data j  = ");
			for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
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

		/* UPDATED in 2023/06/01: add penalty calc in max_data_j */
        for(k = left_band(len1 + 1, len2 + 1, i, band), mdjp = max_data_j + left_band(len1 + 1, len2 + 1, i, band);
            k <= right_band(len1 + 1, len2 + 1, i, band);
            ++ k, ++ mdjp) *mdjp += penalty;

	}
	for(i = 0; i <= len1; ++ i) matrix_set_INT(i, 0, band, i + 1, move_vec, len1 + 1, len2 + 1);
	for(i = 0; i <= len2; ++ i) matrix_set_INT(0, i, band, -(i + 1), move_vec, len1 + 1, len2 + 1);

#if _DEBUG
	reporterr("\n     ");
	for (j = 0; j < len2 + 1; ++j, reporterr("%6d ", j - 1)); reporterr("\n%-6d", 0);
	for (i = 0; i < len1 + 1; ++i, reporterr("\n%-6d", i))
		for (j = 0; j < len2 + 1; ++j, reporterr(" "))
		{
			if (insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%6d", matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("??????");
			}
		}
	reporterr("\n       ");
	for (j = 0; j < len2 + 1; ++j, reporterr("%7d ", j - 1)); reporterr("\n%-6d", 0);
	for (i = 0; i < len1 + 1; ++i, reporterr("\n%-6d", i))
		for (j = 0; j < len2 + 1; ++j, reporterr(" "))
		{
			if (insideband(i, j, band, len1 + 1, len2 + 1))
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
		Aalign_traceback(r1, r2, len1, len2, band, move_vec, val_vec);
	FreeDoubleVec(val_vec);
	FreeDoubleVec(homo_this_line);
	FreeIntVec(move_vec);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
#if DEBUG
	printf("%s %s\n", r1, r2);
#endif
	if(swapped)
	{
		len = len2; len2 = len1; len1 = len;
		s = r1; r1 = r2; r2 = s;
	}

	return score;

}

double MS_Aalign_traceback(char* r1, char* r2, int len1, int len2, int band, int *move_vec, double *last_len1, double *last_len2)
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
        *++ tmpr1 = 'o';
		*++ tmpr2 = *newgapstr;
    }
    h = nowj - first2;
    while(h -- > 0)
    {
        *++ tmpr1 = *newgapstr;
        *++ tmpr2 = 'o';
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
        *++ tmpr1 = 'o';
		*++ tmpr2 = 'o';
		// insert it!
		h = nowi - tari - 1;
        while(-- h >= 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = 'o';
			*++ tmpr2 = *newgapstr;
		}
		h = nowj - tarj - 1;
		while(-- h >= 0)
		{
			*++ tmpr1 = *newgapstr;
			*++ tmpr2 = 'o';
		}
		nowi = tari, nowj = tarj;
		// reporterr("r1 = %s, r2 = %s\n", r1, r2);
		if (nowi <= 0 || nowj <= 0) break;
	}
    // last cycle
	while (--nowi >= 0)
	{
		*++ tmpr1 = 'o';
		*++ tmpr2 = *newgapstr;
	}
	while (--nowj >= 0)
	{
		*++ tmpr1 = *newgapstr;
		*++ tmpr2 = 'o';
	}
	*++ tmpr1 = 0;
	*++ tmpr2 = 0;
	//if(strlen(r1) != strlen(r2)) reporterr("%s %s\n", r1, r2);
}

double MS_Aalign_Kband(int icyc, int jcyc, int len1, int len2, double **cpmx1, double **cpmx2, char *r1, char *r2, int band, double penalty, double score_before)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
	//fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	//printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
	int i, j, h, k;
	static TLS int len, swapped, minlen;
	static TLS char *s;
	swapped = 0;
	if(len1 > len2)
	{
		swapped = 1;
		len = len2; len2 = len1; len1 = len;
		s = r1; r1 = r2; r2 = s;
		minlen = len2;
	}
	else minlen = len1;
	/*
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval;
	static TLS double score = 0, tmpval, max_data_i = -inf, tmpval_row, tmpval_col;
	static TLS char *tmpr1, *tmpr2;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec_before, *val_vec_now, *first_len1, *last_len1, *max_data_j, *mdjp, *homo_this_line, *val_tmp;

	static TLS unsigned long long Kband_len, delta_len;
	delta_len = (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0);
	Kband_len = (unsigned long long)(len1 + 1) * (len2 + 1) - delta_len;

	move_vec = AllocateIntVecLarge(Kband_len + 10);
	val_vec_before = AllocateDoubleVec(len2 + 10);
	val_vec_now = AllocateDoubleVec(len2 + 10);
	first_len1 = AllocateDoubleVec(len1 + 10);
	last_len1 = AllocateDoubleVec(len1 + 10);
	homo_this_line = AllocateDoubleVec(len2 + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);


	//initial varibles
	memset(r1, '?', sizeof(char) * (len1 + len2));
	memset(r2, '?', sizeof(char) * (len1 + len2));

	/* Initial of the matrix */
	val_vec_before[0] = 0;
	for(i = 1; i <= len1; ++ i) first_len1[i]     = penalty * i;
	for(i = 1; i <= len2; ++ i) val_vec_before[i] = penalty * i;
	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = val_vec_before[i];
		*mpjp = 0;
	}
	for(i = 1; i <= len1; ++ i)
	{
		max_data_i = val_vec_before[MAX(i - 1 - band, 0)];
		max_place_i = MAX(i - 1 - band, 0);
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = val_vec_before[0] + penalty; // matrix_query(i - 1, 0, band, val_vec, len1 + 1, len2 + 1) + penalty;
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}
		// pre-load the previous line homo matrix
		copy_protein_score(cpmx1, cpmx2, len1, len2, band, homo_this_line, i - 1);
#if _DEBUG
		for (j = 0; j <= len2; ++j)
		{
			if (insideband(i - 1, j, band, len1 + 1, len2 + 1)) reporterr("%f ", homo_this_line[j]);
			else reporterr("???? ");
		}
		reporterr("\n");
#endif
		val_vec_now[0] = first_len1[i];
		for(j = MAX(left_band(len1 + 1, len2 + 1, i, band), 1); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			tmpval = val_vec_before[j - 1];
			moveval = 0;
			/* Calcuate now the strategy of moving */
			if(insideband(i, j - 1, band, len1 + 1, len2 + 1))
			{
				if((tmpval_row = max_data_i) > tmpval)
				{
					tmpval = tmpval_row;
					moveval = -(j - max_place_i);
#if _DEBUG
					reporterr("Set by i, j - 1, j = %d, max_place_i = %d, val = %f\n", j, max_place_i, tmpval);
#endif
				}

			}
			if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
			{
				if((tmpval_col = *(max_data_j + j - 1)) > tmpval)
				{
					tmpval = tmpval_col;
					moveval = +(i - *(max_place_j + j - 1));
#if _DEBUG
					reporterr("Set by i - 1, j, i = %d, *(max_place_j + j) = %d, val = %f\n", i, *(max_place_j + j), tmpval);
#endif
				}
			}
			val_vec_now[j] = tmpval + homo_this_line[j - 1]; // protein score = homo[i - 1][j - 1]
			matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if _DEBUG
			reporterr("\n(%d, %d) = %d\nmax place i = %d, max data i = %f\nmax place j = ", i, j, matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
			for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
			reporterr("\nmax data j  = ");
			for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
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
		last_len1[i] = val_vec_now[len2];
		/* UPDATED in 2023/06/01: add penalty calc in max_data_j */
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
	for (j = 0; j < len2 + 1; ++j, reporterr("%6d ", j - 1)); reporterr("\n%-6d", 0);
	for (i = 0; i < len1 + 1; ++i, reporterr("\n%-6d", i))
		for (j = 0; j < len2 + 1; ++j, reporterr(" "))
		{
			if (insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%6d", matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("??????");
			}
		}
	reporterr("\n       ");
#endif

	score = memsave_banded_score(len1, len2, band, last_len1, val_vec_before);
	if(score == score_before || score_before == inf)
		MS_Aalign_traceback(r1, r2, len1, len2, band, move_vec, last_len1, val_vec_before);
	FreeDoubleVec(val_vec_before);
	FreeDoubleVec(val_vec_now);
	FreeDoubleVec(first_len1);
	FreeDoubleVec(last_len1);
	FreeDoubleVec(homo_this_line);
	FreeIntVec(move_vec);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
#if DEBUG
	printf("%s %s\n", r1, r2);
#endif
	if(swapped)
	{
		len = len2; len2 = len1; len1 = len;
		s = r1; r1 = r2; r2 = s;
	}

	return score;

}

double twogroupAalign(int icyc, int jcyc, int len1, int len2, double **cpmx1, double **cpmx2, char *r1, char *r2, int band, double penalty, double score_before)
{
#if 1
	if(nevermemsave) return Aalign_Kband(icyc, jcyc, len1, len2, cpmx1, cpmx2, r1, r2, band, penalty, score_before);
	else
#endif
	return MS_Aalign_Kband(icyc, jcyc, len1, len2, cpmx1, cpmx2, r1, r2, band, penalty, score_before);
}


double Aalign( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen )
{
    /*
        Kband Algorithm: align **two** **Profile** sequences
        Return: Alignment Score. The sequence seq1 & seq2 must be aligned.
    */
    // Part -1: loop varibles
    register int i, j;
    /* Part 1: Defining the varibles of KBand */
    static TLS double **amino_dynamicmtx = NULL;
    static TLS int length1, length2, tmplen, band, swapped = 0, swapp, needrerun = 0, mxlength;
	static TLS char **res1, **res2; // Result array
	static TLS char *gap1, *gap2, *swapgap; // Gap string '-' has gap 'o' not gap
	static TLS int gaplen1, gaplen2; // Gap length without '?'
	static TLS double oldval, val, gappenalty;
	static TLS double **cpmx1, **cpmx2; // Every amino in every place
	char *rev, *rev2;

    /* Part 3: Give initial varibles of the varibles of Part 1 */
	band = alignband;
	gappenalty = (double)penalty * 0.5;
    length1 = strlen(seq1[0]);
    length2 = strlen(seq2[0]);
	mxlength = MAX(length1, length2);
	gaplen1 = 0;
	gaplen2 = 0;
    /* Part 4: Exception */
    // No data
    if(length1 == 0 || length2 == 0) return 0.0;

	/* Part 2: Allocing the varibles that using length information */
	length1 += 10;
	length2 += 10;

	cpmx1 = AllocateDoubleMtx(nalphabets, length1);
	cpmx2 = AllocateDoubleMtx(nalphabets, length2);
	gap1 = AllocateCharVec(length1 + length2);
	gap2 = AllocateCharVec(length1 + length2);
	swapgap = AllocateCharVec(length1 + length2);
	res1 = AllocateCharMtx(icyc, length1 + length2 + 80);
	res2 = AllocateCharMtx(jcyc, length1 + length2 + 80);

	length1 -= 10;
	length2 -= 10;

    /* Part 5: KBand algorithm and alloc the matrix */
	cpmx_calc_new(seq1, cpmx1, eff1, length1, icyc);
	cpmx_calc_new(seq2, cpmx2, eff2, length2, jcyc);
	if(alignband != NOTSPECIFIED)
	{
		if(length1 < length2)
			val = twogroupAalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty, inf);
		else
			val = twogroupAalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty, inf);
	}
	else
	{
	band = 10;
	if(length1 < length2)
	{
		if(band > mxlength) val = twogroupAalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty, inf);
		while(band <= mxlength)
		{
			val = twogroupAalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty, val);
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val < oldval) {band >>= 1; needrerun = 1; break;}
			else if(val == oldval) {needrerun = 0; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
		// if(needrerun) twogroupAalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty);
	}
	else
	{
		if(band > mxlength) val = twogroupAalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty, inf);
		while(band <= mxlength)
		{
			val = twogroupAalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty, val);
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val < oldval) {band >>= 1; needrerun = 1; break;}
			else if(val == oldval) {needrerun = 0; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
		// if(needrerun) twogroupAalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty);
	}
	}

	for(rev2 = gap1; ((*rev2) && *rev2 != '?'); ++ rev2) ++ gaplen1;
	for(rev = swapgap; rev2 != gap1; ++ rev)
	{
		*rev = *--rev2;
		//reporterr("%s\n", swapgap);
	}
	*rev = '\0';
	strncpy(gap1, swapgap, gaplen1 * sizeof(char));
	for(rev2 = gap2; ((*rev2) && *rev2 != '?'); ++ rev2) ++ gaplen2;
	for(rev = swapgap; rev2 != gap2; ++ rev)
	{
		//printf("%d %s\n", *tmp, res2);
		*rev = *--rev2;
		//reporterr("%s\n", swapgap);
	}
	*rev = '\0';
	strncpy(gap2, swapgap, gaplen2 * sizeof(char));

	/* According to gap1 & gap2, insert the gap into the profiles */
	if(strchr(gap1, '-'))
	{
		for(i = 0; i < icyc; ++ i) gapireru(res1[i], seq1[i], gap1);
		for(i = 0; i < icyc; ++ i) strcpy(seq1[i], res1[i]);
	}
	if(strchr(gap2, '-'))
	{
		for(j = 0; j < jcyc; ++ j) gapireru(res2[j], seq2[j], gap2);
		for(j = 0; j < jcyc; ++ j) strcpy(seq2[j], res2[j]);
	}

	FreeDoubleMtx(cpmx1);
	FreeDoubleMtx(cpmx2);
	free(gap1);
	free(gap2);
	free(swapgap);
	FreeCharMtx(res1);
	FreeCharMtx(res2);
	return val;
}
