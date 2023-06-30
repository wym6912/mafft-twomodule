#include "mltaln.h"
#include "dp.h"
#include "Kband.h"

#define NOTPRINT
#define DEBUG 0
#define NOTPRINTLEN
#define PENALTY_EX 0

double A_align_traceback(char* r1, char* r2, int len1, int len2, int band, int *move_vec, double *val_vec,
                         double fpenalty_ex, double* og1, double* og2, double* fg1, double* fg2, double* gf1, double* gf2)
{
	/* Traceback */
	char *tmpr1, *tmpr2;
	int nowi, nowj, val, tari, tarj, h, i, first1, first2, j;
	double first_score = matrix_query(len1, len2, band, val_vec, len1 + 1, len2 + 1), this_score;
    // determine the start place
	first1 = len1, first2 = len2;
    gf1 += len1 - 1, gf2 += len2 - 1;
    for(i = len1 - 1; i; -- i) // not compare (0, len2)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            this_score = matrix_query(i, len2, band, val_vec, len1 + 1, len2 + 1) + fpenalty_ex * (len1 - i) +
                         *(og1 + i - 1) * *gf2 + *(fg1 + i) * *(gf2 + 1);
            if(this_score >= first_score)
            {
                first_score = this_score;
                first1 = i;
                first2 = len2;
            }
        }
        else break; // must outside the band
    for(j = len2 - 1; j; -- j) // not compare (len1, 0)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            this_score = matrix_query(len1, j, band, val_vec, len1 + 1, len2 + 1) + fpenalty_ex * (len2 - j) +
                         *(og2 + j - 1) * *gf1 + *(fg2 + j) * *(gf1 + 1);
            if(this_score >= first_score)
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

double A_align_Kband(int p1, int p2, int len1, int len2, int band, double fpenalty_ex,
				  double *og1, double *fg1, double *og2, double *fg2, double *gf1, double *gf2, double **cpmx1, double **cpmx2,
				  char *r1, char *r2, double **ad, double hgp1, double hgp2, int headgp, int tailgp, double score_before)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
#if _DEBUG
	reporterr("In Kband, len1 = %d, len2 = %d, band = %d, penalty = %f\n", len1, len2, band, fpenalty_ex);
# endif
	int i, j, k, len = MAX(len1, len2), minlen = MIN(len1, len2);
	/*
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval;
	static TLS double score = 0, tmpval, max_data_i = -inf, tmpval_col, tmpval_row;
	static TLS char *tmpr1, *tmpr2;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec, *homo_this_line, *max_data_j, *mdjp;

	static TLS unsigned long long Kband_len, delta_len;
	delta_len = (unsigned long long)__max__(minlen - 1 - band, 0) * (unsigned long long)__max__(minlen - band, 0);
	Kband_len = (unsigned long long)(len1 + 1) * (len2 + 1) - delta_len;

	move_vec = AllocateIntVecLarge(Kband_len + 10);
	val_vec = AllocateDoubleVecLarge(Kband_len + 10);
	homo_this_line = AllocateDoubleVec(len2 + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);

	// initial varibles
	memset(r1, '?', sizeof(char) * (len1 + len2));
	memset(r2, '?', sizeof(char) * (len1 + len2));
	/* Initial of the matrix */
	matrix_set(0, 0, band, 0, val_vec, len1 + 1, len2 + 1);
	/* headgap may be 0? */
	if(headgp == 1)
	{
		for(i = 1; i <= len1; ++ i) matrix_set(i, 0, band, fpenalty_ex * i + *og1 * hgp2 + *(fg1 + i - 1) * *gf2, val_vec, len1 + 1, len2 + 1);
		for(i = 1; i <= len2; ++ i) matrix_set(0, i, band, fpenalty_ex * i + *og2 * hgp1 + *(fg2 + i - 1) * *gf1, val_vec, len1 + 1, len2 + 1);
	}

	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(j = 0, mdjp = max_data_j, mpjp = max_place_j; j <= len2; ++ mdjp, ++ mpjp, ++ j) // Initial the max_data_j && max_place_j
	{
		*mdjp = matrix_query(0, j, band, val_vec, len1 + 1, len2 + 1);
		*mpjp = 0;
	}
	for(i = 1; i <= len1; ++ i)
	{
#define LBAND(i) left_band(len1 + 1, len2 + 1, i - 1, band)
		max_data_i = matrix_query(i - 1, LBAND(i), band, val_vec, len1 + 1, len2 + 1);
		max_place_i = LBAND(i);
#undef LBAND
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = matrix_query(i - 1, 0, band, val_vec, len1 + 1, len2 + 1);
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}
		// pre-load the previous line homo matrix
		copy_protein_score(cpmx1, cpmx2, len1, len2, band, homo_this_line, i - 1);
		for(j = MAX(left_band(len1 + 1, len2 + 1, i, band), 1); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1);
			moveval = 0;
			/* Calcuate now the strategy of moving */
            if(i != 1 && j != 1)
            {
                // calc row
                tmpval_row = max_data_i + *(og2 + max_place_i) * *(gf1 + i - 2) + *(fg2 + j - 2) * *(gf1 + i - 1);
#if PENALTY_EX
                tmpval_row += fpenalty_ex * (j - 1 - max_place_i);
#endif
                if(tmpval_row >= tmpval)
                {
                    tmpval = tmpval_row;
                    moveval = -(j - max_place_i);
                }
            }
            else if(j != 1)
            {
                // calc row
                tmpval_row = max_data_i + *(og2 + max_place_i) * hgp1 + *(fg2 + j - 2) * *(gf1 + i - 1);
#if PENALTY_EX
                tmpval_row += fpenalty_ex * (j - 1 - max_place_i);
#endif
                if(tmpval_row >= tmpval)
                {
                    tmpval = tmpval_row;
                    moveval = -(j - max_place_i);
                }
            }
            else tmpval_row = -inf;

            if(i != 1 && j != 1)
            {
                // calc col
                mdjp = max_data_j + j - 1;
                mpjp = max_place_j + j - 1;
                tmpval_col = *mdjp + *(og1 + *mpjp) * *(gf2 + j - 2) + *(fg1 + i - 2) * *(gf2 + j - 1);
#if PENALTY_EX
                tmpval_col += fpenalty_ex * (i - 1 - *mpjp);
#endif
				if(tmpval_col >= tmpval)
				{
					tmpval = tmpval_col;
					moveval = +(i - *(max_place_j + j - 1));
				}
			}
            else if(i != 1)
            {
                // calc col
                mdjp = max_data_j + j - 1;
                mpjp = max_place_j + j - 1;
                tmpval_col = *mdjp + *(og1 + *mpjp) * hgp2 + *(fg1 + i - 2) * *(gf2 + j - 1);
#if PENALTY_EX
                tmpval_col += fpenalty_ex * (i - 1 - *mpjp);
#endif
				if(tmpval_col >= tmpval)
				{
					tmpval = tmpval_col;
					moveval = +(i - *(max_place_j + j - 1));
				}
            }
            else tmpval_col = -inf;
#if _DEBUG
            reporterr("(%d, %d) ? diag: %f, row: %f, col: %f -> val = %f, moveval = %d\n", i, j,
                      matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1),
                      tmpval_row, tmpval_col, tmpval, moveval);
#endif

			matrix_set(i, j, band, tmpval + homo_this_line[j - 1], val_vec, len1 + 1, len2 + 1);
			matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if _DEBUG
			reporterr("\n(%d, %d) = %d\nmax place i = %d, max data i = %f\nmax place j = ", i, j, matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
			for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
			reporterr("\nmax data j  = ");
			for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
			reporterr("\n\n");
#endif
			tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1);
			/* Calcuate the next */
            /* There is no need to add fpenalty_ex in these process */
			if(tmpval >= max_data_i)
			{
                max_data_i = tmpval;
                max_place_i = j - 1;
                // reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
            }
            max_data_i += fpenalty_ex;
			if(tmpval >= *(max_data_j + j - 1))
			{
				*(max_data_j + j - 1) = tmpval;
				*(max_place_j + j - 1) = i - 1;
				// reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
			}
            *(max_data_j + j - 1) += fpenalty_ex;
		}
	}
	for(i = 0; i <= len1; ++ i) matrix_set_INT(i, 0, band, i + 1, move_vec, len1 + 1, len2 + 1);
	for(i = 0; i <= len2; ++ i) matrix_set_INT(0, i, band, -(i + 1), move_vec, len1 + 1, len2 + 1);

#if _DEBUG
	reporterr("\n     ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%5d ", j - 1)); reporterr("\n%-5d", 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%-5d", i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%5d", matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("?????");
			}
		}
	reporterr("\n      ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%7d ", j - 1)); reporterr("\n%-6d", 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%-6d", i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%7.1f", matrix_query(i, j, band, val_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("???????");
			}
		}
	reporterr("\n");
#endif

	score = simple_Salignmm_score(len1, len2, band, val_vec, fpenalty_ex, og1, og2, fg1, fg2, gf1, gf2);
	if(score == score_before || score_before == inf)
		A_align_traceback(r1, r2, len1, len2, band, move_vec, val_vec, fpenalty_ex, og1, og2, fg1, fg2, gf1, gf2);

	FreeDoubleVec(val_vec);
	FreeDoubleVec(homo_this_line);
	FreeIntVec(move_vec);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
#if DEBUG
	printf("%s %s\n", r1, r2);
#endif
	return score;
}

double MSA_align_traceback(char* r1, char* r2, int len1, int len2, int band, int *move_vec, double *last_len1, double *last_len2,
                         double fpenalty_ex, double* og1, double* og2, double* fg1, double* fg2, double* gf1, double* gf2)
{
	/* Traceback */
	char *tmpr1, *tmpr2;
	int nowi, nowj, val, tari, tarj, h, i, first1, first2, j;
	double first_score = last_len2[len1], this_score;
    // determine the start place
	first1 = len1, first2 = len2;
    gf1 += len1 - 1, gf2 += len2 - 1;
    for(i = len1 - 1; i; -- i) // not compare (0, len2)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            this_score = last_len1[i] + fpenalty_ex * (len1 - i) + *(og1 + i - 1) * *gf2 + *(fg1 + i) * *(gf2 + 1);
            if(this_score >= first_score)
            {
                first_score = this_score;
                first1 = i;
                first2 = len2;
            }
        }
        else break; // must outside the band
    for(j = len2 - 1; j; -- j) // not compare (len1, 0)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            this_score = last_len2[j] + fpenalty_ex * (len2 - j) + *(og2 + j - 1) * *gf1 + *(fg2 + j) * *(gf1 + 1);
            if(this_score >= first_score)
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

double MS_A_align_Kband(int p1, int p2, int len1, int len2, int band, double fpenalty_ex,
				  double *og1, double *fg1, double *og2, double *fg2, double *gf1, double *gf2, double **cpmx1, double **cpmx2,
				  char *r1, char *r2, double **ad, double hgp1, double hgp2, int headgp, int tailgp, double score_before)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
#if _DEBUG
	reporterr("In Kband, len1 = %d, len2 = %d, band = %d, penalty = %f\n", len1, len2, band, fpenalty_ex);
# endif
	int i, j, k, len = MAX(len1, len2), minlen = MIN(len1, len2);
	/*
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval;
	static TLS double score = 0, tmpval, max_data_i = -inf, tmpval_col, tmpval_row;
	static TLS char *tmpr1, *tmpr2;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec_before, *val_vec_now, *homo_this_line, *first_len1, *last_len1, *max_data_j, *mdjp, *val_tmp;

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

	// initial varibles
	memset(r1, '?', sizeof(char) * (len1 + len2));
	memset(r2, '?', sizeof(char) * (len1 + len2));
	/* Initial of the matrix */
    val_vec_before[0] = 0;
	/* headgap may be 0? */
	if(headgp == 1)
	{
		for(i = 1; i <= len1; ++ i) first_len1[i] = fpenalty_ex * i + *og1 * hgp2 + *(fg1 + i - 1) * *gf2;
		for(i = 1; i <= len2; ++ i) val_vec_before[i] = fpenalty_ex * i + *og2 * hgp1 + *(fg2 + i - 1) * *gf1;
	}

	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(j = 0, mdjp = max_data_j, mpjp = max_place_j; j <= len2; ++ mdjp, ++ mpjp, ++ j) // Initial the max_data_j && max_place_j
	{
		*mdjp = val_vec_before[j];
		*mpjp = 0;
	}
	for(i = 1; i <= len1; ++ i)
	{
#define LBAND(i) left_band(len1 + 1, len2 + 1, i - 1, band)
		max_data_i = val_vec_before[LBAND(i)];
		max_place_i = LBAND(i);
#undef LBAND
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = val_vec_before[0];
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}
        val_vec_now[0] = first_len1[i];
		// pre-load the previous line homo matrix
		copy_protein_score(cpmx1, cpmx2, len1, len2, band, homo_this_line, i - 1);
		for(j = MAX(left_band(len1 + 1, len2 + 1, i, band), 1); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			tmpval = val_vec_before[j - 1];
			moveval = 0;
			/* Calcuate now the strategy of moving */
            if(i != 1 && j != 1)
            {
                // calc row
                tmpval_row = max_data_i + *(og2 + max_place_i) * *(gf1 + i - 2) + *(fg2 + j - 2) * *(gf1 + i - 1);
#if PENALTY_EX
                tmpval_row += fpenalty_ex * (j - 1 - max_place_i);
#endif
                if(tmpval_row >= tmpval)
                {
                    tmpval = tmpval_row;
                    moveval = -(j - max_place_i);
                }
            }
            else if(j != 1)
            {
                // calc row
                tmpval_row = max_data_i + *(og2 + max_place_i) * hgp1 + *(fg2 + j - 2) * *(gf1 + i - 1);
#if PENALTY_EX
                tmpval_row += fpenalty_ex * (j - 1 - max_place_i);
#endif
                if(tmpval_row >= tmpval)
                {
                    tmpval = tmpval_row;
                    moveval = -(j - max_place_i);
                }
            }
            else tmpval_row = -inf;

            if(i != 1 && j != 1)
            {
                // calc col
                mdjp = max_data_j + j - 1;
                mpjp = max_place_j + j - 1;
                tmpval_col = *mdjp + *(og1 + *mpjp) * *(gf2 + j - 2) + *(fg1 + i - 2) * *(gf2 + j - 1);
#if PENALTY_EX
                tmpval_col += fpenalty_ex * (i - 1 - *mpjp);
#endif
				if(tmpval_col >= tmpval)
				{
					tmpval = tmpval_col;
					moveval = +(i - *(max_place_j + j - 1));
				}
			}
            else if(i != 1)
            {
                // calc col
                mdjp = max_data_j + j - 1;
                mpjp = max_place_j + j - 1;
                tmpval_col = *mdjp + *(og1 + *mpjp) * hgp2 + *(fg1 + i - 2) * *(gf2 + j - 1);
#if PENALTY_EX
                tmpval_col += fpenalty_ex * (i - 1 - *mpjp);
#endif
				if(tmpval_col >= tmpval)
				{
					tmpval = tmpval_col;
					moveval = +(i - *(max_place_j + j - 1));
				}
            }
            else tmpval_col = -inf;
#if _DEBUG
            reporterr("(%d, %d) ? diag: %f, row: %f, col: %f -> val = %f, moveval = %d\n", i, j,
                      val_vec_before[j - 1],
                      tmpval_row, tmpval_col, tmpval, moveval);
#endif
            val_vec_now[j] = tmpval + homo_this_line[j - 1];
			matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if _DEBUG
			reporterr("\n(%d, %d) = %d\nmax place i = %d, max data i = %f\nmax place j = ", i, j, matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
			for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
			reporterr("\nmax data j  = ");
			for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
			reporterr("\n\n");
#endif
			tmpval = val_vec_before[j - 1];
			/* Calcuate the next */
            /* There is no need to add fpenalty_ex in these process */
			if(tmpval >= max_data_i)
			{
                max_data_i = tmpval;
                max_place_i = j - 1;
                // reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
            }
            max_data_i += fpenalty_ex;
			if(tmpval >= *(max_data_j + j - 1))
			{
				*(max_data_j + j - 1) = tmpval;
				*(max_place_j + j - 1) = i - 1;
				// reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
			}
            *(max_data_j + j - 1) += fpenalty_ex;
		}
        last_len1[i] = val_vec_now[len2];
        val_tmp = val_vec_before;
        val_vec_before = val_vec_now;
        val_vec_now = val_tmp;
	}
	for(i = 0; i <= len1; ++ i) matrix_set_INT(i, 0, band, i + 1, move_vec, len1 + 1, len2 + 1);
	for(i = 0; i <= len2; ++ i) matrix_set_INT(0, i, band, -(i + 1), move_vec, len1 + 1, len2 + 1);

	score = memsave_Salignmm_score(len1, len2, band, last_len1, val_vec_before, fpenalty_ex, og1, og2, fg1, fg2, gf1, gf2);
	if(score == score_before || score_before == inf)
		MSA_align_traceback(r1, r2, len1, len2, band, move_vec, last_len1, val_vec_before, fpenalty_ex, og1, og2, fg1, fg2, gf1, gf2);

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
	return score;
}

double twogroupAlign(int p1, int p2, int len1, int len2, int band, double fpen, double *og1, double *fg1, double *og2, double *fg2, double *gf1, double *gf2,
					 double **cpmx1, double **cpmx2, char *r1, char *r2, double **ad, double hgp1, double hgp2, int headgp, int tailgp, double score_before)
{
#if 1
	if(nevermemsave) return A_align_Kband(p1, p2, len1, len2, band, fpen, og1, fg1, og2, fg2, gf1, gf2, cpmx1, cpmx2, r1, r2, ad, hgp1, hgp2, headgp, tailgp, score_before);
	else
#endif
	return MS_A_align_Kband(p1, p2, len1, len2, band, fpen, og1, fg1, og2, fg2, gf1, gf2, cpmx1, cpmx2, r1, r2, ad, hgp1, hgp2, headgp, tailgp, score_before);
}


double A__align( double **n_dynamicmtx, int penalty_l, int penalty_ex_l, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, char *sgap1, char *sgap2, char *egap1, char *egap2, int headgp, int tailgp)
{
	//puts("In A__align11");
	/*
		Kband Algorithm: align **two** **Profile** sequences
		Return: Alignment Score. The sequence seq1 & seq2 must be aligned.
	*/
	// Part -1: loop varibles
	int i, j;
# ifndef NOTPRINT
	// Part 0: Print some arguments
	printf("nalphabets = %d, nscoredalphabets = %d\n", nalphabets, nscoredalphabets);
	for(i = 0; i < nscoredalphabets; ++ i, puts(""))
		for(j = 0; j < nscoredalphabets; ++ j, putchar(' '))
			printf("%f", n_dynamicmtx[i][j]);
	for(i = 0; i < icyc; ++ i) puts(seq1[i]);
	puts("");
	for(j = 0; j < jcyc; ++ j) puts(seq2[j]);
	// the number of Sequences
	printf("icyc = %d, jcyc = %d\n", icyc, jcyc);
	printf("alloclen = %d, headgp = %d, tailgp = %d\n", alloclen, headgp, tailgp);
	// Protein alphabets
	for(i = 0; i < nscoredalphabets; ++ i, putchar(' ')) printf("%c", amino[i]); puts("");
	// Gap char
	printf("Newgapstr = %s\n", newgapstr);
	// Gap penalty
	printf("Gap penalty = %d, penalty_ex = %d\n", penalty_l, penalty_ex_l);
	// Start gap
	printf("Start gap pointer = %p, %p\n", sgap1, sgap2);
	// End gap
	printf("End gap pointer = %p, %p\n", egap1, egap2);
# endif

	/* Part 1: Defining the varibles of KBand */
	static TLS double **amino_dynamicmtx = NULL;
	static TLS int length1, length2, mxlength, tmplen, band, swapped = 0, swapp, needrerun;
	static TLS char **res1, **res2; // Result array
	static TLS char *gap1, *gap2, *swapgap; // Gap string '-' has gap 'o' not gap
	static TLS int gaplen1, gaplen2; // Gap length without '?'
	static TLS double oldval, val, gappenalty, fpenalty_ex;
	static TLS double **cpmx1, **cpmx2; // Every amino in every place
	static TLS double *gapfreq1, *gapfreq2; // Gap frequency
	double *gapf1qp, *gapf2qp;
	char *rev, *rev2;

	static TLS double *opgap1f, *opgap2f; // Opening Gap frequency
	static TLS double *fgap1f, *fgap2f; // Closing gap frequency
	double *opgap1fp, *opgap2fp, *fgap1fp, *fgap2fp;
	static TLS double headgapfreq1, headgapfreq2; // first place penalty

	/* Part 2: Allocing the varible of KBand - amino matrix */
	amino_dynamicmtx = AllocateDoubleMtx(0x100, 0x100);

	/* Part 3: Give initial varibles of the varibles of Part 1 */
	//amino__dynamicmtx: weight of protein sequence
	for(i = 0; i < nalphabets; ++ i)
		for(j = 0; j < nalphabets; ++ j)
			amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];
	//lenght1, length2: length of sequence
	length1 = strlen(seq1[0]);
	length2 = strlen(seq2[0]);

	if (length1 == 0 || length2 == 0)
	{
		if (length1 == 0 && length2 == 0) return 0.0;
		else if (length1 == 0)
		{
			for(i = 0; i < icyc; ++ i)
			{
				for (j = 0; j < length2; ++j) seq1[i][j] = *newgapstr;
				seq1[i][length2] = 0;
			}
			return 0.0;
		}
		else //len2 == 0
		{
			for(i = 0; i < jcyc; ++ i)
			{
				for (j = 0; j < length1; ++j) seq2[i][j] = *newgapstr;
				seq2[i][length1] = 0;
			}
			return 0.0;
		}
	}

	mxlength = MAX(length1, length2);
	oldval = -inf;
	band = alignband;
	gappenalty = (double)penalty;
	fpenalty_ex = (double)penalty_ex;

	/* Part 4: Exception */
	// No data
	if(length1 == 0 || length2 == 0)
	{
		FreeDoubleMtx(amino_dynamicmtx);
		return 0.0;
	}

	/* Part 2.2: Allocing the varibles that using length information */
	length1 += 10;
	length2 += 10;

	cpmx1 = AllocateDoubleMtx(nalphabets, length1);
	cpmx2 = AllocateDoubleMtx(nalphabets, length2);
	gapfreq1 = AllocateDoubleVec(length1);gapf1qp = gapfreq1;
	gapfreq2 = AllocateDoubleVec(length2);gapf2qp = gapfreq2;

	opgap1f = AllocateDoubleVec(length1); opgap1fp = opgap1f;
	opgap2f = AllocateDoubleVec(length2); opgap2fp = opgap2f;
	fgap1f = AllocateDoubleVec(length1);  fgap1fp = fgap1f;
	fgap2f = AllocateDoubleVec(length2);  fgap2fp = fgap2f;

	res1 = AllocateCharMtx(icyc, length1 + length2 + 80);
	res2 = AllocateCharMtx(jcyc, length1 + length2 + 80);
	gap1 = AllocateCharVec(length1 + length2);
	gap2 = AllocateCharVec(length1 + length2);
	swapgap = AllocateCharVec(length1 + length2);

	length1 -= 10;
	length2 -= 10;

	gaplen1 = 0;
	gaplen2 = 0;
	/* Part 5: KBand algorithm and alloc the matrix */

	// Calcuate the frequency of all characters
	cpmx_calc_new(seq1, cpmx1, eff1, length1, icyc);
#ifndef NOTPRINTLEN
	fprintf(stderr, "OK on cpmx_calc to seq1\n");
#endif
	cpmx_calc_new(seq2, cpmx2, eff2, length2, jcyc);
#ifndef NOTPRINTLEN
	fprintf(stderr, "OK on cpmx_calc to seq2\n");
#endif
	// Count open gap && final gap
	st_OpeningGapCount(opgap1fp, icyc, seq1, eff1, length1);
	st_FinalGapCount(fgap1fp, icyc, seq1, eff1, length1);

	st_OpeningGapCount(opgap2fp, jcyc, seq2, eff2, length2);
	st_FinalGapCount(fgap2fp, jcyc, seq2, eff2, length2);

	// get head gap frequency && tail gap frequency
	outgapcount(&headgapfreq1, icyc, sgap1, eff1);
	outgapcount(&headgapfreq2, jcyc, sgap2, eff2);

	outgapcount( gapf1qp + length1, icyc, egap1, eff1 );
	outgapcount( gapf2qp + length2, jcyc, egap2, eff2 );

	if( legacygapcost == 0 )
	{
		// Count common gap frequency
		gapcountf(gapf1qp, seq1, icyc, eff1, length1);
		gapcountf(gapf2qp, seq2, jcyc, eff2, length2);

		// Reverse
		for(i = 0; i <= length1; ++ i) gapf1qp[i] = 1.0 - gapf1qp[i];
		for(i = 0; i <= length2; ++ i) gapf2qp[i] = 1.0 - gapf2qp[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for(i = 0; i <= length1; ++ i) gapf1qp[i] = 1.0;
		for(i = 0; i <= length2; ++ i) gapf2qp[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

	// Calcuate the opening gap && final gap
	for(i = 0; i <= length1; ++ i )
	{
		opgap1f[i] = 0.5 * ( 1.0 - opgap1f[i] ) * gappenalty * ( gapf1qp[i] );
		fgap1f[i] = 0.5 * ( 1.0 - fgap1f[i] ) * gappenalty * ( gapf1qp[i] );
	}
	for(i = 0; i <= length2; ++ i )
	{
		opgap2f[i] = 0.5 * ( 1.0 - opgap2f[i] ) * gappenalty * ( gapf2qp[i] );
		fgap2f[i] = 0.5 * ( 1.0 - fgap2f[i] ) * gappenalty * ( gapf2qp[i] );
	}


#if _DEBUG
	reporterr("\nopening gap1: \n");
	for(i = 0; i <= length1; ++ i, reporterr(" ")) reporterr("%f", *(opgap1f + i));
	reporterr("\nopening gap2: \n");
	for(i = 0; i <= length2; ++ i, reporterr(" ")) reporterr("%f", *(opgap2f + i));
	reporterr("\nfinal gap1: \n");
	for(i = 0; i <= length1; ++ i, reporterr(" ")) reporterr("%f", *(fgap1f + i));
	reporterr("\nfinal gap2: \n");
	for(i = 0; i <= length2; ++ i, reporterr(" ")) reporterr("%f", *(fgap2f + i));
	reporterr("\ngap1 penalty: \n");
	for(i = 0; i <= length1; ++ i, reporterr(" ")) reporterr("%f", *(gapf1qp + i));
	reporterr("\ngap2 penalty: \n");
	for(i = 0; i <= length2; ++ i, reporterr(" ")) reporterr("%f", *(gapf2qp + i));
	reporterr("\nhead gap penalty: %f %f\n", headgapfreq1, headgapfreq2);
#endif
	if(alignband != NOTSPECIFIED)
	{
		if(length1 < length2)
			val = twogroupAlign(icyc, jcyc, length1, length2, band, fpenalty_ex, opgap1f, fgap1f, opgap2f, fgap2f, gapf1qp, gapf2qp,
							    cpmx1, cpmx2, gap1, gap2, amino_dynamicmtx, headgapfreq1, headgapfreq2, headgp, tailgp, inf);
		else
			val = twogroupAlign(jcyc, icyc, length2, length1, band, fpenalty_ex, opgap2f, fgap2f, opgap1f, fgap1f, gapf2qp, gapf1qp,
							    cpmx2, cpmx1, gap2, gap1, amino_dynamicmtx, headgapfreq2, headgapfreq1, headgp, tailgp, inf);

	}
	else
	{
	band = 10;
	if(length1 < length2)
	{
		val = twogroupAlign(icyc, jcyc, length1, length2, band, fpenalty_ex, opgap1f, fgap1f, opgap2f, fgap2f, gapf1qp, gapf2qp,
						    cpmx1, cpmx2, gap1, gap2, amino_dynamicmtx, headgapfreq1, headgapfreq2, headgp, tailgp, inf);
		while(band <= mxlength)
		{
			// According to the paper, band must be larger than abs(length1 - length2)
			val = twogroupAlign(icyc, jcyc, length1, length2, band, fpenalty_ex, opgap1f, fgap1f, opgap2f, fgap2f, gapf1qp, gapf2qp,
							    cpmx1, cpmx2, gap1, gap2, amino_dynamicmtx, headgapfreq1, headgapfreq2, headgp, tailgp, val);
#if DEBUG
			reporterr("Val = %f, band = %d, length1 = %d, length2 = %d\n", val, band, length1, length2);
#endif
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val == oldval) {band >>= 1; needrerun = 1; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
	}
	else
	{
		val = twogroupAlign(jcyc, icyc, length2, length1, band, fpenalty_ex, opgap2f, fgap2f, opgap1f, fgap1f, gapf2qp, gapf1qp,
						    cpmx2, cpmx1, gap2, gap1, amino_dynamicmtx, headgapfreq2, headgapfreq1, headgp, tailgp, inf);
		while(band <= mxlength)
		{
			// According to the paper, band must be larger than abs(length1 - length2)
			val = twogroupAlign(jcyc, icyc, length2, length1, band, fpenalty_ex, opgap2f, fgap2f, opgap1f, fgap1f, gapf2qp, gapf1qp,
							    cpmx2, cpmx1, gap2, gap1, amino_dynamicmtx, headgapfreq2, headgapfreq1, headgp, tailgp, val);
#if DEBUG
			reporterr("Val = %f, band = %d, length1 = %d, length2 = %d\n", val, band, length1, length2);
#endif
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val == oldval) {band >>= 1; needrerun = 1; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
	}
	}

#if DEBUG
	reporterr("%s\n", gap1);
#endif
	for(rev2 = gap1; ((*rev2) && *rev2 != '?'); ++ rev2) ++ gaplen1;
	for(rev = swapgap; rev2 != gap1; ++ rev)
	{
		*rev = *--rev2;
		//reporterr("%s\n", swapgap);
	}
	*rev = '\0';
#if DEBUG
	reporterr("gaplen = %d\n", gaplen1);
#endif
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


	/* Part 6: Free All varibles and return */
	FreeDoubleMtx(amino_dynamicmtx);
	FreeCharMtx(res1);
	FreeCharMtx(res2);
	FreeDoubleMtx(cpmx1);
	FreeDoubleMtx(cpmx2);
	free(swapgap);
	free(gap1);
	free(gap2);
	FreeDoubleVec(gapfreq1);
	FreeDoubleVec(gapfreq2);
#if 1
	FreeDoubleVec(fgap1f);
	FreeDoubleVec(fgap2f);
	FreeDoubleVec(opgap1f);
	FreeDoubleVec(opgap2f);
#endif
#ifndef NOTPRINTLEN
	reporterr("In A__align, [");
	for(i = 0; i < icyc; ++ i) reporterr("%d ", strlen(seq1[i]));
	reporterr(",");
	for(j = 0; j < jcyc; ++ j) reporterr("%d ", strlen(seq2[j]));
	reporterr("]\n");
#endif
	return val;
}
