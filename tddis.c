#include "mltaln.h"
#include "Kband.h"

#define DEBUG 0

void cpmx_calc( char **seq, double **cpmx, double *eff, int lgth, int clus )
{
	int  i, j, k;
	double totaleff = 0.0;

	for( i=0; i<clus; i++ ) totaleff += eff[i];
	for( i=0; i<nalphabets; i++ ) for( j=0; j<lgth; j++ ) cpmx[i][j] = 0.0;
	for( j=0; j<lgth; j++ ) for( k=0; k<clus; k++ )
			cpmx[(int)amino_n[(unsigned char)seq[k][j]]][j] += (double)eff[k] / totaleff;
}


int fastconjuction_noname( int *memlist, char **seq, char **aseq, double *peff, double *eff, char *d, double mineff, double *oritotal  )
{
	int m, k, dln;
	char b[BLEN];
	double total;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 );
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		if( eff[m] < mineff )
			peff[k] = mineff;
		else
			peff[k] = eff[m];

		total += peff[k];
	}
	if( oritotal ) *oritotal = total;
#if 1
	for( m=0; m<k; m++ )
	{
//		fprintf( stderr, "Apr17   peff[%d] = %20.10f\n", m, peff[m] );
		peff[m] /= total;
	}
#endif
	return( k );
}
void cpmx_calc_new( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double *cpmxpt, **cpmxptpt;
	char *seqpt;

	j = nalphabets;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = lgth;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
		// fprintf( stderr, "seqpt = %s, seqlen = %d, lgth = %d\n", seqpt, strlen(seqpt), lgth );
		for( j=0; j<lgth; j++ )
		{
			cpmx[(unsigned char)amino_n[(unsigned char)*seqpt++]][j] += feff;
		}
	}
}
void MScpmx_calc_new( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double **cpmxptpt, *cpmxpt;
	char *seqpt;

	j = lgth;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = nalphabets;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
		cpmxptpt = cpmx;
		j = lgth;
		while( j-- )
			(*cpmxptpt++)[(int)amino_n[(unsigned char)*seqpt++]] += feff;
	}
#if 0
	for( j=0; j<lgth; j++ ) for( i=0; i<nalphabets; i++ ) cpmx[j][i] = 0.0;
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		for( j=0; j<lgth; j++ )
			cpmx[j][(int)amino_n[(int)seq[k][j]]] += feff;
	}
#endif
}

void cpmx_calc_add( char **seq, double **cpmx, double *eff, int lgth, int clus ) // lastmem = newmem; summ eff must be 1.0
{
	double neweff, orieff;
	int newmem, i, j;

	newmem = clus-1;
	neweff = eff[clus-1];
	orieff = 1.0-neweff;
#if 1 // TESTING Feb/1/18:00
	for( j=0; j<lgth; j++ )
	{
		for( i=0;i<nalphabets; i++ ) cpmx[i][j] *= orieff;
		cpmx[(unsigned char)amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
	}
#else // possibly faster?
	for( i=0;i<nalphabets; i++ )
	{
		for( j=0; j<lgth; j++ ) cpmx[i][j] *= orieff;
	}
	for( j=0; j<lgth; j++ ) cpmx[(unsigned char)amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
#endif
}

void copy_protein_score(double **cpmx1, double **cpmx2, int len1, int len2, int band, double *arr, int i)
{
	int j, k, l;
	double *scarr = AllocateDoubleVec(nalphabets + 5), *s, tmp;
	memset(arr, 0, sizeof(arr) * len2);
	// copy only one line !!
	for(s = scarr, j = 0; j < nalphabets; ++ j, ++ s) *s = cpmx1[j][i]; // NOTICE: cpmx[alphabet][place]
	for(j = left_band(len1 + 1, len2 + 1, i, band); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
	{
		for(s = scarr, k = 0; k < nalphabets; ++ k, ++ s)
			if (fabs(*s) > eps)
			{
				for(l = 0; l < nalphabets; ++ l)
				{
					if(fabs(tmp = cpmx2[l][j]) > eps)
					{
#if _DEBUG
						reporterr("(%d, %d) = (%c, %c), val = %f\n", i, j, amino[k], amino[l], n_dis_consweight_multi[k][l]);
#endif
						arr[j] += *s * tmp * n_dis_consweight_multi[k][l];
					}
				}
			}
#if _DEBUG
        reporterr("sum: %f\n", arr[j]);
#endif
	}
	FreeDoubleVec(scarr);
	//fprintf(stderr, "OK protein score, band = %d\n", band);
}

double simple_banded_score(int len1, int len2, int band, double *val_vec)
{
	double score = matrix_query(len1, len2, band, val_vec, len1 + 1, len2 + 1), this_score;
    int i, j;
    for(i = len1 - 1; ~i; -- i)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            if((this_score = matrix_query(i, len2, band, val_vec, len1 + 1, len2 + 1) + penalty * (len1 - i)) > score)
                score = this_score;
        }
        else break; // must outside the band
    for(j = len2 - 1; ~j; -- j)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            if((this_score = matrix_query(len1, j, band, val_vec, len1 + 1, len2 + 1) + penalty * (len2 - j)) > score)
                score = this_score;
        }
        else break; // must outside the band
    return score;
}

double memsave_banded_score(int len1, int len2, int band, double *last_len1, double *last_len2)
{
	double score = last_len2[len2], this_score;
    int i, j;
    for(i = len1 - 1; ~i; -- i)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            if((this_score = last_len1[i] + penalty * (len1 - i)) > score)
                score = this_score;
        }
        else break; // must outside the band
    for(j = len2 - 1; ~j; -- j)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            if((this_score = last_len2[j] + penalty * (len2 - j)) > score)
                score = this_score;
        }
        else break; // must outside the band
    return score;
}

double simple_Salignmm_score(int len1, int len2, int band, double *val_vec, double fpenalty_ex, double* og1, double* og2, double* fg1, double* fg2, double* gf1, double* gf2)
{
	double score = matrix_query(len1, len2, band, val_vec, len1 + 1, len2 + 1), this_score;
    int i, j;
    gf1 += len1 - 1, gf2 += len2 - 1;
    for(i = len1 - 1; i; -- i) // not compare (0, len2)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            this_score = matrix_query(i, len2, band, val_vec, len1 + 1, len2 + 1) + fpenalty_ex * (len1 - i) +
                         *(og1 + i - 1) * *gf2 + *(fg1 + i) * *(gf2 + 1);
            if(this_score > score) score = this_score;
        }
        else break; // must outside the band
    for(j = len2 - 1; j; -- j) // not compare (len1, 0)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            this_score = matrix_query(len1, j, band, val_vec, len1 + 1, len2 + 1) + fpenalty_ex * (len2 - j) +
                         *(og2 + j - 1) * *gf1 + *(fg2 + j) * *(gf1 + 1);
            if(this_score > score) score = this_score;
        }
        else break; // must outside the band
    return score;
}

double memsave_Salignmm_score(int len1, int len2, int band, double *last_len1, double *last_len2, double fpenalty_ex, double* og1, double* og2, double* fg1, double* fg2, double* gf1, double* gf2)
{
    double score = last_len2[len2], this_score;
    int i, j;
    gf1 += len1 - 1, gf2 += len2 - 1;
    for(i = len1 - 1; i; -- i) // not compare (0, len2)
        if(insideband(i, len2, band, len1 + 1, len2 + 1))
        {
            this_score = last_len1[i] + fpenalty_ex * (len1 - i) + *(og1 + i - 1) * *gf2 + *(fg1 + i) * *(gf2 + 1);
            if(this_score > score) score = this_score;
        }
        else break; // must outside the band
    for(j = len2 - 1; j; -- j) // not compare (len1, 0)
        if(insideband(len1, j, band, len1 + 1, len2 + 1))
        {
            this_score = last_len2[j] + fpenalty_ex * (len2 - j) + *(og2 + j - 1) * *gf1 + *(fg2 + j) * *(gf1 + 1);
            if(this_score > score) score = this_score;
        }
        else break; // must outside the band
    return score;
}
