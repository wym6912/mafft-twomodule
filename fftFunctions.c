#include "mltaln.h"

#define SEGMENTSIZE 150
#define TMPTMPTMP 0

#define DEBUG 0

void keika( char *str, int current, int all )
{
	if( current == 0 )
		fprintf( stderr, "%s :         ", str );

		fprintf( stderr, "\b\b\b\b\b\b\b\b" );
		fprintf( stderr, "%3d /%3d", current+1, all+1 );

	if( current+1 == all )
		fprintf( stderr, "\b\b\b\b\b\b\b\bdone.     \n" );
}

double maxItch( double *soukan, int size )
{
	int i;
	double value = 0.0;
	double cand;
	for( i=0; i<size; i++ )
		if( ( cand = soukan[i] ) > value ) value = cand;
	return( value );
}

void calcNaiseki( Fukusosuu *value, Fukusosuu *x, Fukusosuu *y )
{
	value->R =  x->R * y->R + x->I * y->I;
	value->I = -x->R * y->I + x->I * y->R;
}

Fukusosuu *AllocateFukusosuuVec( int l1 )
{
	Fukusosuu *value;
	value = (Fukusosuu *)calloc( l1, sizeof( Fukusosuu ) );
	if( !value )
	{
		fprintf( stderr, "Cannot allocate %d FukusosuuVec\n", l1 );
		return( NULL );
	}
	return( value );
}

Fukusosuu **AllocateFukusosuuMtx( int l1, int l2 )
{
	Fukusosuu **value;
	int j;
//	fprintf( stderr, "allocating %d x %d FukusosuuMtx\n", l1, l2 );
	value = (Fukusosuu **)calloc( l1+1, sizeof( Fukusosuu * ) );
	if( !value )
	{
		fprintf( stderr, "Cannot allocate %d x %d FukusosuuVecMtx\n", l1, l2 );
		exit( 1 );
	}
	for( j=0; j<l1; j++ )
	{
		value[j] = AllocateFukusosuuVec( l2 );
		if( !value[j] )
		{
			fprintf( stderr, "Cannot allocate %d x %d FukusosuuVecMtx\n", l1, l2 );
			exit( 1 );
		}
	}
	value[l1] = NULL;
	return( value );
}

Fukusosuu ***AllocateFukusosuuCub( int l1, int l2, int l3 )
{
	Fukusosuu ***value;
	int i;
	value = calloc( l1+1, sizeof( Fukusosuu ** ) );
	if( !value ) ErrorExit( "Cannot allocate Fukusosuu" );
	for( i=0; i<l1; i++ ) value[i] = AllocateFukusosuuMtx( l2, l3 );
	value[l1] = NULL;
	return( value );
}

void FreeFukusosuuVec( Fukusosuu *vec )
{
	free( (void *)vec );
}

void FreeFukusosuuMtx( Fukusosuu **mtx )
{
	int i;

	for( i=0; mtx[i]; i++ )
		free( (void *)mtx[i] );
	free( (void *)mtx );
}

int getKouho( int *kouho, int nkouho, double *soukan, int nlen2 )
{
	int i, j;
	int nlen4 = nlen2 / 2;
	double max;
	double tmp;
	int ikouho = 0; // by D.Mathog, iinoka?
	for( j=0; j<nkouho; j++ )
	{
		max = -9999.9;
		for( i=0; i<nlen2; i++ )
		{
			if( ( tmp = soukan[i] ) > max )
			{
				ikouho = i;
				max = tmp;
			}
		}
#if 0
		if( max < 0.15 )
		{
			break;
		}
#endif
#if 0
		fprintf( stderr, "Kouho No.%d, pos=%d, score=%f, lag=%d\n", j, ikouho, soukan[ikouho], ikouho-nlen4 );
#endif
		soukan[ikouho] = -9999.9;
		kouho[j] = ( ikouho - nlen4 );
	}
	return( j );
}

void zurasu2( int lag, int    clus1, int    clus2,
                       char  **seq1, char  **seq2,
		 			   char **aseq1, char **aseq2 )
{
	int i;
#if 0
	fprintf( stderr, "### lag = %d\n", lag );
#endif
	if( lag > 0 )
	{
		for( i=0; i<clus1; i++ ) aseq1[i] = seq1[i];
		for( i=0; i<clus2; i++ ) aseq2[i] = seq2[i]+lag;
	}
	else
	{
		for( i=0; i<clus1; i++ ) aseq1[i] = seq1[i]-lag;
		for( i=0; i<clus2; i++ ) aseq2[i] = seq2[i];
	}
}

void zurasu( int lag, int    clus1, int    clus2,
                      char  **seq1, char  **seq2,
					  char **aseq1, char **aseq2 )
{
	int i;
#if DEBUG
	fprintf( stderr, "lag = %d\n", lag );
#endif
	if( lag > 0 )
	{
		for( i=0; i<clus1; i++ ) strcpy( aseq1[i], seq1[i] );
		for( i=0; i<clus2; i++ ) strcpy( aseq2[i], seq2[i]+lag );
	}
	else
	{
		for( i=0; i<clus1; i++ ) strcpy( aseq1[i], seq1[i]-lag );
		for( i=0; i<clus2; i++ ) strcpy( aseq2[i], seq2[i] );
	}
}


int alignableReagion( int    clus1, int    clus2,
					   char  **seq1, char  **seq2,
					   double *eff1, double *eff2,
					   Segment *seg )
{
	int i, j, k;
	int status, starttmp = 0; // by D.Mathog, a gess
	double score;
	int value = 0;
	int len, maxlen;
	int length = 0; // by D.Mathog, a gess
	static TLS double *stra = NULL;
	static TLS int alloclen = 0;
	double totaleff;
	double cumscore;
	static TLS double threshold;
	static TLS double *prf1 = NULL;
	static TLS double *prf2 = NULL;
	static TLS int *hat1 = NULL;
	static TLS int *hat2 = NULL;
	int pre1, pre2;
#if 0
	char **seq1pt;
	char **seq2pt;
	double *eff1pt;
	double *eff2pt;
#endif

#if 0
	fprintf( stderr, "### In alignableRegion, clus1=%d, clus2=%d \n", clus1, clus2 );
	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "eff1[0] = %f\n", eff1[0] );
	fprintf( stderr, "eff2[0] = %f\n", eff2[0] );
#endif

	if( clus1 == 0 )
	{
		if( stra ) FreeDoubleVec( stra ); stra = NULL;
		if( prf1 ) FreeDoubleVec( prf1 ); prf1 = NULL;
		if( prf2 ) FreeDoubleVec( prf2 ); prf2 = NULL;
		if( hat1 ) FreeIntVec( hat1 ); hat1 = NULL;
		if( hat2 ) FreeIntVec( hat2 ); hat2 = NULL;
		alloclen = 0;
		return( 0 );
	}

	if( prf1 == NULL )
	{
		prf1 = AllocateDoubleVec( nalphabets );
		prf2 = AllocateDoubleVec( nalphabets );
		hat1 = AllocateIntVec( nalphabets+1 );
		hat2 = AllocateIntVec( nalphabets+1 );
	}

	len = MIN( strlen( seq1[0] ), strlen( seq2[0] ) );
	maxlen = MAX( strlen( seq1[0] ), strlen( seq2[0] ) ) + fftWinSize;
	if( alloclen < maxlen )
	{
		if( alloclen )
		{
			FreeDoubleVec( stra );
		}
		else
		{
			threshold = (int)fftThreshold / 100.0 * 600.0 * fftWinSize;
		}
		stra = AllocateDoubleVec( maxlen );
		alloclen = maxlen;
	}


	totaleff = 0.0;
	for( i=0; i<clus1; i++ ) for( j=0; j<clus2; j++ ) totaleff += eff1[i] * eff2[j];
	for( i=0; i<len; i++ )
	{
		/* make prfs */
		for( j=0; j<nalphabets; j++ )
		{
			prf1[j] = 0.0;
			prf2[j] = 0.0;
		}
#if 0
		seq1pt = seq1;
		eff1pt = eff1;
		j = clus1;
		while( j-- ) prf1[amino_n[(*seq1pt++)[i]]] += *eff1pt++;
#else
		for( j=0; j<clus1; j++ ) prf1[amino_n[(unsigned char)seq1[j][i]]] += eff1[j];
#endif
		for( j=0; j<clus2; j++ ) prf2[amino_n[(unsigned char)seq2[j][i]]] += eff2[j];

		/* make hats */
		pre1 = pre2 = nalphabets;
		for( j=25; j>=0; j-- )
		{
			if( prf1[j] )
			{
				hat1[pre1] = j;
				pre1 = j;
			}
			if( prf2[j] )
			{
				hat2[pre2] = j;
				pre2 = j;
			}
		}
		hat1[pre1] = -1;
		hat2[pre2] = -1;

		/* make site score */
		stra[i] = 0.0;
		for( k=hat1[nalphabets]; k!=-1; k=hat1[k] )
			for( j=hat2[nalphabets]; j!=-1; j=hat2[j] )
//				stra[i] += n_dis[k][j] * prf1[k] * prf2[j];
				stra[i] += n_disFFT[k][j] * prf1[k] * prf2[j];
		stra[i] /= totaleff;
	}

	(seg+0)->skipForeward = 0;
	(seg+1)->skipBackward = 0;
	status = 0;
	cumscore = 0.0;
	score = 0.0;
	for( j=0; j<fftWinSize; j++ ) score += stra[j];

	// determine the first block can be determined
	if( score > threshold )
	{
		status = 1; // status must be 0
		starttmp = 0;
		length = 0;
		cumscore = score;
	}

	for( i=1; i<len-fftWinSize; i++ )
	{
		score = score - stra[i-1] + stra[i+fftWinSize-1];
#if TMPTMPTMP
		fprintf( stderr, "%d %10.0f   ? %10.0f\n", i, score, threshold );
#endif

		if( score > threshold )
		{
#if 0
			seg->start = i;
			seg->end = i;
			seg->center = ( seg->start + seg->end + fftWinSize ) / 2 ;
			seg->score = score;
			status = 0;
			value++;
#else
			if( !status )
			{
				status = 1;
				starttmp = i;
				length = 0;
				cumscore = score - stra[i + fftWinSize - 1];
			}
			length++;
			cumscore += stra[i + fftWinSize - 1];
#endif
		}
		if( score <= threshold || length > SEGMENTSIZE )
		{
			if( status )
			{
				if( length > fftWinSize )
				{
					seg->start = starttmp;
					seg->end = i;
					seg->center = ( seg->start + seg->end + fftWinSize ) / 2 ;
					seg->score = cumscore;
#if 0
					fprintf( stderr, "%d-%d length = %d, score = %f, value = %d\n", seg->start, seg->end, length, cumscore, value );
#endif
					if( length > SEGMENTSIZE )
					{
						(seg+0)->skipForeward = 1;
						(seg+1)->skipBackward = 1;
					}
					else
					{
						(seg+0)->skipForeward = 0;
						(seg+1)->skipBackward = 0;
					}
					value++;
					seg++;
				}
				length = 0;
				cumscore = 0.0;
				status = 0;
				starttmp = i;
				if( value > MAXSEG - 3 ) ErrorExit( "TOO MANY SEGMENTS!");
			}
		}
	}
	if( status && length > fftWinSize )
	{
		seg->end = i;
		seg->start = starttmp;
		seg->center = ( starttmp + i + fftWinSize ) / 2 ;
		seg->score = cumscore;
#if 0
fprintf( stderr, "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
		value++;
	}
#if TMPTMPTMP
	exit( 0 );
#endif
//	fprintf( stderr, "returning %d\n", value );
	return( value );
}


static int permit( Segment *seg1, Segment *seg2 )
{
	return( 0 );
	if( seg1->end >= seg2->start ) return( 0 );
	if( seg1->pair->end >= seg2->pair->start ) return( 0 );
	else return( 1 );
}

void blockAlign2( int *cut1, int *cut2, Segment **seg1, Segment **seg2, double **ocrossscore, int *ncut )
{
	int i, j, k, shift, cur1, cur2, count, klim;
	static TLS int crossscoresize = 0;
	static TLS int *result1 = NULL;
	static TLS int *result2 = NULL;
	static TLS int *ocut1 = NULL;
	static TLS int *ocut2 = NULL;
	double maximum;
	static TLS double **crossscore = NULL;
	static TLS int **track = NULL;
	static TLS double maxj, maxi;
	static TLS int pointj, pointi;

	if( cut1 == NULL)
	{
		if( result1 )
		{
			if( result1 ) free( result1 ); result1 = NULL;
			if( result2 ) free( result2 ); result2 = NULL;
			if( ocut1 ) free( ocut1 ); ocut1 = NULL;
			if( ocut2 ) free( ocut2 ); ocut2 = NULL;
			if( track ) FreeIntMtx( track ); track = NULL;
	        	if( crossscore ) FreeDoubleMtx( crossscore ); crossscore = NULL;
		}
		crossscoresize = 0;
		return;
	}

	if( result1 == NULL )
	{
		result1 = AllocateIntVec( MAXSEG );
		result2 = AllocateIntVec( MAXSEG );
		ocut1 = AllocateIntVec( MAXSEG );
		ocut2 = AllocateIntVec( MAXSEG );
	}

    if( crossscoresize < *ncut+2 )
    {
        crossscoresize = *ncut+2;
		if( fftkeika ) fprintf( stderr, "allocating crossscore and track, size = %d\n", crossscoresize );
		if( track ) FreeIntMtx( track );
        if( crossscore ) FreeDoubleMtx( crossscore );
		track = AllocateIntMtx( crossscoresize, crossscoresize );
        crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
    }

#if 0
	for( i=0; i<*ncut-2; i++ )
		fprintf( stderr, "%d.start = %d, score = %f\n", i, seg1[i]->start, seg1[i]->score );

	for( i=0; i<*ncut; i++ )
		fprintf( stderr, "i=%d, cut1 = %d, cut2 = %d\n", i, cut1[i], cut2[i] );
	for( i=0; i<*ncut; i++ )
	{
		for( j=0; j<*ncut; j++ )
			fprintf( stderr, "%#4.0f ", ocrossscore[i][j] );
		fprintf( stderr, "\n" );
	}
#endif

	for( i=0; i<*ncut; i++ ) for( j=0; j<*ncut; j++ )  /* mudadanaa */
		crossscore[i][j] = ocrossscore[i][j];
	for( i=0; i<*ncut; i++ )
	{
		ocut1[i] = cut1[i];
		ocut2[i] = cut2[i];
	}

	for( i=1; i<*ncut; i++ )
	{
#if 0
		fprintf( stderr, "### i=%d/%d\n", i,*ncut );
#endif
		for( j=1; j<*ncut; j++ )
		{
			pointi = 0; maxi = 0.0;
			klim = j-2;
			for( k=0; k<klim; k++ )
			{
/*
				fprintf( stderr, "k=%d, i=%d\n", k, i );
*/
				if( !k || k>=*ncut-1 || j>=*ncut-1 ) continue;
				if( crossscore[i-1][k] > maxi )
				{
					pointi = k;
					maxi = crossscore[i-1][k];
				}
			}

			pointj = 0; maxj = 0.0;
			klim = i-2;
			for( k=0; k<klim; k++ )
			{
				if( !k || k>=*ncut-1 || i>=*ncut-1 ) continue;
				if( crossscore[k][j-1] > maxj )
				{
					pointj = k;
					maxj = crossscore[k][j-1];
				}
			}

			maxi += penalty;
			maxj += penalty;

			maximum = crossscore[i-1][j-1];
			track[i][j] = 0;

			if( maximum < maxi )
			{
				maximum = maxi ;
				track[i][j] = j - pointi;
			}

			if( maximum < maxj )
			{
				maximum = maxj ;
				track[i][j] = pointj - i;
			}

			crossscore[i][j] += maximum;
		}
	}
#if 0
	for( i=0; i<*ncut; i++ )
	{
		for( j=0; j<*ncut; j++ )
			fprintf( stderr, "%3d ", track[i][j] );
		fprintf( stderr, "\n" );
	}
#endif


	result1[MAXSEG-1] = *ncut-1;
	result2[MAXSEG-1] = *ncut-1;

	for( i=MAXSEG-1; i>=1; i-- )
	{
		cur1 = result1[i];
		cur2 = result2[i];
		if( cur1 == 0 || cur2 == 0 ) break;
		shift = track[cur1][cur2];
		if( shift == 0 )
		{
			result1[i-1] = cur1 - 1;
			result2[i-1] = cur2 - 1;
			continue;
		}
		else if( shift > 0 )
		{
			result1[i-1] = cur1 - 1;
			result2[i-1] = cur2 - shift;
		}
		else if( shift < 0 )
		{
			result1[i-1] = cur1 + shift;
			result2[i-1] = cur2 - 1;
		}
	}

	count = 0;
	for( j=i; j<MAXSEG; j++ )
	{
		if( ocrossscore[result1[j]][result2[j]] == 0.0 ) continue;

		if( result1[j] == result1[j-1] || result2[j] == result2[j-1] )
			if( ocrossscore[result1[j]][result2[j]] > ocrossscore[result1[j-1]][result2[j-1]] )
				count--;

		cut1[count] = ocut1[result1[j]];
		cut2[count] = ocut2[result2[j]];

		count++;
	}

	*ncut = count;
#if 0
	for( i=0; i<*ncut; i++ )
		fprintf( stderr, "i=%d, cut1 = %d, cut2 = %d\n", i, cut1[i], cut2[i] );
#endif
}
