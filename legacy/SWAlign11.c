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

#define inf 1e100
static const double eps = 1e-5;

static void match_calc_mtx(double** mtx, double* match, char** s1, char** s2, int i1, int lgth2)
{
	char* seq2 = s2[0];
	double* doubleptr = mtx[(unsigned char)s1[0][i1]];

	while (lgth2--)
		*match++ = doubleptr[(unsigned char)*seq2++];
}


static local_align_pair Atracking__local_no_align(double* lasthorizontalw, double* lastverticalw, int lgth1, int lgth2, int** ijp)
{
	static TLS int i, j, l, iin, jin, ifi, jfi, k, limk;
	static TLS local_align_pair pair;
	pair.start = pair.end = 0;
	char* gap = newgapstr;
	static TLS double wm, g;
	double fpenalty = (double)penalty;
	double fpenalty_ex = (double)penalty_ex;


#if 0
	for (i = 0; i < lgth1; i++)
	{
		fprintf(stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i]);
	}
#endif

	for (i = 0; i < lgth1 + 1; i++)
	{
		ijp[i][0] = i + 1;
	}
	for (j = 0; j < lgth2 + 1; j++)
	{
		ijp[0][j] = -(j + 1);
	}


	wm = lasthorizontalw[lgth2 - 1] - 1.0; // lasthorizontalw[lgth2-1] yori kanarazu chiisai.
	pair.end = lgth2 - 1; // initial value for shorter sequence

	// find near of center to split
	for (j = lgth2 - 2; j > (lgth2 - 2) >> 1; j--)
	{
		if ((g = lasthorizontalw[j]) >= wm)
		{
			wm = g;
			iin = lgth1 - 1; jin = j;
			ijp[lgth1][lgth2] = -(lgth2 - j);
			pair.end = j;
		}
	}

	for (j = (lgth2 - 2) >> 1; j >= 0; j--)
	{
		if ((g = lasthorizontalw[j]) > wm)
		{
			wm = g;
			iin = lgth1 - 1; jin = j;
			ijp[lgth1][lgth2] = -(lgth2 - j);
			pair.end = j;
		}
	}

	for (i = lgth1 - 2; i > (lgth1 - 2) >> 1; i--)
	{
		if ((g = lastverticalw[i]) >= wm)
		{
			wm = g;
			iin = i; jin = lgth2 - 1;
			ijp[lgth1][lgth2] = +(lgth1 - i);
		}
	}

	for (i = (lgth1 - 2) >> 1; i >= 0; i--)
	{
		if ((g = lastverticalw[i]) > wm)
		{
			wm = g;
			iin = i; jin = lgth2 - 1;
			ijp[lgth1][lgth2] = +(lgth1 - i);
			pair.end = jin;
		}
	}

	if (lasthorizontalw[lgth2 - 1] > wm)  // score ga onaji baai erabarenai
	{
		wm = lasthorizontalw[lgth2 - 1];
		iin = lgth1 - 1; jin = lgth2 - 1;
		ijp[lgth1][lgth2] = 0;
	}
	iin = lgth1; jin = lgth2;
	limk = lgth1 + lgth2 + 1;
	for (k = 0; k < limk; k++)
	{
		if (ijp[iin][jin] < 0)
		{
			ifi = iin - 1; jfi = jin + ijp[iin][jin];
		}
		else if (ijp[iin][jin] > 0)
		{
			ifi = iin - ijp[iin][jin]; jfi = jin - 1;
		}
		else
		{
			ifi = iin - 1; jfi = jin - 1;
		}

		l = iin - ifi;
		while (--l > 0)
		{
			//*--mseq1[0] = seq1[0][ifi + l];
			//*--mseq2[0] = *gap;
			k++;
		}
		l = jin - jfi;
		while (--l > 0)
		{
			//*--mseq1[0] = *gap;
			//*--mseq2[0] = seq2[0][jfi + l];
			k++;
		}
		if (iin <= 0 || jin <= 0) break;
		//*--mseq1[0] = seq1[0][ifi];
		//*--mseq2[0] = seq2[0][jfi];
		k++;
		iin = ifi; jin = jfi;
	}
	if (iin < 0 && jin < 0) ErrorExit("Error: can not use local alignment to determine the sequence");
	pair.start = jin;
	return pair;
}

local_align_pair SWAlign11(double** n_dynamicmtx, char** seq1, char** seq2, int alloclen)
// return the first place of sequence 2
// find pattern from sequence 1 to sequence 2
{

	//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lastj;
	int lgth1, lgth2;
	int resultlen, start_longer;
	double wm, wmo;   /* int ?????? */
	double g;
	double* currentw, * previousw;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	double fpenalty_tmp;

	double* wtmp;
	int* ijppt;
	double* mjpt, * prept, * curpt;
	int* mpjpt;

	static TLS double mi = 0.0;
	static TLS double* m = NULL;
	static TLS int** ijp = NULL;
	static TLS int mpi = 0;
	static TLS int* mp = NULL;
	static TLS double* w1 = NULL;
	static TLS double* w2 = NULL;
	static TLS double* match = NULL;
	static TLS double* initverticalw = NULL;    /* kufuu sureba iranai */
	static TLS double* lastverticalw = NULL;    /* kufuu sureba iranai */
	static TLS char** mseq1 = NULL;
	static TLS char** mseq2 = NULL;
	static TLS char** mseq = NULL;
	static TLS int** intwork = NULL;
	static TLS double** doublework = NULL;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double** amino_dynamicmtx = NULL; // ??

	double* wmrecords = NULL;
	double* prevwmrecords = NULL;
	double curm = 0.0;
	double* wmrecordspt, * wmrecords1pt, * prevwmrecordspt;

	local_align_pair ans = { 0, 0 };

 	lgth1 = strlen(seq1[0]);
	lgth2 = strlen(seq2[0]);

	if (lgth1 <= 0 || lgth2 <= 0)
	{
		fprintf(stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2);
	}

	if (lgth1 == 0 && lgth2 == 0)
		return ans;

	if (lgth1 == 0 || lgth2 == 0)
		return ans;


	wm = 0.0;

	if (orlgth1 == 0)
	{
		mseq1 = AllocateCharMtx(2, 0); 
		mseq2 = AllocateCharMtx(2, 0);
	}


	if (lgth1 > orlgth1 || lgth2 > orlgth2)
	{
		int ll1, ll2;

		ll1 = (int)(1.3 * lgth1) + 100;
		ll2 = (int)(1.3 * lgth2) + 100;

#if DEBUG
		fprintf(stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2);
#endif

		w1 = AllocateFloatVec(ll2 + 2);
		w2 = AllocateFloatVec(ll2 + 2);
		match = AllocateFloatVec(ll2 + 2);

		initverticalw = AllocateFloatVec(ll1 + 2);
		lastverticalw = AllocateFloatVec(ll1 + 2);

		m = AllocateFloatVec(ll2 + 2);
		mp = AllocateIntVec(ll2 + 2);

		mseq = AllocateCharMtx(2, ll1 + ll2); // 2020/Apr


		doublework = AllocateFloatMtx(nalphabets, MAX(ll1, ll2) + 2);
		intwork = AllocateIntMtx(nalphabets, MAX(ll1, ll2) + 2);
		amino_dynamicmtx = AllocateDoubleMtx(0x100, 0x100);

#if DEBUG
		fprintf(stderr, "succeeded\n");
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}

	for (i = 0; i < nalphabets; i++) for (j = 0; j < nalphabets; j++)
		amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];


	mseq1[0] = mseq[0];
	mseq2[0] = mseq[1];


#if DEBUG
		fprintf(stderr, "\n\ntrying to allocate %dx%d matrices ... ", orlgth1 + 1, orlgth2 + 1);
#endif

	commonIP = AllocateIntMtx(orlgth1 + 10, orlgth2 + 10);

#if DEBUG
		fprintf(stderr, "succeeded\n\n");
#endif


	ijp = commonIP;


#if 0
	for (i = 0; i < lgth1; i++)
		fprintf(stderr, "ogcp1[%d]=%f\n", i, ogcp1[i]);
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx(amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1);
	match_calc_mtx(amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2);


	for (j = 1; j < lgth2 + 1; ++j)
	{
		m[j] = currentw[j - 1]; mp[j] = 0;
	}

	if (lgth2 == 0)
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2 - 1]; // lgth2==0 no toki error

	lasti = lgth1;
	lastj = lgth2 + 1;

#if DEBUG
	fprintf(stderr, "currentw = \n");
	for (i = 0; i < lgth1 + 1; i++)
	{
		fprintf(stderr, "%5.2f ", currentw[i]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "initverticalw = \n");
	for (i = 0; i < lgth2 + 1; i++)
	{
		fprintf(stderr, "%5.2f ", initverticalw[i]);
	}
	fprintf(stderr, "\n");
#endif

	for (i = 1; i < lasti; i++)
	{
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i - 1];

		match_calc_mtx(amino_dynamicmtx, currentw, seq1, seq2, i, lgth2);
#if DEBUG
		fprintf(stderr, "\n");
		fprintf(stderr, "i=%d\n", i);
		fprintf(stderr, "currentw = \n");
		for (j = 0; j < lgth2; j++)
		{
			fprintf(stderr, "%5.2f ", currentw[j]);
		}
		fprintf(stderr, "\n");
#endif
#if DEBUG
		fprintf(stderr, "\n");
		fprintf(stderr, "i=%d\n", i);
		fprintf(stderr, "currentw = \n");
		for (j = 0; j < lgth2; j++)
		{
			fprintf(stderr, "%5.2f ", currentw[j]);
		}
		fprintf(stderr, "\n");
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;


		for (j = 1; j < lastj; j++)
		{

			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf(stderr, "%5.0f->", wm);
#endif
#if 0
			fprintf(stderr, "%5.0f?", g);
#endif
			if ((g = mi + fpenalty) > wm)
			{
				if (g >= 0)
				{
					wm = g;
					*ijppt = -(j - mpi);
				}
				else wm = 0;
			}
			if ((g = *prept) >= mi)
			{
				mi = g;
				mpi = j - 1;
			}

#if 0 
			fprintf(stderr, "%5.0f?", g);
#endif
			if ((g = *mjpt + fpenalty) > wm)
			{
				wm = g;
				if (g >= 0)
				{
					wm = g;
					*ijppt = +(i - *mpjpt);
				}
				else wm = 0;
			}
			if ((g = *prept) >= *mjpt)
			{
				*mjpt = g;
				*mpjpt = i - 1;
			}

			if (wm == 0) *ijppt = 0;

			* curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2 - 1]; // lgth2==0 no toki error
	}

	ans = Atracking__local_no_align(currentw, lastverticalw, lgth1, lgth2, ijp);

	free(mseq1);                     mseq1 = NULL;
	free(mseq2);                     mseq2 = NULL;
	FreeFloatVec(w1);                w1 = NULL;
	FreeFloatVec(w2);                w2 = NULL;
	FreeFloatVec(match);             match = NULL;
	FreeFloatVec(initverticalw);     initverticalw = NULL;
	FreeFloatVec(lastverticalw);     lastverticalw = NULL;
	FreeFloatVec(m);                 m = NULL;
	FreeIntVec(mp);                  mp = NULL;
	FreeCharMtx(mseq);               mseq = NULL;
	FreeFloatMtx(doublework);        doublework = NULL;
	FreeIntMtx(intwork);             intwork = NULL;
	FreeDoubleMtx(amino_dynamicmtx); amino_dynamicmtx = NULL;
	FreeIntMtx(commonIP);            commonIP = NULL;

	orlgth1 = orlgth2 = 0; 
	commonAlloc1 = commonAlloc2 = 0;

	return ans;
}
