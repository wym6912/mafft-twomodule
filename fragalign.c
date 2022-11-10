#include "mltaln.h"
#include "threadpool.h"
#include <time.h>


#define REPORTCOSTS 1

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SKIP 1

#define ITERATIVECYCLE 2

#define END_OF_VEC -1


#if 0
#define PLENFACA 0.0123
#define PLENFACB 10252
#define PLENFACC 10822
#define PLENFACD 0.5
#define DLENFACA 0.01
#define DLENFACB 2445
#define DLENFACC 2412
#define DLENFACD 0.1
#else
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0
#endif

#ifdef enablemultithread
typedef struct _Falign_arg
{
	char **mseq1, **mseq2;
	int clus1, clus2, alloclen, fftlog;
	double *eff;
} Falign_arg;

typedef struct _merge_arg
{
	int id, depth, njob, alloclen;
	char **centerseq, **bseq, **common;
} merge_arg;
#endif

int nmax_shift_factor, longest_id;
void print_help_message()
{
	reporterr("Fragment Align %d.%d.%d.%d%s\n", VER_MAJOR, VER_MINOR, VER_RELEASE_FRAG, VER_BUILD, VERSION);
	reporterr("This program is based on MAFFT, only give one file to the program:\n");
	reporterr("- The fragment file, which the file has more than two sequences.\n");
	reporterr("We align this file using FFT local alignment.\n");
	reporterr("Arguments: \n");
	reporterr("-i: file with FASTA format\n");
	reporterr("-f, -g, -h: ppenalty, ppenalty_ex(not used), poffset(not used)\n");
	reporterr("-Q, -V: penalty_shift_factor (not used), ppenalty_dist\n");
	reporterr("-b: BLOSUM Matrix\n");
	reporterr("-j: use jtt/kimura model, pamN is needed\n");
	reporterr("-m: use tm model, pamN is needed\n");
	reporterr("-a: 0 is default model, 1 is fmodel, -1 is raw model\n");
	reporterr("-D, -P: -D the sequence is DNA, -P the sequence is Protein\n");
	reporterr("-S: Calcuate SP Scores after alignment\n");
	reporterr("-d: Print Score Matrix fits the common file and exit\n");
	reporterr("-z, -w: FFT align arguments: fftthreshold, fftWinSize\n");
	reporterr("-B: Kband in calcuating DP-matrix during the alignment. If not defined this varible, use auto double\n");
	reporterr("-T: Use T threads to run this program\n");
	reporterr("-s: nmax shift factor, it means the times of the max length\n");
	reporterr("-L: Use legacy gap cost in order to get less gaps\n");
	reporterr("-v: show program version\n");
	reporterr("-H, -?: Print this help message and exit\n");
}

void print_version()
{
	reporterr("fragalign %d.%d.%d.%d%s\n", VER_MAJOR, VER_MINOR, VER_RELEASE_FRAG, VER_BUILD, VERSION);
}

void arguments( int argc, char *argv[] )
{
	int c;
	nthread = 1;
	inputfile = NULL;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	devide = 0;
	use_fft = 1;
	fftscore = 1;
	utree = 1;
	tbutree = 1;
	refine = 0;
	check = 1;
	cut = 0.0;
	disp = 0;
	outgap = 1;
	alg = 'A';
	mix = 0;
	kobetsubunkatsu = 0;
	sueff_global = 0.1;
	scoremtx = 1;
	dorp = NOTSPECIFIED;
	ppenalty_dist = NOTSPECIFIED;
	ppenalty = -1530;
	ppenalty_ex = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	spscoreout = 0;
	nmax_shift_factor = 1;
	legacygapcost = 0;
	alignband = NOTSPECIFIED;

	while( --argc > 0 && (*++argv)[0] == '-' )
	{
		while ( (c = *++argv[0]) )
		{
			switch( c )
			{
				case 'i':
					inputfile = *++argv;
					reporterr("inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'V':
					ppenalty_dist = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					reporterr(       "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					reporterr(       "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					reporterr(       "tm %d\n", pamN );
					--argc;
					goto nextoption;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
				case 'd':
					disp = 1;
					break;
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'B':
					alignband = myatoi(* ++ argv);
					-- argc;
					goto nextoption;
				case 'T':
					nthread = myatoi(*++argv);
					-- argc;
					goto nextoption;
				case 's':
					nmax_shift_factor = myatoi(*++ argv);
					-- argc;
					goto nextoption;
				case 'L':
					legacygapcost = 1;
					break;
				case 'H':
				case '?':
					print_help_message();
					exit(0);
				case 'v':
					print_version();
					exit(0);
				default:
					reporterr(       "illegal option %c\n", c );
					argc = 0;
					break;
			}
		}
		nextoption:
			;
	}
	if( argc == 1 )
	{
		cut = atof( (*argv) );
		argc--;
	}
	if( argc != 0 ) 
	{
		reporterr(       "options: Check source file !\n" );
		exit( 1 );
	}
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
	reporterr(       "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
	if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "fragalign method\n" );
   	fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );


	if( use_fft )
	{
		fprintf( fp, "FFT on\n" );
		if( dorp == 'd' )
			fprintf( fp, "Basis : 4 nucleotides\n" );
		else
		{
			if( fftscore )
				fprintf( fp, "Basis : Polarity and Volume\n" );
			else
				fprintf( fp, "Basis : 20 amino acids\n" );
		}
		fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
		fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
	}
	else
		fprintf( fp, "FFT off\n" );
	fflush( fp );
}

/* list start from 0 */
void insertgaplist(char *res, char *src, int *list)
{
	int i = 0, l;
	char ch; char gapchar = *newgapstr;
	while(ch = *src)
	{
		if((l = *list) == i)
		{
			*res ++ = gapchar;
			++ list;
		}
		else
		{
			*res ++ = *src ++;
			++ i;
		}
	}
	while((l = *list) != -1)
	{
		*res ++ = gapchar;
		++ list;
	}
	*res = 0;
}

#ifdef enablemultithread
int cnt_OK;
void* dispatch_Falign(void* arg)
{
#define F(X) (((Falign_arg *)arg) -> X)
	Falign(NULL, NULL, n_dis_consweight_multi, F(mseq1), F(mseq2), F(eff), F(eff), NULL, NULL, 1, 1, F(alloclen), &F(fftlog), NULL, 0, NULL);
#undef F

	++ cnt_OK;
	if (cnt_OK % 100 == 0 || cnt_OK >= njob - 1 || cnt_OK == 1) reporterr("\rSTEP %d / %d ", cnt_OK, njob - 1);
	return NULL;
}
#endif


char** mergeallresult(char** resultseq, char** common, char** center, int njob, int alloclen)
{
	// insert the first two sequences
	strcpy(resultseq[0], center[0]); strcpy(resultseq[1], common[0]);
	int i, j, k, len1, len2, p1, p2, l, edit1, restore;
	char* tmp, ** swap;
	int* list1, * list2;
	tmp = AllocateCharVec(alloclen);
	list1 = AllocateIntVec(alloclen);
	list2 = AllocateIntVec(alloclen);
	for (i = 1; i < njob; ++i) // resultseq has i sequences
	{
		// merge common[i] && center[i] to resultseq
		len1 = strlen(resultseq[0]); len2 = strlen(center[i]);
		j = 0, k = 0;
		p1 = -1, p2 = -1;
		edit1 = 0;
		restore = 0;
		while (j < len1 && k < len2)
		{
			if (resultseq[0][j] == center[i][k])
			{
				++j;
				++k;
			}
			else if (resultseq[0][j] == *newgapstr)
			{
				++j;
				list2[++p2] = k;
			}
			else if (center[i][k] == *newgapstr)
			{
				list1[++p1] = j;
				++k;
				++edit1;
			}
			else
			{
				++j;
				++k;
			}
		}
		// add tail gaps
		while (j < len1)
		{
			++j;
			list2[++p2] = len2;
		}
		while (k < len2)
		{
			++k;
			list1[++p1] = len1;
			++edit1;
		}

		list1[++p1] = -1; list2[++p2] = -1;
#if DEBUG
		reporterr("%d %d %d %d\n", j, len1, k, len2);
		reporterr("seqs = \n%s\ncenter = %s\ncommon = %s\n", resultseq[0], center[i], common[i]);
		reporterr("list1 = "); for (l = 0; l <= p1; ++l, reporterr(" ")) reporterr("%d", list1[l]); reporterr("\n");
		reporterr("list2 = "); for (l = 0; l <= p2; ++l, reporterr(" ")) reporterr("%d", list2[l]); reporterr("\n");
#endif
		if (len1 + p1 >= alloclen || len2 + p2 >= alloclen)
		{
			reporterr("Reallocing... ");
			swap = AllocateCharMtx(njob + 1, alloclen + 1);
			for (l = 0; l <= i; ++l) strcpy(swap[l], resultseq[l]);
			alloclen = MAX(len1 + p1 + 10, len2 + p2 + 10) << 1;
			ReallocateCharMtx(resultseq, njob + 2, alloclen);
			for (l = 0; l <= i; ++l) strcpy(resultseq[l], swap[l]);
			free(tmp);
			tmp = AllocateCharVec(alloclen);
			FreeCharMtx(swap);
			reporterr("done.\n");
			restore = 1;
		}
		if (edit1)
		{
			for (l = 0; l <= i; ++l)
			{
				strcpy(tmp, resultseq[l]);
				insertgaplist(resultseq[l], tmp, list1);
			}
		}
		insertgaplist(resultseq[i + 1], common[i], list2);
#if DEBUG
		reporterr("After align, resultseq = \n");
		for (l = 0; l <= i + 1; ++l) reporterr("%s\n", resultseq[l]);
#endif
		if (restore)
		{
			FreeIntVec(list1);
			FreeIntVec(list2);
			list1 = AllocateIntVec(alloclen);
			list2 = AllocateIntVec(alloclen);
		}
		if (i == 1 || i == njob - 1 || i % 100 == 0) reporterr("\rSTEP %d / %d ", i, njob - 1);
	}
	FreeIntVec(list1);
	FreeIntVec(list2);
	free(tmp);
	return resultseq;
}


int main(int argc, char** argv)
{
	int* nlen = NULL; // the length of common file and center file	
	char** name = NULL, ** seq = NULL, ** centerseq = NULL;
	char** bseq = NULL;
	int i, j;
	FILE* infp = NULL;
	char c;
	int alloclen, maxlen;
	int* grpseq = NULL;
	char* tmpseq = NULL;
	double* eff = NULL;
	int fftlog;
	char b[B];

	arguments(argc, argv);

	if (inputfile)
	{
		infp = fopen(inputfile, "rb");
		if (!infp)
		{
			reporterr("Cannot open sequence file %s\n", inputfile);
			exit(1);
		}
	}
	else
	{
		reporterr("In fragalign, you must specify sequence file!\n");
		exit(1);
	}


	getnumlen(infp);
	rewind(infp);

	/* fragment aligner support for more than 2 fragments alignment */
	if (njob == 0)
	{
		reporterr("Error: Sequence file %s has no sequences. Please check this file.\n", inputfile);
		exit(1);
	}
	else if (njob == 1)
	{
		reporterr("Info: Sequence file %s has only one sequence. There is no need to use Fragment alignment.\n", inputfile);
		Filecopy(infp, stdout);
		exit(0);
	}

#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize((unsigned long long)1024 * 1024 * 1024); // topolorder() de ookime no stack wo shiyou.
#endif
	if ((long)(nlenmax << nmax_shift_factor) > (int)(1e7))
	{
		reporterr("nmax_shift_factor is to long. please make it smaller. \n");
		exit(1);
	}
	nlenmax <<= nmax_shift_factor;
	seq = AllocateCharMtx(njob, nlenmax + 1);
#if SAFE
	mseq1 = AllocateCharMtx(1, 0);
	mseq2 = AllocateCharMtx(1, 0);
#endif
	eff = AllocateDoubleVec(1);
	bseq = AllocateCharMtx(njob + 2, nlenmax + 10);
	centerseq = AllocateCharMtx(njob + 2, nlenmax + 10);

	name = AllocateCharMtx(njob, B + 1);
	nlen = AllocateIntVec(njob);
	readData_pointer(infp, name, nlen, seq);
	fclose(infp);
	constants(njob, seq);

	
#if 0 //maybe useless now (2021/09/14)
	initSignalSM();

	initFiles();

	WriteOptions( trap_g );
#endif
	c = seqcheck( seq );
	if( c )
	{
		reporterr("Illegal character %c on sequence file.\n", c);
		exit( 1 );
	}
	reporterr("\n");
	alloclen = nlenmax;
	*eff = 1.0;
	reporterr("Fragment alignment, the file has %d sqeuences\n", njob );

	//find the longest fragment sequence
	longest_id = -1;
	for (i = 0; i < njob; ++i)
		if (strlen(seq[i]) == (nlenmax >> nmax_shift_factor))
		{
			longest_id = i;
			reporterr("found longest sequence %d\n", longest_id);
			break;
		}
	if (longest_id == -1) { reporterr("Not found longest sequence???"); exit(1); }

	// copy the longest fragment sequence to center sequence, and make center star alignment!
	for (i = 0; i < njob; ++i) strcpy(centerseq[i], seq[longest_id]);

#if REPORTCOSTS
	time_t starttime, startclock;
	starttime = time(NULL);
	startclock = clock();
#endif
	reporterr("Calcuating the pairwise alignment... \n");
#ifdef enablemultithread
	// Star alignment in multi thread
	cnt_OK = 0;
	Falign_arg *_Falign_arg_;
	_Falign_arg_ = malloc(njob * sizeof(Falign_arg));
	threadpool_t tp;
	threadpool_init(&tp, nthread);
	for(i = 0; i < njob; ++ i)
	{
		if (i == longest_id) continue;
		_Falign_arg_[i].mseq1 = &seq[i];
		_Falign_arg_[i].mseq2 = &centerseq[i];
		_Falign_arg_[i].eff = eff;
		_Falign_arg_[i].alloclen = alloclen;
		threadpool_add_task(&tp, dispatch_Falign, _Falign_arg_ + i);
	}
	threadpool_destroy(&tp);
	free(_Falign_arg_);
#else
	// Star alignment in one thread
	for(i = 0; i < njob; ++ i)
	{
		if (i == longest_id) continue;
		// reporterr("seq1 = %s, centerseq = %s\n", seq[i], centerseq[i]);
#if SAFE
		strcpy(mseq1[0], seq[i]);  strcpy(mseq2[0], centerseq[i]);
		Falign(NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, eff, eff, NULL, NULL, 1, 1, alloclen, &fftlog, NULL, 0, NULL);
		strcpy(seq[i], mseq1[0]);  strcpy(centerseq[i], mseq2[0]);
#else
		Falign(NULL, NULL, n_dis_consweight_multi, &seq[i], &centerseq[i], eff, eff, NULL, NULL, 1, 1, alloclen, &fftlog, NULL, 0, NULL);
#endif
		// reporterr("After align, seq1 = %s, center = %s\n", seq[i], centerseq[i]);
	}
#endif
	reporterr("done. OK on aligning.\n");
#if 0 // REPORTCOSTS
//		use_getrusage();
		reporterr( "\nAfter Falign, real = %f min\n", (float)(time(NULL) - starttime)/60.0 );
		reporterr( "After Falign, user = %f min\n", (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
	// merge the result
	reporterr("Start merging...\n");
#if 1
#if DEBUG
	reporterr("centerseq = "); for(i = 0; i < njob; ++ i) reporterr("%s ", centerseq[i]);
	reporterr("\ncommonseq = "); for(i = 0; i < njob; ++ i) reporterr("%s ", seq[i]);
	reporterr("\n");
#endif
	bseq = mergeallresult(bseq, seq, centerseq, njob, alloclen);
#else
#ifdef ALGORITHM1
	cnt_OK = 0;
	static TLS merge_arg *_merge_arg_;
	_merge_arg_ = malloc(njob * sizeof(merge_arg));
	int lognj = ceil(log2(njob));
	for(i = 1; i <= lognj; ++ i)
	{
		threadpool_init(&tp, nthread);
		// reporterr("i = %d, j = %d, njob = %d\n", i, njob / lognj, njob);
		for(j = 0; j < njob; j += 1 << i)
		{
			_merge_arg_[j].id = j;
			_merge_arg_[j].centerseq = centerseq;
			_merge_arg_[j].depth = i;
			_merge_arg_[j].njob = njob;
			_merge_arg_[j].common = seq;
			_merge_arg_[j].bseq = bseq;
			_merge_arg_[j].alloclen = alloclen;
			// reporterr("ID = %d, depth = %d, segment = [%d, %d]\n", j, i, j, MIN(j + (1 << i) - 1, njob - 1));
			threadpool_add_task(&tp, dispatch_mergecenter, (void *)(_merge_arg_ + j));
			// dispatch_mergecenter((void *)(_merge_arg_ + j));
		}
		threadpool_destroy(&tp);
		reporterr("\rSTEP %d / %d ", i, lognj);
	}

	free(_merge_arg_);
	
	strcpy(bseq[0], centerseq[0]);
	// for(i = 0; i < njob; ++ i) gapireru(bseq[i + 1], seq[i], centerseq[i]);
#else	
	bseq = mergeallresult(bseq, seq, centerseq, njob, alloclen);
#endif
#endif
	reporterr("done. \n");
	strcpy(bseq[longest_id + 1], centerseq[0]);
	// in mergeallresult(), the first sequence of bseq is the center sequence, ignore it when processing the following steps
#if REPORTCOSTS
//	use_getrusage();
	reporterr( "\nfragalign, real = %f min\n", (float)(time(NULL) - starttime)/60.0 );
	reporterr( "fragalign, user = %f min\n", (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
	reporterr(       "\ndone.\n\n" );

	int len0 = strlen(bseq[0]), len1;
	for (i = 1; i <= njob; ++i)
	{
		len1 = strlen(bseq[i]);
		if (len1 == 0)
		{
			for (j = len1; j < len0; ++j) bseq[i][j] = *newgapstr;
			bseq[i][len0] = 0;
		}
		else if (len1 > len0)
		{
			reporterr("ERROR: NOT ALIGNED. Please concat the author and submit your sequences.\n");
			exit(1);
		}
	}
	writeData_pointer( stdout, njob, name, nlen, &bseq[1] );

	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob + 1, bseq ) );
	reporterr( "%s (%s, %d-bit) Version" , progName( argv[0] ), (dorp=='d')?"nuc":((nblosum==-2)?"text":"aa"), sizeof(int *) * 8 );
	reporterr( "%d.%d.%d.%d", VER_MAJOR, VER_MINOR, VER_RELEASE_FRAG, VER_BUILD );
	reporterr( "%s \nalg=%c, model=%s, amax=%3.1f\n%d thread(s)\n\n", VERSION, alg, modelname, specificityconsideration, nthread );

	// FreeCharMtx( bseq );
	FreeCharMtx(centerseq);
	FreeCharMtx( name );
	FreeIntVec( nlen );
	FreeDoubleVec(eff);	
	FreeCharMtx( seq );
#if SAFE
	free( mseq1 );
	free( mseq2 );
#endif
	freeconstants();
#if 0
	closeFiles();
#endif
	FreeCommonIP();

	return 0;
}
