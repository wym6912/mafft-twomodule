#include "mltaln.h"
#include "threadpool.h"
#include <time.h>

#define REPORTCOSTS 1

#define DEBUG 0
#define TESTCONST 0

#define ITERATIVECYCLE 2

#define END_OF_VEC -1

static int tuplesize, nunknown, exitval, aligncases;

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

#ifdef enablemultithread
typedef struct _dispatch_SWalign_local_arg
{
	char** seq0, ** seq2;
	pthread_mutex_t* lock, * lockout;
	int alloclen, id, maxseq1;
} dispatch_SWalign_local_arg;
#endif

void print_help()
{
	reporterr("Find similar fragments %d.%d.%d.%d%s Help:\n", VER_MAJOR, VER_MINOR, VER_RELEASE_SIM, VER_BUILD, VERSION);
	reporterr("-i: need to find fragment (destination) file with FASTA format\n");
	reporterr("-p: fragment file with FASTA format\n");
	reporterr("-T: use T threads to run this program\n");
	reporterr("-f, -g, -h: ppenalty, ppenalty_ex(not used), poffset(not used)\n");
	reporterr("-Q, -V: penalty_shift_factor(not used), ppenalty_dist\n");
	reporterr("-b: BLOSUM Matrix\n");
	reporterr("-j: use jtt/kimura model, pamN is needed\n");
	reporterr("-m: use tm model, pamN is needed\n");
	reporterr("-d: Print Score Matrix fits the destination file and exit\n");
	reporterr("-D, -P: -D the sequence is DNA, -P the sequence is Protein\n");
	reporterr("-S: Calcuate SP Scores after alignment\n");
	reporterr("-v: show program version and exit\n");
	reporterr("-H, -?: Print help message and exit\n");
}

void print_version()
{
	reporterr("findsim %d.%d.%d.%d%s\n", VER_MAJOR, VER_MINOR, VER_RELEASE_SIM, VER_BUILD, VERSION);
}

void arguments( int argc, char *argv[] )
{
	int c;
	nthread = 1;
	outnumber = 0;
	nevermemsave = 0;
	inputfile = NULL;
	targetfile = NULL;
	nblosum = 62;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = -1530;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	TMorJTT = JTT;
	disp = 0;
	
	while( --argc > 0 && (*++argv)[0] == '-' )
	{
		while ( (c = *++argv[0]) )
		{
			switch( c )
			{
				case 'i': // common sequence file per line, must not have space(s)
					inputfile = *++argv;
					reporterr( "destination sequence file = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'p':
					targetfile = *++ argv;
					reporterr( "target sequence file = %s\n", targetfile );
					-- argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( myatof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = myatof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( myatof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( myatof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					reporterr(       "blosum %d / kimura 200 \n", nblosum );
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
				case 'T':
					nthread = myatoi( *++argv );
					reporterr(       "nthread = %d\n", nthread );
					-- argc; 
					goto nextoption;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'H':
				case '?':
					print_help();
					exit(0);
				case 'v':
					print_version();
					exit(0);
				case 'd':
					disp = 1;
					break;
				default:
					reporterr(       "illegal option %c\n", c );
					argc = 0;
					break;
			}
		}
		nextoption:
			;
	}
	if( argc != 0 ) 
	{
		reporterr( "options: Check source file !\n" );
		exit( 1 );
	}
}

#ifdef enablemultithread
int finished;
void *dispatch_SWalign11_local(void *arg)
{
#define F(X) (((dispatch_SWalign_local_arg *)arg) -> X)
	char **seq0 = F(seq0), **seq2 = F(seq2);
	int alloclen = F(alloclen), id = F(id), maxseq1 = F(maxseq1);
	pthread_mutex_t* lock = F(lock), *lockout = F(lockout);
#undef F
	static local_align_pair lp;
	static int i, whichseq;
	whichseq = -1;
	for(i = 0; i < maxseq1; ++ i)
		if (!pthread_mutex_trylock(lock + i))
		{
			whichseq = i;
			break;
		}
	if (whichseq == -1) ErrorExit("Error: cannot make Smith-Waterman alignment!");

	lp = SWAlign11(n_dis_consweight_multi, seq0 + whichseq, seq2, strlen(seq2[0]));
	pthread_mutex_unlock(lock + whichseq);

	pthread_mutex_lock(lockout);
	fprintf(stdout, "%d %d %d\n", id, lp.start, lp.end);

	++ finished;
	if(finished % 10 == 0 || finished == njob) reporterr("\r   %d / %d", finished, njob);
	pthread_mutex_unlock(lockout);
	
	return NULL;
}
#endif

int main(int argc, char **argv)
{
	int nlen1, *nlen2 = NULL;
	char **seq1 = NULL, **name1 = NULL, **seq2 = NULL, **name2 = NULL, **seq0;
	int alloclen, i, j, k;
	FILE *infp = NULL, *tafp = NULL;
	arguments(argc, argv);
	if(inputfile && targetfile)
	{
		infp = fopen(inputfile, "rb");
		if(! infp)
		{
			reporterr("Error: cannot open destination sequence file %s\n", inputfile);
			exit(1);
		}
		tafp = fopen(targetfile, "rb");
		if(! tafp)
		{
			reporterr("Error: cannot open target sequence file %s\n", targetfile);
			exit(1);
		}
	}
	else
	{
		ErrorExit("In find similar program, you must specify the destination sequence file and the target sequence file. ");
	}
#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize( (unsigned long long)1024 * 1024 * 1024 ); // topolorder() de ookime no stack wo shiyou.
#endif

	/* Part 1: read target file and destination file */
	getnumlen(tafp);
	if (njob != 1)
	{
		fclose(tafp);
		fclose(infp);
		reporterr("Error: the target file %s must have only 1 sequence. Please check this file and try again.\n", targetfile);
		exit(1);
	}
	rewind(tafp);
	
	name1 = AllocateCharMtx(1, BLEN + 5);
	seq1 = AllocateCharMtx(1, nlenmax + 5);
	readData_pointer(tafp, name1, &nlen1, seq1);
	fclose(tafp);

	getnumlen(infp);
	if (njob == 0)
	{
		fclose(infp);
		reporterr("Warning: no sequences found in %s. Program will exit\n", inputfile);
		exit(0);
	}
	rewind(infp);

	name2 = AllocateCharMtx(njob, BLEN + 5);
	seq2 = AllocateCharMtx(njob, nlenmax + 5);
	nlen2 = AllocateIntVec(njob);
	readData_pointer(infp, name2, nlen2, seq2);
	fclose(infp);
	constants(njob, seq2);
#if REPORTCOSTS
		time_t starttime, startclock;
		starttime = time(NULL);
		startclock = clock();
#endif

	/* Part 2: Use Smith-Waterman Algorithm to find the most similar fragment */
	reporterr("Finding head for destination...\n");

#ifdef enablemultithread
	seq0 = AllocateCharMtx(MIN(njob, nthread), (int)(strlen(seq1[0])) + 1);
	pthread_mutex_t *lock, lockout;
	lock = malloc(sizeof(pthread_mutex_t) * (MIN(njob, nthread) + 1));
	if(lock == NULL) ErrorExit("Error: can't allocate enough space for multi-thread mode\n");
	for (i = 0; i < MIN(njob, nthread); ++i)
	{
		pthread_mutex_init(lock + i, NULL);
		strcpy(seq0[i], seq1[0]);
	}
	pthread_mutex_init(&lockout, NULL);
	threadpool_t tp;
	threadpool_init(&tp, nthread);
	dispatch_SWalign_local_arg *_SWalign_arg_;
	_SWalign_arg_ = malloc(sizeof(dispatch_SWalign_local_arg) * njob);
	if (_SWalign_arg_ == NULL) ErrorExit("Error: can't allocate enough space for multi-thread mode\n");
	finished = 0;
	for(i = 0; i < njob; ++ i)
	{
		_SWalign_arg_[i].alloclen = (int)(strlen(seq1[0]));
		_SWalign_arg_[i].id = i;
		_SWalign_arg_[i].lock = lock;
		_SWalign_arg_[i].lockout = &lockout;
		_SWalign_arg_[i].seq0 = seq0;
		_SWalign_arg_[i].seq2 = &seq2[i];
		_SWalign_arg_[i].maxseq1 = MIN(njob, nthread);
		threadpool_add_task(&tp, dispatch_SWalign11_local, _SWalign_arg_ + i);
	}
	threadpool_destroy(&tp);
	for (i = 0; i < MIN(njob, nthread); ++ i) pthread_mutex_destroy(lock + i);
	pthread_mutex_destroy(&lockout);
	free(_SWalign_arg_);
	free(lock);
	FreeCharMtx(seq0);
	reporterr("\ndone. \n");

#else
	local_align_pair lp;
	for (i = 0; i < njob; ++i)
	{
		lp = SWAlign11(n_dis_consweight_multi, seq1, &seq2[i], strlen(seq2[i]));
		fprintf(stdout, "%d %d %d\n", i, lp.start, lp.end);
		reporterr("\r   %d / %d", i + 1, njob);
	}
	reporterr("\ndone. \n");
#endif

	/* Part end: free and show info */
#if REPORTCOSTS
//		use_getrusage();
		reporterr( "\nfindsim, real = %f min\n", (float)(time(NULL) - starttime)/60.0 );
		reporterr( "findsim, user = %f min\n", (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
	reporterr( "%s (%s, %d-bit) Version " , progName( argv[0] ), (dorp=='d')?"nuc":((nblosum==-2)?"text":"aa"), sizeof(int *) * 8 );
	reporterr( "%d.%d.%d.%d", VER_MAJOR, VER_MINOR, VER_RELEASE_SIM, VER_BUILD );
	reporterr( "%s \nmodel=%s, amax=%3.1f\n%d thread(s)\n\n", VERSION, modelname, specificityconsideration, nthread );
#ifndef enablemultithread
	reporterr("...but NOT used multi threads in find similar.\n\n\n");
#endif
	freeconstants();
	FreeCharMtx(seq1);
	FreeCharMtx(seq2);
	FreeCharMtx(name1);
	FreeCharMtx(name2);
	FreeIntVec(nlen2);
	return 0;
}
