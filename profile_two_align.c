#include "mltaln.h"
#include <time.h>

#define REPORTCOSTS 1

#define DEBUG 0

#define END_OF_VEC -1

static int nunknown, exitval, aligncases, print_to_two_files;

char *profilename1, *profilename2;

void print_help()
{
	reporterr("2-profiles alignment %d.%d.%d.%d%s Help:\n", VER_MAJOR, VER_MINOR, VER_RELEASE_PROF, VER_BUILD, VERSION);
	reporterr("make alignment between profile 1 and profile 2, and write result to profile 1 without -M argument.\n");
	reporterr("-p: profile 1 file name\n");
	reporterr("-q: profile 2 file name\n");
	reporterr("-f, -g, -h: ppenalty, ppenalty_ex(not used), poffset(not used)\n");
	reporterr("-Q, -V: penalty_shift_factor(not used), ppenalty_dist\n");
	reporterr("-b: BLOSUM Matrix\n");
	reporterr("-j: use jtt/kimura model, pamN is needed\n");
	reporterr("-m: use tm model, pamN is needed\n");
	reporterr("-D, -P: -D the sequence is DNA, -P the sequence is Protein\n");
	reporterr("-S: Calcuate SP Scores after alignment\n");
	reporterr("-z, -w: FFT align arguments: fftthreshold, fftWinSize\n");
	reporterr("-B: Kband in calcuating DP-matrix during the alignment\n");
	reporterr("-A: Use Aalign to align sequences\n");
	reporterr("-F: Use FFT align to align sequences\n");
	reporterr("-v: show program version and exit\n");
	reporterr("-M: print profile result to two profiles\n");
	reporterr("-H, -?: Print help message and exit\n");
}

void print_version()
{
	reporterr("2-profiles align %d.%d.%d.%d%s\n", VER_MAJOR, VER_MINOR, VER_RELEASE_PROF, VER_BUILD, VERSION);
}

void arguments( int argc, char *argv[] )
{
	int c;
	alg = 'A';
	nthread = 1;
	outnumber = 0;
	nevermemsave = 0;
	inputfile = NULL;
	nblosum = 62;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = -1530;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	spscoreout = 0;
	aligncases = 1;
	alignband = NOTSPECIFIED;
	print_to_two_files = 0;

	while( --argc > 0 && (*++argv)[0] == '-' )
	{
		while ( (c = *++argv[0]) )
		{
			switch( c )
			{
				case 'p':
					profilename1 = *++ argv;
					reporterr( "profile 1 file = %s\n", profilename1 );
					--argc;
					goto nextoption;
				case 'q':
					profilename2 = *++ argv;
					reporterr( "profile 2 file = %s\n", profilename2 );
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
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kimura model, kappa = %d\n", kimuraR );
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
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
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
					alignband = myatoi( *++argv );
					-- argc;
					goto nextoption;
				case 'A':
					aligncases = 1; // A__align11
					reporterr("Use Aalign\n");
					break;
				case 'F':
					aligncases = 0; // Falign
					reporterr("Use FFT Align\n");
					break;
				case 'M':
					print_to_two_files = 1;
					break;
				case 'H':
				case '?':
					print_help();
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
	if( argc != 0 )
	{
		reporterr( "options: Check source file !\n" );
		exit( 1 );
	}
}

int main(int argc, char **argv)
{
	int *nlen = NULL;
	char **seq = NULL, **name = NULL, **seq2 = NULL, **name2 = NULL;
	int maxlen, alloclen, i, j, k, fftlog, centerseqs;
	int f1seq, f2seq; // sequence numbers
	int f1len, f2len; // max length of two sequences
	int *grpseq = NULL, **pointt = NULL, *table1 = NULL, *nogaplen = NULL, *nlen22;
	char *tmpseq = NULL, *align1 = NULL, *align2 = NULL;
	double *eff = NULL, **mtx = NULL, **nlen2 = NULL, *eff2 = NULL;
	char *sgap1 = NULL, *sgap2 = NULL, *egap1 = NULL, *egap2 = NULL;
	char b[BLEN];
	FILE *prof1 = NULL, *prof2 = NULL;
	arguments(argc, argv);

#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize( (unsigned long long)1024 * 1024 * 1024 ); // topolorder() de ookime no stack wo shiyou.
#endif

	if(profilename1 == NULL || profilename2 == NULL)
	{
		reporterr("Error: Missing one profile name. Please check your arguments -p and -q.\n");
		exit(1);
	}

	/* use changed FFT (changed to calcuate density instead of sequence) to align profile data */
	reporterr("Aligning two profiles %s and %s...\n", profilename1, profilename2);
	// There's no need to use multithread mode. Single thread on profile alignment
	prof1 = fopen(profilename1, "rb");
	prof2 = fopen(profilename2, "rb");
	if(prof1 == NULL || prof2 == NULL)
	{
		reporterr("Error: can not open profile %s or %s.\n", profilename1, profilename2);
		exit(1);
	}
	getnumlen(prof1);
	rewind(prof1);
	f1seq = njob;
	maxlen = nlenmax;
	getnumlen(prof2);
	rewind(prof2);
	f2seq = njob;
	if (f1seq == 0 || f2seq == 0)
	{
		reporterr("Warning: one profile has only 0 sequences. There is no need to make 2-profiles align.\n");
		if(! print_to_two_files && f1seq == 0)
		{
			// move f2 to f1
			reporterr("Move profile %s to %s...\n", profilename2, profilename1);
			fclose(prof1); fclose(prof2);
			prof1 = fopen(profilename1, "wb");
			prof2 = fopen(profilename2, "rb");
			if(prof1 == NULL || prof2 == NULL)
			{
				reporterr("Error: can not move profile %s to %s.\n", profilename2, profilename1);
				exit(1);
			}
			Filecopy(prof2, prof1);
			fclose(prof1); fclose(prof2);
			prof2 = fopen(profilename2, "wb");
			if(prof2 == NULL)
			{
				reporterr("Warning: can not clean profile %s. Please clean in manually.\n", profilename2);
				exit(0);
			}
			// do nothing for cleaning prof2
			fclose(prof2);
		}
		exit(0);
	}
	maxlen = MAX(maxlen, nlenmax);
	nlenmax = maxlen;
	seq = AllocateCharMtx(f1seq, nlenmax << 1);
	name = AllocateCharMtx(f1seq, BLEN);
	nlen = AllocateIntVec(f1seq);
	seq2 = AllocateCharMtx(f2seq, nlenmax << 1);
	name2 = AllocateCharMtx(f2seq, BLEN);
	nlen22 = AllocateIntVec(f2seq);
	readData_pointer2(prof1, f1seq, name, nlen, seq);
	readData_pointer2(prof2, f2seq, name2, nlen22, seq2);
	fclose(prof1);
	fclose(prof2);
	constants(njob, seq);
	// Calcuate the effect value: average of all the sequence
	eff = AllocateDoubleVec(f1seq);
	eff2 = AllocateDoubleVec(f2seq);
	if(f1seq > 0) for(j = 0; j < f1seq; ++ j) eff[j] = 1.0 / f1seq;
	else { reporterr("Error: the sequence file %s has no sequences! It may coursed by the smaller value of fftWinsize, please make it larger. The arugment of fftWinsize is -w.\n", profilename1); exit(1); }
	if(f2seq > 0) for(j = 0; j < f2seq; ++ j) eff2[j] = 1.0 / f2seq;
	else { reporterr("Error: the sequence file %s has no sequences! It may coursed by the smaller value of fftWinsize, please make it larger. The arugment of fftWinsize is -w.\n", profilename2); exit(1); }
	f1len = nlen[0];
	for(j = 1; j < f1seq; ++ j) f1len = MAX(nlen[j], f1len);
	f2len = nlen22[0];
	for(j = 1; j < f2seq; ++ j) f2len = MAX(nlen22[j], f2len);

#if REPORTCOSTS
	time_t starttime, startclock;
	starttime = time(NULL);
	startclock = clock();
#endif

	alloclen = f1len + f2len + 10;
	if(aligncases == 1)
	{
		sgap1 = AllocateCharVec(f1seq + 10);
		sgap2 = AllocateCharVec(f2seq + 10);
		egap1 = AllocateCharVec(f1seq + 10);
		egap2 = AllocateCharVec(f2seq + 10);
		memset(sgap1, 'o', f1seq * sizeof(char));
		memset(sgap2, 'o', f2seq * sizeof(char));
		memset(egap1, 'o', f1seq * sizeof(char));
		memset(egap2, 'o', f2seq * sizeof(char));
		A__align(n_dis_consweight_multi, penalty, penalty_ex, seq, seq2, eff, eff2, f1seq, f2seq, alloclen, sgap1, sgap2, egap1, egap2, 1, 1);
		free(sgap1);
		free(sgap2);
		free(egap1);
		free(egap2);
	}
	else if(aligncases == 0)
		Falign(NULL, NULL, n_dis_consweight_multi, seq, seq2, eff, eff2, NULL, NULL, f1seq, f2seq, alloclen, &fftlog, NULL, 0, NULL);
	else ErrorExit("ERROR: aligncases is error. Please check your command.\n");
	if(! print_to_two_files) reporterr("Writing alignment to %s...\n", profilename1);
	else reporterr("Writing alignment into %s and %s ...\n", profilename1, profilename2);
	prof1 = fopen(profilename1, "wb");
	if(prof1 == NULL) { reporterr("ERROR: can not write answer into %s.\n", profilename1); exit(1); }
	writeData_pointer(prof1, f1seq, name, nlen, seq);
	if(! print_to_two_files) writeData_pointer(prof1, f2seq, name2, nlen22, seq2);
	fclose(prof1);
	prof2 = fopen(profilename2, "wb");
	if(prof1 == NULL) { reporterr("ERROR: can not write answer into %s.\n", profilename2); exit(1); }
	if(print_to_two_files) writeData_pointer(prof2, f2seq, name2, nlen22, seq2);
	fclose(prof2);
	FreeDoubleVec(eff);
	FreeDoubleVec(eff2);

	reporterr("\ndone. \n");

#if REPORTCOSTS
//		use_getrusage();
		reporterr( "\n2-profiles align, real = %f min\n", (float)(time(NULL) - starttime)/60.0 );
		reporterr( "2-profiles align, user = %f min\n", (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
	if(spscoreout)
	{
		if(seq) FreeCharMtx(seq);
		if(name) FreeCharMtx(name);
		if(nlen) FreeIntVec(nlen);
		prof1 = fopen(profilename1, "rb");
		getnumlen(prof1); rewind(prof1);
		seq = AllocateCharMtx(njob, nlenmax + 10);
		name = AllocateCharMtx(njob, BLEN + 10);
		nlen = AllocateIntVec(njob);
		readData_pointer(prof1, name, nlen, seq);
		fclose(prof1);
		reporterr("SP Scores = %.6f\n", sumofpairsscore(njob, seq));
		FreeIntVec(nlen);
		FreeCharMtx(name);
		FreeCharMtx(seq);
	}
	reporterr( "%s (%s, %d-bit) Version" , progName( argv[0] ), (dorp=='d')?"nuc":((nblosum==-2)?"text":"aa"), sizeof(int *) * 8 );
	reporterr( "%d.%d.%d.%d", VER_MAJOR, VER_MINOR, VER_RELEASE_PROF, VER_BUILD );
	reporterr( "%s \nalg=%c, model=%s, amax=%3.1f\n%d thread\n\n", VERSION, alg, modelname, specificityconsideration, nthread );
	return 0;
}
