/* SPARROW - version 1.0 - July 2009, Francesco Bettella. */

#define STAGE_I 0
#define STAGE_II 1
#define DEBUG_PROTEIN -1
#define VERSION "SPARROW v1.0"

#include "define.h"
#include "aatype.h"
#include "parameters.h"
#include "myClassifier.h"
#include "scoreGear.h"
#include "toolz.h"

using namespace std;

char rgbcolor[10][12] = {
 "1.0,0.0,1.0",
 "1.0,0.1,1.0",
 "0.9,0.2,0.9",
 "0.8,0.3,0.8",
 "0.5,0.4,0.5",
 "0.4,0.5,0.4",
 "0.4,0.6,0.4",
 "0.1,0.7,0.1",
 "0.1,0.8,0.1",
 "0.0,0.9,0.0"
};

class coresGrid
{
	int Nseeds;
	int *pipstart;
	int *pipmotif;
	coresGrid(int pstage = STAGE_I);
};


/** reads the list of starter amino acid positions */
void readStarters(int nprot, const char *filename, Protein *pept, int **(&mpips), int **(&pips), int *(&n));

/** writes the list of starter amino acid positions */
void writeStarters(uint nprot, const char *filename, Protein *pept, int **(&mpips), int **(&pips), int *(&n));

/** fixes eventual clashes in the raw predicted secondary structure */
int ambiguity(int *(&myguess), int protlen, int k);

/** fixes undersized portions of raw predicted secondary structures */
void adjustLength(int *(&g), double **score, int protlen, int k);

/** postprocessing tools container/manager */
void processChain(Protein *(&chain), int *(&myguess), double **myscore);

/** summarizes the material produced by analyseChain() */
void summarizeChain(Protein *(&chain), double **s, double *avs, myClassifier *CE, int *(&conflictType));

/** protein chain secondary structure prediction tool (standard) */
void analyseChain(Protein *(&chain), scoreGear *(&SG), myClassifier *CE, int amode, bool use_auxiliary);

/** protein chain prediction manager */
void chaingamble(scoreGear *(&SG), Protein *(&pept), uint n, int amode, bool use_aux);



int main(int argc, char *argv[])
{
	int dummy = 0;
	int amode = 0; // analysis mode: normal or no-transitions
	int zone = CENTER; // just another dummy variable in this module
	bool auxmode = false; // decides whether to use auxiliary scorers or not
	string profilefilename = ""; // holds the name of the file containing the psiblast profile
	inputMode imode = FASTA_MODE;	// input file type
	scoreGear *myGear = NULL;
	Protein *pept = NULL;
	time_t t0, t1;
	t0 = time(0);
	checkDirectory();
	aminoacid::init();
	cout.setf(ios::fixed, ios::floatfield); // floatfield set to fixed
	dummy = init_pattern(argc, argv, zone, dummy);
	for (int i = 1; i<argc; i++)
	{
		if (strncmp(argv[i], "-", 1)==0)
		{
			if (strncmp(argv[i], "-aux", strlen("-aux"))==0) auxmode = true;
			else if (strncmp(argv[i], "-notrans", strlen("-notrans"))==0) amode = NO_TRANSITIONS;
			else if (strncmp(argv[i], "-psiblast", strlen("-psiblast"))==0) imode = PSIBLAST_MODE;
			else if (strncmp(argv[i], "-fasta", strlen("-fasta"))==0) imode = FASTA_MODE;
			else if (strncmp(argv[i], "-verbosity=", strlen("-verbosity="))==0)
				sscanf(argv[i], "-verbosity=%d", &verbosity);
			else if (strncmp(argv[i], "-width=", strlen("-width="))==0)
				sscanf(argv[i], "-width=%u", &linewidth);
			else if (strncmp(argv[i], "-dsspmode=", strlen("-dsspmode="))==0)
			{
				string buffer = argv[i];
				encoding = buffer.substr(strlen("-dsspmode="));
				cout<<"DSSP reduction scheme set to: "<<encoding<<endl;
			}
			else cout<<"undefined option: "<<argv[i]<<endl;
		}
		else
		{
			proteinfilename = argv[i];
			profilefilename = proteinfilename;
			predfilename = proteinfilename + ".out";
			netinputfilename = proteinfilename + ".nni";
		}
	}
	if (imode==FASTA_MODE)
		profilefilename += ".psiblast";
	setShifts();
	setPowerBaseUnits();
	xsetSize = 1;
	chainsCount = xsetSize;
	pept = readChain(imode, proteinfilename, profilefilename);
	myGear = new scoreGear[patternCount];
	for (int p = 0; p<patternCount; p++)
	{
		myGear[p].Pattern = p;
		myGear[p].setup(auxmode);
	}
	if (verbosity>1)
	{
		cout<<"sequence:"<<endl<<pept->sequence<<endl;
		cout<<"profile:"<<endl<<pept->profile<<endl;
	}
	chaingamble(myGear, pept, xsetSize, amode, auxmode);
	if (pept!=NULL) delete pept;
	for (int p = 0; p<patternCount; p++)
	{
		myGear[p].eraseScoreTools();
		myGear[p].eraseSuperCoefficients();
	}
	delete[] myGear;
	aminoacid::kill();
	eraseStuff();
	t1 = time(0);
	return 0;
}




void chaingamble(scoreGear *(&SG), Protein *(&pept), uint nprot, int amode, bool use_auxiliary)
{
	int nr = 0;
	myClassifier *confidenceEstimator = NULL;
	confidenceEstimator = new myClassifier(2, computeTotalResidues(pept), StructCount, "myConfInd");
	remove("scores.dat");
	confidenceEstimator->readSolutions();
	cout<<"proceeding with the analysis...."<<endl;
	nr += pept->sequence.length();
	if (verbosity>1) cout<<"chain "<<"profile "<<pept->profile<<endl;
	analyseChain(pept, SG, confidenceEstimator, amode, use_auxiliary);
	cout<<"total number of residues: "<<nr<<endl;
	if (confidenceEstimator!=NULL) delete confidenceEstimator;
}




void analyseChain(Protein *(&chain), scoreGear *(&SG), myClassifier *CE, int amode, bool use_auxiliary)
{
	int k = 0, h = 0; // current positions
	int protlen = chain->sequence.length(); // length of the protein chain
	int *mymotif = new int[protlen]; // predicted motif per position
	int *conflictType = new int[protlen]; // stores the type of conflict the residue gives rise to (0 = no conflict)
	double **myscore = new double*[protlen]; // score per motif per position
	double *avscore = new double[StructCount]; // the average score for each secondary structure on the protein
	for (k = 0; k<protlen; k++)
	{
		conflictType[k] = 0;
		myscore[k] = new double[StructCount];
		initResidue(k, myscore, mymotif);
	}
	k = 0;
	h = k-1;
	if (k<indent)
		k = indent;
	if (h>protlen-indent-1)
		h = protlen-indent-1;
	if (h==k) h = k-1;
	while (k<protlen-indent)
	{
		if (verbosity>1) cout<<"position "<<k<<endl;
		int positivesCount = scoreResidue(k, 1, chain, protlen, SG, myscore, mymotif, amode, STD_MODE);
		if (positivesCount==0)
		{
			if (use_auxiliary)
				positivesCount = scoreResidue(k, 1, chain, protlen, SG, myscore, mymotif, amode, AUXN_MODE);
			conflictType[k] = -1;
		}
		else if (positivesCount>1)
		{
			if (use_auxiliary)
				positivesCount = scoreResidue(k, 1, chain, protlen, SG, myscore, mymotif, amode, AUXA_MODE);
			conflictType[k] = 1;
		}
		if (use_auxiliary && positivesCount!=1)
			scoreResidue(k, 1, chain, protlen, SG, myscore, mymotif, amode, STD_MODE);
		k++;
	}
	while (h>=0)
	{
		if (verbosity>1) cout<<"position "<<h<<endl;
		int positivesCount = scoreResidue(h, -1, chain, protlen, SG, myscore, mymotif, amode, STD_MODE);
		if (positivesCount==0)
		{
			if (use_auxiliary)
				positivesCount = scoreResidue(h, -1, chain, protlen, SG, myscore, mymotif, amode, AUXN_MODE);
			conflictType[h] = -1;
		}
		else if (positivesCount>1)
		{
			if (use_auxiliary)
				positivesCount = scoreResidue(h, -1, chain, protlen, SG, myscore, mymotif, amode, AUXA_MODE);
			conflictType[h] = 1;
		}
		if (use_auxiliary && positivesCount!=1)
			scoreResidue(h, -1, chain, protlen, SG, myscore, mymotif, amode, STD_MODE);
		h--;
	}
	summarizeChain(chain, myscore, avscore, CE, conflictType);
	delete[] conflictType;
	delete[] mymotif;
	for (k = 0; k<protlen; k++)
		delete[] myscore[k];
	delete[] myscore;
	delete[] avscore;
}




int ambiguity(int *(&g), int protlen, int k)
{
	if (k>0 && k<protlen)
	{
		if (detail<NO_REGIONS)
		{
			if ((g[k]<SHORT_STRAND || (g[k]>RANDOM_COIL && g[k]<TURN)) && (g[k-1]>=SHORT_STRAND && g[k-1]<RANDOM_COIL))
			{
				if (segmentSize(g, protlen, k)>segmentSize(g, protlen, k-1)) return 1;
				else return 0;
			}
			if ((g[k-1]<SHORT_STRAND || (g[k-1]>RANDOM_COIL && g[k-1]<TURN)) && (g[k]>=SHORT_STRAND && g[k]<RANDOM_COIL))
			{
				if (segmentSize(g, protlen, k-1)>segmentSize(g, protlen, k)) return 0;
				else return 1;
			}
		}
		else if (detail<DSSP_PLUS_SHORTS)
		{
			if ((g[k]<BETA_STRAND || (g[k]>BETA_BRIDGE && g[k]<GENERIC_COIL)) && (g[k-1]>=BETA_STRAND && g[k-1]<=BETA_BRIDGE))
			{
				if (segmentSize(g, protlen, k)>segmentSize(g, protlen, k-1)) return 1;
				else return 0;
			}
			if ((g[k-1]<BETA_STRAND || (g[k-1]>BETA_BRIDGE && g[k-1]<GENERIC_COIL)) && (g[k]>=BETA_STRAND && g[k]<=BETA_BRIDGE))
			{
				if (segmentSize(g, protlen, k-1)>segmentSize(g, protlen, k)) return 0;
				else return 1;
			}
		}
		else if (detail<DSSP_NO_EXOTIC)
		{
			if ((g[k]==DSSP_HELIX || g[k]==DSSP_T_HELIX || g[k]==DSSP_P_HELIX) && (g[k-1]==DSSP_STRAND || g[k-1]==DSSP_BRIDGE))
			{
				if (segmentSize(g, protlen, k)>segmentSize(g, protlen, k-1)) return 1;
				else return 0;
			}
			if ((g[k-1]==DSSP_HELIX || g[k-1]==DSSP_T_HELIX || g[k-1]==DSSP_P_HELIX) && (g[k]==DSSP_STRAND || g[k]==DSSP_BRIDGE))
			{
				if (segmentSize(g, protlen, k-1)>segmentSize(g, protlen, k)) return 0;
				else return 1;
			}
		}
		else
		{
			if (g[k]==_HELIX && g[k-1]==_STRAND)
			{
				if (segmentSize(g, protlen, k)>segmentSize(g, protlen, k-1)) return 1;
				else return 0;
			}
			if (g[k-1]==_HELIX && g[k]==_STRAND)
			{
				if (segmentSize(g, protlen, k-1)>segmentSize(g, protlen, k)) return 0;
				else return 1;
			}
		}
	}
	return -1;
}




void adjustLength(int *(&g), double **score, int protlen, int k)
{
	int i = k, j = k;
	int new_length = 1;	// holds the progressive residue count in the secondary structure element
	int current_length = 1;
	int *gp = new int[protlen];
	for (int h = 0; h<protlen; h++) gp[h] = g[h];
	bool *used = new bool[StructCount];
	for (int p = 0; p<StructCount; p++)
	{
		if (p==g[k]) used[p] = true;
		else used[p] = false;
	}
	do
	{
		if (new_length>current_length) current_length = new_length; // we succeeded in increasing the original motif's length
		else // if it's not possible to extend the motif in any direction we try by changing the motif
		{
			new_length = 0;
			for (int h = 0; h<protlen; h++) gp[h] = g[h];
			double maxscore = -11.;
			for (int p = 0; p<StructCount; p++)
			{
				if (!used[p])
				{
					if (score[k][p]>maxscore)
					{
						gp[k] = p;
						maxscore = score[k][p];
						current_length = segmentSize(gp, protlen, k);
						new_length = 1; // new_length starts always from 1
					}
				}
			}
			used[gp[k]] = true;
			i = k; j = k; // i and j get initialized here
		}
		while (--i>=0 && gp[k]==gp[i]) new_length++; // new_length gets incremented with the number of matching residues behind....
		if (i>=0 && new_length<=minlength[gp[k]]) // .... if necessary we try extending backwards.
		{
			gp[i] = gp[k];
			int s = segmentSize(gp, protlen, i-1);
			if (ambiguity(gp, protlen, i)>=0 || (s<=minlength[gp[i-1]] && s>=current_length))
			{
				gp[i] = g[i]; // if the change generates incompatible neighbours or affects existing significant motifs we restore the previous motif
				i = -1; // we set the index to the minimum in order to prevent further attempts in this direction
			}
			else new_length++; // otherwise the new_length of the motif actually increases of 1
		}
		while (++j<protlen && gp[k]==gp[j]) new_length++; // new_length gets incremented with the number of matching residues ahead....
		if (j<protlen && new_length<=minlength[gp[k]]) // .... if necessary we try extending forwards.
		{
			gp[j] = gp[k];
			int s = segmentSize(gp, protlen, j+1);
			if (ambiguity(gp, protlen, j+1)>=0 || (s<=minlength[gp[j+1]] && s>=current_length))
			{
				gp[j] = g[j]; // if the change generates incompatible neighbours or affects existing significant motifs we restore the previous motif
				j = protlen; // we set the index to the maximum in order to prevent further attempts in this direction
			}
			else new_length++; // otherwise the new_length of the motif actually increases of 1
		}
	}
	while (new_length>0 && new_length<=minlength[gp[k]]);
	for (int h = 0; h<protlen; h++) g[h] = gp[h];
	delete[] used;
	delete[] gp;
}




void processChain(Protein *(&chain), int *(&myguess), double **myscore)
{
	double tempscore = 0.;
	int protlen = chain->sequence.length(); // protein length
	for (int k = indent; k<protlen-indent; k++)
	{
		double goldscore = -11.;
		double silverscore = -11.;
		for (int p = 0; p<StructCount; p++)
		{
			tempscore = (avwgh*((k<protlen-1 && k>0) ? myscore[k+1][p]+myscore[k-1][p] : 0.) + myscore[k][p]) / (1.+2.*avwgh);
			if (tempscore>=silverscore)
			{
				if (tempscore>=goldscore)
				{
					myguess[k] = p;
					silverscore = goldscore;
					goldscore = tempscore;
				}
				else silverscore = tempscore;
			}
		}
	//	if the best score is below the safeval treshold or
	// the difference between the best and the second best is
	// below the safediff treshold the guess for the relative
	// residue position is set to the value of StructCount, which
	// causes the residue to be ignored in the final count.
		if (goldscore<safeval || fabs(goldscore-silverscore)<safediff) myguess[k] = StructCount;
	}
	if (postprocessing)
	{
		for (int k = indent; k<protlen-indent; k++)
		{
			int a = ambiguity(myguess, protlen, k);
			if (k>0 and a>=0) myguess[k-a] = coil();
		}
		for (int k = indent; k<protlen-indent; k++)
		if (myguess[k]!=StructCount && segmentSize(myguess, protlen, k)<=minlength[myguess[k]])
		myguess[k] = coil();
	//	adjustLength(myguess, myscore, protlen, k);
	}
}




void summarizeChain(Protein *(&chain), double **myscore, double *avscore, myClassifier *CE, int *(&conflictType))
{
	int protlen = chain->sequence.length(); // protein length
	int *myguess = new int[protlen];
	for (int k = 0; k<protlen; k++) myguess[k] = coil();
	double *confidenceFeats = new double[StructCount];
	fstream ps(predfilename.c_str(), ofstream::out);
	fstream ss(netinputfilename.c_str(), ofstream::out);
	if (!ps.is_open())
	{
		cout<<"\nError: couldn't open prediction file. Aborted."<<endl;
		exit(1);
	}
	if (!ss.is_open())
	{
		cout<<"\nError: couldn't open scores file. Aborted."<<endl;
		exit(1);
	}
	ps<<"sequence > "<<chain->sequence<<endl;
	ps<<"sparrow >> "<<flush;
	processChain(chain, myguess, myscore);
	ss<<"#domain "<<chain->sequence<<endl;
	for (int k = 0; k<protlen; k++)
	{
		for (int p = 0; p<StructCount; p++)
		ss<<setprecision(9)<<myscore[k][p]<<'\t';
		ss<<coil()<<'\t'<<k<<'\t'<<chain->sequence[k]<<endl;
		if (k<indent or k>=protlen-indent) ps<<"-";
		else ps<<sstruc2char(myguess[k]);
	}
	ps<<endl<<"confidence "<<flush;
	for (int k = 0; k<protlen; k++)
	{
		for (int p = 0; p<StructCount; p++)
			confidenceFeats[p] = myscore[k][p];
		double confidence = CE->getClassScore(1, confidenceFeats);
		ps<<fastaTransform(confidence);
	}
	ps<<endl;
	ps.close();
	ss.close();
	delete[] confidenceFeats;
	delete[] myguess;
}
