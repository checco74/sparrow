/* parameters.h - version 1.57prof
	- note to version 1.57prof:	* different neighbour numbers for auxiliary functions. */

#ifndef _PARAMETERS
#define _PARAMETERS
#include "define.h"
#include <string>

using namespace std;

detailLevel detail = Q3_STANDARD; // secondary structure definition level
unsigned int SPARROW_WINDOW_LENGTH = 15;
unsigned int SPARROW_SCORE_TYPES = 1;
unsigned int SPARROW_NEIGHBOURS = 7;
string encoding = "loose"; // DSSP reduction scheme
string tablename = "checco1"; // DBInterface table containing the target domains
string proteinfilename = "target.dat"; // domain coordinates in the DBInterface table
string predfilename = "mygamble.txt"; // secondary structure prediction output file
string netinputfilename = "myscores.dat"; // secondary structure scores to be input by the neural network
string logfile = "lsm.log"; // program messages output file
string pipsfile = "pipsfile.dat";
string servername = "coulomb";
string separator = "eu";
string commentchars = "#>"; // "^°$%&?@+*~#<>|-";
char seqbreakchar = '!';
int subdim = 1; // secondary structure pattern size (in number of residues involved)
int StructCount = NUM_OF_STRUCTURE_TYPES; // actual number of secondary structure types (single residue)
int totalStructCount = NUM_OF_STRUCTURE_TYPES; // total number of secondary structure types (single residue)
int patternCount = StructCount; // total number of secondary structue pattern types (multiple residues)
int pairsCount = -1; // number of secondary structure pattern pairs
unsigned int setSize = 0;
unsigned int xsetSize = 0;
unsigned int scorePlotSize = 0;
unsigned int bufferSize = 1;
unsigned int fieldSize = 1; // the field is a sort of buffer in which all possible sequences of interest can fit.
int **neighboursCount = NULL; // number of neighbours used in super scoring functions per secondary structure pattern and pair type.
int **auxnneighboursCount = NULL; // number of neighbours used in aux.n super scoring functions per secondary structure pattern andpair type.
int **auxaneighboursCount = NULL; // number of neighbours used in aux.a super scoring functions per secondary structure pattern and pair type.
int zoneCount = CENTER+1;
int CtermSize = 0;
int NtermSize = 0;
int indent = 0;
int *minlength = NULL;
int gapsize = 0; // the gap between secondary structure components in higher dimensional secondary structure patterns
int *shift = NULL; // relative position of the secondary structure components in higher dimensional secondary structure patterns
int *powFactor = NULL; // secondary structure algebra power base units
bool isotropy = true; // decides whether to differentiate NC from CN directions
bool transcoding = false;
bool preprocessing = false;
bool postprocessing = false;
double thresh = 0.5;
double safeval = -1.0;
double safediff = 0.0;
double criticval = 0.0;
double lambdalin = 0.00001;	// regularization parameter for linear terms
double lambdaquad = 0.00001;	// regularization parameter for quadratic terms
double lambdaquadmax = 0.1;	// gives the maximum lambdaquad after dsq reweighting
double lambdaDeltaSquare = 10.;	// gives the square half height width of the lambdaquad strength
double *super_reweight = NULL;	// structure reweighting for the super scoring functions
double avwgh = 0.0;	// used to average subsequent signals in the prediction routine.
double cheatthresh = 0.5;
double normpower = 0.;	// to renormalize the parameter vectors
double devpower = 1.;	// used in the computation of the zeta score....
double oneoverdevpower = 1./devpower;	// ....also used in the computation of the zeta score
double confidenceStep = 0.1; // confidence interval for the confidence/accuracy correlation arrays
unsigned int confidenceMaxLevel = (unsigned int)(1.5/confidenceStep); // maximum confidence level for the confidence/accuracy correlation arrays
unsigned int linewidth = 40; // number of residues per line of output
// Additional parameters for score evaluation
bool scorePlot = true;
int EvaRandomSeed = 248;
int chainsCount = 10;
int verbosity = 0;
string LaTeXHeader = "\
\\documentclass[twoside]{report}\n\
\\usepackage[T1]{fontenc}\n\
\\usepackage[latin1]{inputenc}\n\
\\usepackage{courier}\n\
\\usepackage{dsfont}\n\
\\usepackage{fullpage}\n\
\\usepackage{rotating}\n\
\\usepackage{amsmath}\n\
\\usepackage{amssymb}\n\
\\usepackage{tabls}\n\
\\usepackage{multirow}\n\
\\usepackage{graphicx}\n\
\\usepackage{wrapfig}\n\
\\usepackage{subfig}\n\
\\usepackage{color}\n\
\\makeatletter\n\
\\begin{document}\
";

#endif
