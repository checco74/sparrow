/* toolz.h - version 1.7prof
some tools for the LSM programs
- note to version 1.7prof:	* getScoreFactors() was changed;
	* different neighbour numbers for auxiliary functions. */

#ifndef _TOOLZ
#define _TOOLZ

#include <string>

#include "define.h"
#include "aatype.h"
#include "scoreGear.h"

// see parameters.h
extern detailLevel detail;
extern string encoding;
extern unsigned int SPARROW_WINDOW_LENGTH;
extern unsigned int SPARROW_SCORE_TYPES;
extern unsigned int SPARROW_NEIGHBOURS;
extern int subdim;
extern int StructCount;
extern int totalStructCount;
extern int patternCount;
extern int pairsCount;
extern unsigned int setSize;
extern unsigned int xsetSize;
extern unsigned int scorePlotSize;
extern unsigned int bufferSize;
extern unsigned int fieldSize;
extern int **neighboursCount;
extern int **auxnneighboursCount;
extern int **auxaneighboursCount;
extern int zoneCount;
extern int CtermSize;
extern int NtermSize;
extern int indent;
extern int *minlength;
extern int gapsize;
extern int *shift;
extern int *powFactor;
extern bool isotropy;
extern bool transcoding;
extern bool preprocessing;
extern bool postprocessing;
extern double thresh;
extern double safeval;
extern double safediff;
extern double criticval;
extern double lambdalin;
extern double lambdaquad;
extern double lambdaquadmax;
extern double lambdaDeltaSquare;
extern double *super_reweight;
extern double avwgh;
extern double normpower;
extern double devpower;
extern double oneoverdevpower;
extern double confidenceStep;
extern unsigned int confidenceMaxLevel;
extern unsigned int linewidth;
extern string tablename;
extern string servername;
extern string proteinfilename;
extern string predfilename;
extern string netinputfilename;
extern string pipsfile;
extern string logfile;
extern string separator;
extern string commentchars;
extern char seqbreakchar;
// score evaluation parameters
extern bool scorePlot;
extern int EvaRandomSeed;
extern int chainsCount;
extern int verbosity;

static const short ENCODE_RANGE = 20;
static const short OFFSET = 70;

/** @brief containes the specifics of a <code>Sample</code> type */
class SampleType
{
	public:
		/** @brief window length */
	uint length;
	/** @brief window key position */
	uint key;

	/** @brief initializes the values <code>length</code> and <code>key</code> */
	void init(int L = 1);
};

/** @brief holds the results for a single class */
class resultsTable
{
	public:
		/** @brief total occurences of the class */
	int cnt;
	/** @brief total predictions of the class */
	int mycnt;
	/** @brief total correct predictions of the class */
	int good;
	/** @brief total missed predictions of the class */
	int miss;
	/** @brief column of the confusion matrix per amino acid type, corresponding to the pertaining class */
	unsigned long **confusionVector;

	resultsTable();

	~resultsTable();

	/** @brief allocates memory for the confusion vector */
	void init(int classesCount);
};

/** erases some working arrays */
void eraseStuff();

/** checks the presence of working directories */
void checkDirectory();

/** @return number of different secondary structure types depending on <code>detailLevel</code> */
int getStructCount();

/** projects the secondary structure raw code, as given by <code>DBInterface</code>
@param p secondary structure raw <code>StructureType</code> code;
@param definitionLevel secondary structure level of definition.
@return secondary structure code depending on <code>detailLevel</code> */
int project(StructureType p, int definitionLevel = -1);

/** @return code corresponding to random coil depending on <code>detailLevel</code> */
int coil();

/** projects the secondary structure code into the Q3 numerical code
@param p secondary structure's general numerical code (@see <code>detailLevel</code>);
@param level level of detail in secondary structure definition (@see <code>detailLevel</code>).
@return numerical Q3-code: 0 for helix, 1 for strand, 2 for coil */
Category getQ3(int p, detailLevel level);

/** projects any alphabetical code into the Q3 alphabetical code
@param c secondary structure's general alphabetical code (@see <code>detailLevel</code>).
@return alphabetical Q3-code: 'H' for helix, 'e' for strand, '.' for coil */
char sstrucQ3(char c);

/** transforms the secondary structure code into an alphabetical code
@param p secondary structure's general numerical code (@see <code>detailLevel</code>);
@param definitionLevel secondary structure level of definition.
@return alphabetical code */
char sstruc2char(int p, int definitionLevel = -1);

/** converts chain zone code from alphabetical to numerical
@param c chain zone alphabetical code.
@return numerical chain zone code */
int zcode(char c);

/** converts chain zone code from numerical to alphabetical
@param z chain zone numerical code.
@return alphabetical chain zone code */
char zname(int z);

/** sets the target secondary structure pattern (specific or general)
@param argc number of command line arguments;
@param argv comman line arguments;
@param region protein zone of interest: is set by this function;
@param antipattern eventual antagonist secondary structure pattern in 1vs1-modes: also set by this function.
@return target secondary structure pattern or a dummy value */
int init_pattern(int argc, char *argv[], int &region, int &antipattern);

/** computes the pair identifier, given two secondary structure motifs identifiers
@param c0 first secondary structure pattern identifier;
@param c1 second secondary structure pattern identifier.
@return secondary structure pair ID. */
int getMotifsPair(int c0 = -1, int c1 = -1);

/** sets the values of the shifts, depending on <code>gapsize</code>, from the centre of the secondary structure pattern */
void setShifts();

/** sets the values of the power base units, depending on <code>subdim</code> */
void setPowerBaseUnits();

/** constructs a filename reflecting the current settings
@param name filename holding array;
@param ext extension and (hidden) directory;
@param sm scoring function traits container;
@param d secondary structure pattern dimension;
@param zone protein region;
@param dir direction along the protein chain;
@param opt decides which type of scoring function is to be assumed.
@return 1 if any errors occured while trying to open the file with the corresponding name 0 otherwise. */
int buildFilename(char *name, char *ext, scoreMachine *sm, int d, int zone, int dir = 1, readMode opt = STD_MODE);

/** writes a text file containing a square matrix's float entries contained in a one-dimensional array
@param file file name;
@param Nw matrix array;
@param size number of rows. */
void writeMatrix(char *file, double *Nw, uint size);

/** writes a text file containing a square matrix's float entries contained in a two-dimensional array
@param file file name;
@param Nw matrix array;
@param size number of rows. */
void writeMatrix(char *file, double **Nw, uint size);

/** writes a text file containing a square matrix's integer entries contained in a two-dimensional array
@param file file name;
@param Nw matrix array;
@param size number of rows. */
void writeMatrix(char *file, int **Nw, uint size);

/** writes a text file containing a vector's entries
@param file file name;
@param x vector array;
@param size array's size. */
void writeStats(char *file, double *x, uint size);

/** writes a binary file containing a matrix's float entries followed by a long integer
@param file file name;
@param Nw matrix array;
@param size number of rows;
@param N long integer. */
void storeMatrix(char *file, double **Nw, unsigned int size, unsigned int &N);

/** writes a binary file containing a vector's float entries followed by a long integer
@param file file name;
@param x vector array;
@param size array's size;
@param N long integer. */
void storeStats(char *file, double *x, unsigned int size, unsigned int &N);

/** writes a binary file containing a vector's float entries followed by one more float
@param file file name;
@param b float;
@param w vector array;
@param size array's size */
void storeSolutions(const char *file, double b, double *w, unsigned int size);

/** retrieves a matrix's float entries from a binary file
@param file file name;
@param mtx will store the matrix's entries;
@param size number of rows;
@param N optional sequence count. */
int getMatrix(char *file, double **(&mtx), uint size, unsigned int *N = NULL);

/** retrieves a vector's float entries from a binary file
@param file file name;
@param x will store the vector's entries;
@param size vector's size;
@param N optional sequence count. */
int getStats(char *file, double *(&x), uint size, unsigned int *N = NULL);

/** retrieves a vector's float entries followed by one more float from a binary file
@param file file name;
@param b float;
@param w will store the vector's entries;
@param size vector's size. */
int getSolutions(const char *file, double &b, double *(&w), uint size);

/** retrieves support vector machine parameters (SVMlite) */
int getSVMSolutions(char *file, double *b, double *w, uint size);

/** secondary structure element length (alphabetical)
@param s secondary structure string;
@param k position along the chain.
@return computed length */
int segmentSize(const char *(&s), int k);

/** secondary structure element length (numerical)
@param s secondary structure code array;
@param length length of the code array;
@param k position along the chain.
@return computed length */
int segmentSize(int *(&s), int length, int k);

/** computes the effective number of scorers
@param pattern secondary structure pattern of interest;
@param NS array containg the partial number of scorers per secondary structure pair;
@param neighbours array containing the number of neighbours per secondary structure and secondary structure pair.
@return effective number of scorers. */
int computeEffectiveNS(int pattern, int *NS, int **neighbours);

/** computes the number of up to second order scorers combinations
@param totalNS number of scorers.
@return number of combinations. */
int computeNQuadratic(int totalNS);

/** counts the types of <code>Sample</code>s that are going to be needed out of the scorer types
@param sm array containing the needed <code>ScoreMahine</code>;
@param NS array containing the numbers of scorers;
@param ST array containing the needed <code>Sample</code> types.
@return number of <code>Sample</code> types needed. */
int countTypes(scoreMachine ***(&sm), int **NS, SampleType *(&ST));

/** checks whether a line conforms to the profile standard and returns the correspondingly formatted profile
@param line line to be matched against the format;
@return the formatted profile line on success, an empty string otherwise. */
string matchProfileFormat(string line);

/** transforms a raw profile string into an array of doubles
@param rawProfile original profile string;
@param profLength profile length;
@param dir direction along the protein chain;
@param pos starting position along the sequence.
@return a pointer to the array of doubles. */
double *getProfile(string &rawProfile, int profLength, int dir, int pos = 0);

/** collapses a raw profile string and transforms it into an array of doubles
@param sequence original profile string;
@param length chain length;
@param winLength sequence window length;
@param profLength profile length;
@param dir direction along the protein chain;
@param pos starting position along the sequence.
@return a pointer to the array of doubles. */
double *getCollapsedProfile(string &sequence, int length, int winLength, int profLength, int dir, int pos = 0);

/** computes a profile's entry out of the raw <code>DBInterface</code> entry
@param rawEntry encoded profile value.
@return double equivalent of the given character. */
double getProfileEntry(char rawEntry);

/** examines the partial motif for the NC direction, depending on previous predictions or initializations
@param prnID current secondary structure pattern under examination;
@param mymotif array of current secondary structure predictions.
@param pos chain target position
@return the one-residue SSE for the target position if the previous predictions
		are compatible with the current secondary structure pattern, -1 otherwise. */
int getNCPartial(int prnID, int *mymotif, int pos);

/** examines the partial motif for the CN direction, depending on previous predictions or initializations
@param prnID current secondary structure pattern under examination;
@param mymotif array of current secondary structure predictions.
@param pos chain target position
@param protlen chain length
@return the one-residue SSE for the target position if the previous predictions
		are compatible with the current secondary structure pattern, -1 otherwise. */
int getCNPartial(int prnID, int *mymotif, int pos, int protlen);

/** gives the residue position that marks the starting point of the protein profile portion to be examined
@param sm current <code>ScoreMachine</code>;
@param centre secondary structure pattern centre position;
@param dir direction along the protien chain.
@return starting position. */
int getStartingPosition(scoreMachine &sm, int centre, int dir);

/** initializes scores and temporary motif for a given residue
@param k residue position;
@param score array of pattern scores;
@param motif array of motifs. */
void initResidue(int k, double **(&score), int *(&motif));

/** computes the scores and adjusts the draft prediction array for the given residue position
@param k chain target residue position;
@param dir direction along the protein chain;
@param chain protein chain;
@param protlen protein length;
@param SG scoring function tools;
@param myscore array containing scores per position and secondary structure pattern;
@param mymotif array of draft secondary structure predictions;
@param amode analysis mode (if set to <code>NO_TRANSITIONS</code> the behaviour changes);
@param rmode determines the scorers used.
@return number of positively evaluated (@see <code>thresh</code>) secondary structure patterns. */
int scoreResidue(int k, int dir, Protein *chain, int protlen, scoreGear *(&SG), double **(&myscore), int *(&mymotif),
	int amode, readMode rmode);

/** normalizes the scores for list compilation
@param scores array of scores to be normalized;
@param n number of scores in the array. */
void normalize_scores(double *(&scores), int n);

/** arranges the "super"-scoring function input factors
@param prn secondary structure pattern;
@param NSSC number of super scoring function parameters;
@param partNS array containing the partial numbers of scorers;
@param totalNS total number of scorers;
@param sm array of scorers;
@param neighbours array containing the number of neighbours per secondary structure and secondary structure pair;
@param sequence amino acid sequence string;
@param profile profile string;
@param length protein chain length;
@param region protein region of interest;
@param centre secondary structure pattern centre position
@param dir direction along the protien chain.
@return array of raw scores. */
double *getScoreFactors(int prn, int NSSC, int *partNS, int totalNS, scoreMachine **(&sm), int **neighbours, string &sequence, string &profile,
	int length, int region, int centre, int dir);

/** computes a raw score to enter the super scoring functions
@param sm scoring function traits container;
@param sequence amino acid sequence string;
@param profile profile string;
@param length protein chain length;
@param profLength length of a profile;
@param region protein region of interest;
@param dir direction along the protien chain;
@param start starting residue position (NC-direction);
@return raw score pertaining to the given scoring function. */
double computeRawScore(scoreMachine &sm, string &sequence, string &profile, int length, int profLength, int region, int dir,
	int start = 0);

/** computes the actual super scoring function input factors and stores them in an array
@param rawscore array of raw scores;
@param NS total number of scorers;
@param NSSC number of super scoring function parameters;
@return array of super scoring function input factors. */
double *computeScoreFactors(double *(&rawscore), int NS, int NSSC);

/** computes the error the solutions to a linear equations system are affected by and frees debug arrays memory
@param dim number of linear equations;
@param dbgmtx matrix of coefficients;
@param dbgq vector of constant terms;
@param wv solution vector. */
void computeSolutionError(uint dim, double **(&dbgmtx), double *(&dbgq), double *wv);

/** updates the table of results for the single class prediction evaluation
@param k current sample id;
@param score score obtained by the current secondary structure pattern for the current sample;
@param rt table holding the results;
@param pattern current secondary structure pattern;
@param realpattern effective secondary structure pattern for the current sample;
@param bs file stream eventually receiving the output in case this is requested (@see <code>scorePlotSize</code>. */
void updateTable(uint k, double score, resultsTable &rt, int pattern, int realpattern, ofstream &bs);

/** writes a double to an output stream
@param os output stream;
@param v double number to be output;
@param c string containing the separator for the decimals. */
void writeScore(ostream &os, double v, const char *c = ".");

/** writes to a file the results for the single class prediction evaluation
@param total number of samples in the data set;
@param rt table holding the results. */
void writeSummaryTable(int total, resultsTable &rt);

int** memresInt(long unsigned int dimen);

double** memresFloat(long unsigned int dimen);

double* memresVec(long unsigned int dimen);

void memrelease(int **(&v), long unsigned int dimen);

void memrelease(double **(&v), long unsigned int dimen);

void memrelease(double *(&v), long unsigned int dimen);

void reInit(double **(&v), long unsigned int dim);

void reInit(double *(&v), long unsigned int dim);

void writeVector(const char *filename, int *v, uint size, char format = ' ');

uint readVector(const char *filename, int *v, uint size);

void writeVector(const char *filename, double *v, uint size, char format = ' ');

uint readVector(const char *filename, double *v, uint size);

void copyArray(double **sourceArray, double **(&destinationArray), int dim);

void copyArray(double *sourceArray, double *(&destinationArray), int dim);

void writeHitCounter(uint **array);

void writeAminoAcidsStats(resultsTable *(&rt));

string encodeProfile(string &profline);

/** reads a <code>Protein</code> object from the input files provided
@param imode input type;
@param proteinfilename file containing the amino acid sequence;
@param profilename file containing the sequence profile;
@return pointer to the <code>Protein</code> object. */
Protein *readChain(inputMode imode, string proteinfilename, string profilename);

/** computes the factorial.
@param number factorial base. */
int factorial(int number);

/** computes the binomial coefficient.
@param n total number of elements;
@param k number of elements per grouping. */
int binomialCoefficient(int n, int k);

/** computes the binomial probability distribution factor.
@param n total number of elements;
@param k number of elements per grouping;
@param p base probability. */
double binomialProbability(int n, int k, double p);

/** computes the determinant */
double determinant(double **a, int n);

/** computes the generalized Matthews correlation coefficient, given the results table */
double computeCorrelationCoefficient(resultsTable *(&rC), int catalogSize);

/** computes the generalized Matthews correlation coefficient, given the confusion matrix */
double computeCorrelationCoefficient(double **(&confusionMatrix), int catalogSize);

/** prints the confusion matrix */
void printConfusionMatrix(double **(&confusionMatrix), int catalogSize);

/** computes the total number of residues in an array of <code>Protein</code>s */
uint computeTotalResidues(Protein *(&ivec));

/** compares two double precision numbers;
@param a pointer to the first number;
@param b pointer to the second number.
@return the difference between the first and the second number. */
int comparefloat (const void *a, const void *b);

/** projects an input in [0,1] onto integers in [0,9];
@param input the double value to be transformed.
@return an unsigned integer between 0 and 9 proportional to the input; 0 if given negative numbers, 9 if given number greater than 1. */
uint fastaTransform(double input);

/**
 * @brief Encodes a number between -20 and 20 as a character between '2' and 'Z'.
 *
 * @param x A short integer between -20 and 20 that should be encoded
 * @return A character encoding for the given integer
 */
inline char TRencode(short x)
{
  assert(x >= -ENCODE_RANGE && x <= ENCODE_RANGE);
  char encoded = (char)(OFFSET + x);
  return encoded;
}

/**
 * @brief Decodes a character between '2' and 'Z' into a short integer between -20 and 20.
 *
 * @param x A character between '2' and 'Z' that should be decoded
 * @return A short integer representing the decoded character
 */
inline short TRdecode(char x)
{
  assert(x >= (OFFSET - ENCODE_RANGE) && x <= (OFFSET + ENCODE_RANGE));
  short decoded = (short)x - OFFSET;
  return decoded;
}

/** rounds off number to closest integer */
inline double roundoff(double number)
{
	return (int)(number+0.5);
}

inline bool isCommentLine ( string &buffer )
{
        size_t cpos = buffer.find_first_of ( commentchars );
        size_t pos = buffer.find_first_not_of ( commentchars );
        if ( pos == std::string::npos or cpos < pos ) return true;
        else return false;
}

void waitChild(int pid);

#endif

