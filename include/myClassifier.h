/* - note to version 0.1:	* a set of tools to determine protein CATH classes. */
	

#ifndef _CATH_CLASSIFIER
#define _CATH_CLASSIFIER

#include "define.h"
#include "aatype.h"
#include "toolz.h"
#include "LLt.h"

using namespace std;


class myClassifier
{
	public:
	
		bool isInit;
		bool isUsable;
		uint *partCount;
	uint *cclass;
	uint classesCount;
	uint featsCount;
	double **statsVector;
	double ***statsMatrix;
	double ***featsVector;
	double **featsWeightVector;
	double **partWeight;
	double lambda;
	string prompt;
	string cfgfilename;
	string solutionfilename;
	
	resultsTable *results;

	/** (multi)-classifier constructor
	@param NC number of classes: should be set to 2 for yes/no type of classification;
	@param nsamples number of samples from the sample space;
	@param SC number of features;
	@param myprompt output messages prompt. */
	myClassifier(int NC, int nsamples, int SC, string myprompt);

	/** standard (multi)-classifier destructor */
	~myClassifier();

	/** (multi)-classifier initializer
	@param nsamples number of samples from the sample space. It's used here as a check value. */
	void init(int nsamples);

	/** builds up statistics and solves linear equations systems */
	void buildScoringFunctions();

	/** adds an instance of a certain feature to build up statistics for a certain class
	@param cc current class of interest;
	@param featID feature coordinate in the feature vector;
	@param feat actual instance value of the feature to be added. */
	void addFeatureInstance(int cc, uint featID, double feat);

	/** reads solution parameter vectors */
	void readSolutions();

	/** writes solution parameter vectors */
	void writeSolutions();

	/** increments of one the progressive sample coordinate identifier for a certain class
	@param cc current class of interest. */
	void incrementSampleID(int cc);

	/** updates the results table
	@param myClass predicted class;
	@param realClass actual class. */
	void updateTable(int myClass, int realClass);

	/** gets a prediction of the class, given a certain feature vector
	@param targetVector array of features.
	@return predicted class. */
	int getClass(double *targetVector);

	/** gets the score for a certain class, given a certain feature vector
	@param cclass class we want the score for;
	@param targetVector array of features.
	@return score for the given class. */
	double getClassScore(int cclass, double *targetVector);

	/** computes the averages that enter the linear equations system coefficient matrix
	@param feats array of features;
	@param NF number of features;
	@param size size of the Sample;
	@param Nw statistic matrix;
	@param Qw statistic vector;
	@param Nlin number of linear terms;
	@param lin tells whether to use linear terms or not;
	@param pow tells whether to use quadratic terms or not. */
	static int learnStats(double **feats, uint NF, uint size, double **(&Nw), double *(&Qw), uint Nlin = 0, bool lin = true, bool pow = false);
	
	protected:

		uint *sampleID;

		/** reads parameters from a configuration file */
		void readParameters();
};

#endif

