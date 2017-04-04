/* computeConfInd.C - version 1.0
	This program computes the confidence indicators out of a list of scores. */

#define VERSION "computeConfInd - version 1.0"

#include "define.h"
#include "aatype.h"
#include "parameters.h"
#include "myClassifier.h"
#include "toolz.h"

/** reads the guess classes (correct or incorrect) based on the scores and the actual structure types */
bool **readClasses ( ifstream &is, unsigned int nd, unsigned int *nr );

/** reads the scores and stores them into a three-dimensional array */
double ***readScores ( ifstream &is, unsigned int &nd, unsigned int *(&nr) );

/** computes the confidence applying a magic formula to an array of scores */
double getConfidence ( double *myScores );

/** deduces the predictor's suggestion from an array of scores */
unsigned int deduceGuess ( double *myScores );

/** computes the feature ID based on the given parameters */
int computeID ( int t, unsigned short int window, int sstruct );

int main (int argc, char *argv[] )
{
	ifstream is;
	enum program_mode
	{
		LEARN_RULES,
		APPLY_RULES,
		APPLY_FORMULA
	};
	bool sortmode = false;
	std::string filename = "";
	program_mode mode = APPLY_RULES;
	unsigned int sampleSize = 0;
	unsigned int N_domains = 0;
	unsigned int *N_residues = NULL;
	unsigned short int window = 1;
	double ***myScores = NULL;
	double defaultfeat = -1.;
	double totalconf = 0.;
	double effconf = 0.;
	StructCount = getStructCount();
	std::cout.setf ( ios::fixed, ios::floatfield ); // floatfield set to fixed
	for ( int i = 1; i < argc; i++ )
	{
		if ( strncmp ( argv[i], "-learn", strlen ( "-learn" ) ) == 0 )
		{
			mode = LEARN_RULES;
			std::cout << "program set to learning mode." << std::endl;
		}
		else if ( strncmp ( argv[i], "-mf", strlen ( "-mf" ) ) == 0 )
		{
			mode = APPLY_FORMULA;
			std::cout << "program set to magic formula mode." << std::endl;
		}
		else if ( strncmp ( argv[i], "-sort", strlen ( "-sort" ) ) == 0 )
		{
			sortmode = true;
			std::cout << "sorting mode enabled." << std::endl;
		}
		else if ( strncmp ( argv[i], "-l=", strlen ( "-l=" ) ) == 0 )
		{
			sscanf ( argv[i], "-l=%hu", &window );
			std::cout << "window length set to: " << window << std::endl;
		}
		else if ( strncmp ( argv[i], "-verbosity=", strlen ( "-verbosity=" ) ) == 0 )
		{
			sscanf ( argv[i], "-verbosity=%d", &verbosity );
			std::cout << "verbosity set to level: " << verbosity << std::endl;
		}
		else if ( !is.is_open() )
		{
			is.open ( argv[i] );
			if ( !is.is_open() )
				std::cout << "couldn't open file " << argv[i] << std::endl;
			else
			{
				filename = argv[i];
				filename += ".cnf"; // filename is used to store output file name
				std::cout << "reading scores...." << std::endl;
				myScores = readScores ( is, N_domains, N_residues );
				for ( unsigned int k = 0; k < N_domains; k++ )
					sampleSize += N_residues[k];
				std::cout << "total domain chains: " << N_domains << std::endl;
				std::cout << "total residues: " << sampleSize << std::endl;
			}
		}
	}
	if ( is.is_open() )
	{
		stringstream errorMessages ( "" );
		unsigned int NFeats = window * StructCount;
		bool **myGuessClass = NULL;
		myGuessClass = readClasses ( is, N_domains, N_residues );
		myClassifier *confidenceEstimator = NULL;
		confidenceEstimator = new myClassifier ( 2, sampleSize, NFeats, ".conf" );
		switch ( mode )
		{
			case LEARN_RULES:
			{
				std::cout << "allocating " << sizeof(double)*sampleSize << " bytes for the guess classes...." << std::endl;
				for ( unsigned int k = 0; k < N_domains; k++ )
					for ( unsigned int h = 0; h < N_residues[k]; h++ )
						confidenceEstimator->partCount[myGuessClass[k][h]]++;
				confidenceEstimator->init(sampleSize);
				for ( unsigned int k = 0; k < N_domains; k++ )
				{
					for ( unsigned int h = 0; h < N_residues[k]; h++ )
					{
						if ( sortmode && window==1 )
							qsort ( myScores[k][h], (size_t)StructCount, sizeof(double), comparefloat );
						for ( int t = -((int)(window/2)); t <= (int)(window/2); t++ )
						{
							for ( int p = 0; p < StructCount; p++ )
							{
								int featID = computeID ( t, window, p );
								if ( featID >= 0 and featID < (int)confidenceEstimator->featsCount )
								{
									if ( h + t < N_residues[k] and h + t >= 0 )
										confidenceEstimator->addFeatureInstance ( myGuessClass[k][h], featID, myScores[k][h + t][p] );
									else confidenceEstimator->addFeatureInstance ( myGuessClass[k][h], featID, defaultfeat );
								}
								else
								{
									std::cout << "invalid feature ID: " << featID << std::endl;
									exit(1);
								}
							}
						}
						confidenceEstimator->incrementSampleID ( myGuessClass[k][h] );
					}
				}
				if ( sampleSize > confidenceEstimator->featsCount )
				{
					confidenceEstimator->buildScoringFunctions();
					confidenceEstimator->writeSolutions();
				}
				else std::cout << "\nWarning: not enough samples to build a confidence estimator" << std::endl;
				break;
			}
			case APPLY_RULES:
			{
				double *feats = NULL;
				feats = new double[NFeats];
				for ( unsigned int featID = 0; featID < NFeats; featID++ ) feats[featID] = 0.;
				unsigned int **confidenceHitCounter = NULL;
				confidenceHitCounter = new unsigned int*[3*StructCount];
				for(int p = 0; p < 3*StructCount; p++)
				{
					confidenceHitCounter[p] = NULL;
					confidenceHitCounter[p] = new unsigned int[confidenceMaxLevel+1];
					for(unsigned int clevel = 0; clevel < confidenceMaxLevel+1; clevel++)
						confidenceHitCounter[p][clevel] = 0;
				}
				confidenceEstimator->readSolutions();
				ofstream cfile ( filename.c_str() );
				for ( unsigned int k = 0; k < N_domains; k++ )
				{
					cfile << "#domain chain " << k << std::endl << "confidence " << std::flush;
					for ( unsigned int h = 0; h < N_residues[k]; h++ )
					{
						double *backupScores = NULL;
						backupScores = new double[StructCount];
						for ( int p = 0; p < StructCount; p++ )
							backupScores[p] = myScores[k][h][p];
						if ( sortmode && window==1 )
							qsort ( myScores[k][h], (size_t)StructCount, sizeof(double), comparefloat );
						for ( int t = -((int)(window/2)); t <= (int)(window/2); t++ )
						{
							for ( int p = 0; p < StructCount; p++ )
							{
								int featID = computeID ( t, window, p );
								if ( featID >= 0 and featID < (int)confidenceEstimator->featsCount )
								{
									if ( h + t < N_residues[k] and h + t >= 0 )
										feats[featID] = myScores[k][h + t][p];
									else feats[featID] = defaultfeat;
								}
								else
								{
									std::cout << "invalid feature ID: " << featID << std::endl;
									exit(1);
								}
							}
						}
						double confidence = confidenceEstimator->getClassScore ( 1, feats );
						effconf += myGuessClass[k][h] ? confidence : -confidence;
						totalconf += confidence;
						if ( verbosity > 0 )
						{
							std::cout << "residue " << h << " of domain chain " << k << std::flush;
							std::cout << " predicted with confidence " << confidence << std::flush;
							std::cout << " - scores: " << std::flush;
							for ( int p = 0; p < StructCount; p++ )
								std::cout << backupScores[p] << " " << std::flush;
							std::cout << std::endl;
						}
						cfile << fastaTransform(confidence) << std::flush;
						unsigned int myGuess = deduceGuess ( backupScores );
						int confidenceLevel = (int)(confidence/confidenceStep);
						if(confidenceLevel>(int)confidenceMaxLevel) confidenceLevel = confidenceMaxLevel;
						if(confidenceLevel<0) confidenceLevel = 0;
						if(myGuessClass[k][h]) confidenceHitCounter[myGuess][confidenceLevel]++;
						confidenceHitCounter[StructCount + myGuess][confidenceLevel]++;
						if ( backupScores != NULL ) delete[] backupScores;
					}
					if ( verbosity > 0 )
						std::cout << std::endl;
					cfile << std::endl;
				}
				writeHitCounter(confidenceHitCounter);
				std::cout << "cumulative confidence: " << effconf / totalconf << std::endl;
				for(int p = 0; p<3*StructCount; p++)
					if(confidenceHitCounter[p]!=NULL) delete[] confidenceHitCounter[p];
				if(confidenceHitCounter!=NULL) delete[] confidenceHitCounter;
				if ( cfile.is_open() ) cfile.close();
				if ( feats != NULL ) delete[] feats;
				break;
			}
			case APPLY_FORMULA:
			{
				unsigned int **confidenceHitCounter = NULL;
				double maxconfidence = -11., minconfidence = 11.;
				confidenceHitCounter = new unsigned int*[3*StructCount];
				for(int p = 0; p < 3*StructCount; p++)
				{
					confidenceHitCounter[p] = NULL;
					confidenceHitCounter[p] = new unsigned int[confidenceMaxLevel+1];
					for(unsigned int clevel = 0; clevel < confidenceMaxLevel+1; clevel++)
						confidenceHitCounter[p][clevel] = 0;
				}
				ofstream cfile ( filename.c_str() );
				for ( unsigned int k = 0; k < N_domains; k++ )
				{
					for ( unsigned int h = 0; h < N_residues[k]; h++ )
					{
						double confidence = getConfidence ( myScores[k][h] );
						if ( confidence > maxconfidence ) maxconfidence = confidence;
						if ( confidence < minconfidence ) minconfidence = confidence;
					}
				}
				for ( unsigned int k = 0; k < N_domains; k++ )
				{
					cfile << "#domain chain " << k << std::endl << "conf > " << std::flush;
					for ( unsigned int h = 0; h < N_residues[k]; h++ )
					{
						double confidence = getConfidence ( myScores[k][h] );
						confidence = ( confidence - minconfidence ) / ( maxconfidence - minconfidence );
						double cc = myGuessClass[k][h] ? confidence : -confidence;
						effconf += cc;
						totalconf += confidence;
						cfile << fastaTransform(confidence) << std::flush;
						unsigned int myGuess = deduceGuess ( myScores[k][h] );
						int confidenceLevel = (int)(confidence/confidenceStep);
						if(confidenceLevel>(int)confidenceMaxLevel) confidenceLevel = confidenceMaxLevel;
						if(confidenceLevel<0) confidenceLevel = 0;
						if(myGuessClass[k][h]) confidenceHitCounter[myGuess][confidenceLevel]++;
						confidenceHitCounter[StructCount + myGuess][confidenceLevel]++;
						if ( verbosity > 0 )
						{
							std::cout << "residue " << h << " of domain chain " << k << " predicted " << std::flush;
							std::cout << myGuess << "(" << myGuessClass[k][h] << ") with confidence " << confidence << std::endl;
							std::cout << "CRAI contribution: " << cc << std::endl;
						}
					}
					if ( verbosity > 0 )
						std::cout << std::endl;
					cfile << std::endl;
				}
				writeHitCounter(confidenceHitCounter);
				std::cout << "cumulative confidence: " << effconf / totalconf << std::endl;
				for(int p = 0; p<3*StructCount; p++)
					if(confidenceHitCounter[p]!=NULL) delete[] confidenceHitCounter[p];
				if(confidenceHitCounter!=NULL) delete[] confidenceHitCounter;
				if ( cfile.is_open() ) cfile.close();
				break;
			}
			default:
				errorMessages << "unknown program mode: " << mode << std::endl;
				std::cout << ( errorMessages.str().c_str() );
				exit ( 1 );
		}
		for ( unsigned int k = 0; k < N_domains; k++ )
		{
			if ( myGuessClass[k] != NULL ) delete[] myGuessClass[k];
			for ( unsigned int h = 0; h < N_residues[k]; h++ )
				if ( myScores[k][h] != NULL ) delete[] myScores[k][h];
			if ( myScores[k] != NULL ) delete[] myScores[k];
		}
		if ( myScores != NULL ) delete[] myScores;
		if ( myGuessClass != NULL ) delete[] myGuessClass;
		if ( N_residues != NULL ) delete[] N_residues;
		if ( confidenceEstimator != NULL )
			delete confidenceEstimator;
	}
	return 0;
}



bool **readClasses ( ifstream &is, unsigned int nd, unsigned int *nr )
{
	std::string buffer;
	unsigned int k = 0;
	unsigned int h = 0;
	unsigned int n_of_CommentLines = 0;
	unsigned int good = 0;
	unsigned int total = 0;
	bool **myGuessClass = NULL;
	double mean = 0.;
	double stdev = 0.;
	myGuessClass = new bool*[nd];
	for ( k = 0; k < nd; k++ )
	{
		myGuessClass[k] = new bool[nr[k]];
		for ( h = 0; h < nr[k]; h++ )
			myGuessClass[k][h] = false;
	}
	k = 0;
	h = 0;
	is.clear();
	is.seekg ( 0 );
	while ( getline ( is, buffer ) )
	{
		if ( !isCommentLine ( buffer ) )
		{
			double tmp = 0., max = -11.;
			unsigned int index = 0;
			unsigned int maxindex = (unsigned int)StructCount - 1;
			stringstream linebuffer ( buffer );
			while ( linebuffer >> tmp )
			{
				if ( index < (unsigned int)StructCount and tmp > max )
				{
					max = tmp;
					maxindex = index;
				}
				index++;
			}
			if ( k < nd and h < nr[k] )
			{
				total++;
				// the last entry held by tmp is the the real structure code
				if ( maxindex == (unsigned int)tmp )
				{
					myGuessClass[k][h] = true; // correct predictions
					good++;
				}
				else myGuessClass[k][h] = false; // incorrect predictions
				if ( verbosity > 1 )
				{
					std::cout << "chain " << k << ", residue " << h;
					std::cout << " - predicted: " << maxindex << "; actual: " << (unsigned int)tmp;
					std::cout << " --> [" << myGuessClass[k][h] << "]" << std::endl;
				}
			}
			else
			{
				std::cout << "readClasses: domain indices out of bounds --" << std::flush;
				std::cout << " domain ID = " << k << "[" << nd << "]," << std::flush;
				std::cout << " residue ID = " << h << "[" << nr[k] << "]" << std::endl;
				exit ( 1 );
			}
			h++;
		}
		else
		{
			if ( n_of_CommentLines > 0 )
			{
				k++;
				h = 0;
			}
			n_of_CommentLines++;
		}
	}
	mean = double(good) / double(total);
	for ( k = 0; k < nd; k++ )
	{
		for ( h = 0; h < nr[k]; h++ )
		{
			double val = (double)(myGuessClass[k][h]);
			stdev += ( val - mean ) * ( val - mean );
		}
	}
	stdev /= (double)total;
	stdev = sqrt ( stdev );
	std::cout << "score file report:" << std::endl;
	std::cout << " mean: " << mean << std::endl;
	std::cout << " stdev: " << stdev << std::endl;
	return myGuessClass;
}


double ***readScores ( ifstream &is, unsigned int &nd, unsigned int *(&nr) )
{
	std::string buffer;
	unsigned int k = 0;
	unsigned int h = 0;
	unsigned int n_of_CommentLines = 0;
	double ***myScores = NULL;
	nd = 0;
	nr = NULL;
	is.clear();
	is.seekg ( 0 );
	while ( getline ( is, buffer ) )
	{
		if ( isCommentLine ( buffer ) )
		{
			if ( verbosity > 1 )
				std::cout << "added chain " << nd << ": " << buffer << std::endl;
			nd++;
		}
	}
	nr = new unsigned int[nd];
	for ( k = 0; k < nd; k++ ) nr[k] = 0;
	k = 0;
	h = 0;
	is.clear();
	is.seekg ( 0 );
	n_of_CommentLines = 0;
	while ( getline ( is, buffer ) )
	{
		if ( isCommentLine ( buffer ) )
		{
			if ( n_of_CommentLines > 0 )
			{
				if ( k < nd - 1 ) // k is incremented one time less than n_of_CommentLines
				{
					if ( verbosity > 1)
						std::cout << "total residues in domain chain " << k << ": " << h << std::endl;
					nr[k] = h;
					h = 0;
					k++;
				}
				else
				{
					std::cout << "readScores[count]: domain index out of bounds --" << std::flush;
					std::cout << " domain ID = " << k << "[" << nd - 1 << "]" << std::endl;
					exit ( 1 );
				}
			}
			n_of_CommentLines++;
		}
		else h++;
	}
	if ( verbosity > 1)
		std::cout << "total residues in domain chain " << nd - 1 << ": " << h << std::endl;
	nr[nd - 1] = h;
	myScores = new double**[nd];
	for ( k = 0; k < nd; k++ )
	{
		myScores[k] = NULL;
		myScores[k] = new double*[nr[k]];
		for ( h = 0; h < nr[k]; h++ )
		{
			myScores[k][h] = NULL;
			myScores[k][h] = new double[StructCount];
			for ( int p = 0; p < StructCount; p++ )
				myScores[k][h][p] = 0.;
		}
	}
	k = 0;
	h = 0;
	is.clear();
	is.seekg ( 0 );
	n_of_CommentLines = 0;
	while ( getline ( is, buffer ) )
	{
		if ( isCommentLine ( buffer ) )
		{
			if ( n_of_CommentLines > 0 )
			{
				if ( k < nd - 1 ) // k is incremented one time less than n_of_CommentLines
				{
					h = 0;
					k++;
				}
				else
				{
					std::cout << "readScores[check]: domain index out of bounds --" << std::flush;
					std::cout << " domain ID = " << k << "[" << nd - 1 << "]" << std::endl;
					exit ( 1 );
				}
			}
			n_of_CommentLines++;
		}
		else
		{
			double tmp = 0.;
			unsigned int index = 0;
			stringstream linebuffer ( buffer );
			while ( linebuffer >> tmp )
			{
				if ( index < (unsigned int)StructCount )
				{
					if ( k < nd and h < nr[k] )
					{
						myScores[k][h][index] = tmp;
						if ( verbosity > 1 )
						{
							std::cout << "chain " << k << ", residue " << h << ", score n. " << index;
							std::cout << ": " << myScores[k][h][index] << std::endl;
						}
					}
					else
					{
						std::cout << "readScores[read]: domain indices out of bounds --" << std::flush;
						std::cout << " domain ID = " << k << "[" << nd << "]," << std::flush;
						std::cout << " residue ID = " << h << "[" << nr[k] << "]" << std::endl;
						exit ( 1 );
					}
				}
				index++;
			}
			h++;
		}
	}
	double tmp = 0.;
	unsigned int index = 0;
	stringstream linebuffer ( buffer );
	while ( linebuffer >> tmp )
	{
		if ( index < (unsigned int)StructCount )
		{
			if ( h < nr[nd - 1] )
			{
				myScores[nd - 1][h][index] = tmp;
				if ( verbosity > 1 )
				{
					std::cout << "chain " << nd - 1 << ", residue " << h << ", score n. " << index;
					std::cout << myScores[nd - 1][h][index] << std::endl;
				}
			}
			else
			{
				std::cout << "readScores[last] domain indices out of bounds --" << std::flush;
				std::cout << " domain ID = " << nd - 1 << "[" << nd << "]," << std::flush;
				std::cout << " residue ID = " << h << "[" << nr[nd - 1] << "]" << std::endl;
				exit ( 1 );
			}
		}
		index++;
	}
	return myScores;
}


double getConfidence ( double *myScores )
{
	double goldScore = -11.;
	double silverScore = -11.;
	for ( unsigned int p = 0; p < (unsigned int)StructCount; p++ )
	{
		if ( myScores[p] > silverScore )
		{
			if ( myScores[p] > goldScore )
			{
				silverScore = goldScore;
				goldScore = myScores[p];
			}
			else
			{
				silverScore = myScores[p];
			}
		}
	}
	double rawConfidence = goldScore - silverScore;
	return rawConfidence;
}


unsigned int deduceGuess ( double *myScores )
{
	unsigned int pmax = (unsigned int)StructCount - 1;
	double maxScore = -11.;
	for ( unsigned int p = 0; p < (unsigned int)StructCount; p++ )
	{
		if ( myScores[p] > maxScore )
		{
			maxScore = myScores[p];
			pmax = p;
		}
	}
	return pmax;
}


int computeID ( int t, unsigned short int window, int sstruct )
{
	int featID = ( t + (int)( window/2 ) ) * StructCount + sstruct;
// 	int featID = ( t + (int)( window/2 ) ) * ( StructCount - 1 ) + sstruct;
	return featID;
}

