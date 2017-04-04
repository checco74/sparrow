/* ScoreGear.h - version 1.22prof
	- note to version 1.22prof:	* scoreMachine and scoreGear have been merged. */

#ifndef _SCORE_GEAR
#define _SCORE_GEAR

#include "define.h"
#include "aatype.h"

using namespace std;

class scoreMachine
{
	public:
		/** @brief secondary structure pattern */
		int Pattern;
		/** @brief eventual antagonist secondary structure pattern */
		int Antipattern;
		/** @brief tells whether array space has been allocated */
		bool isSet;
		/** @brief tells whether the object has been initialized */
		bool isInit;
		/** @brief tells whether the corresponding scorer is available */
		bool isAvailable;
		/** @brief switches the quadratic-form terms on/off */
		bool power;
		/** @brief switches the linear-form terms on/off */
		bool linear;
		/** @brief switches the memory saving mode (i.e. categories) for linear terms on/off; */
		bool linsave;
		/** @brief switches the memory saving mode (i.e. categories) for linear terms on/off; */
		bool collapse;
		/** @brief sequence window length */
		uint winlen;
		/** @brief sequence window center */
		uint winkey;
		/** @brief number of linear terms of the scoring function */
		uint Ntot;
		/** @brief rank of the bilinear form of the scoring function */
		uint Mtot;
		/** @brief number of parameters of the scoring function */
		uint dim;
		/** @brief the pattern reweighting factor */
		double *reweight;
		/** @brief relative frequencies for each protein zone */
		double **Xdx, **Xsx;
		/** @brief solution parameters of the linear systems for each protein CATH class and zone */
		double **Wdx, **Wsx;
		/** @brief constant terms for each protein CATH class and zone */
		double *bdx, *bsx;

		scoreMachine();

		~scoreMachine();

		scoreMachine &operator=(const scoreMachine &right)
		{
			this->power = right.power;
			this->linear = right.linear;
			this->linsave = right.linsave;
			this->winlen = right.winlen;
			this->winkey = right.winkey;
			return *this;
		}

		/** initializes the scoring machine
		@param first secondary structure pattern;
		@param second eventual antagonist secondary structure pattern. */
		void init(int first = -1, int second = -1);

		/** allocates arrays space */
		void setArraysSpace();

		/** releases arrays space */
		void releaseArraysSpace();

		/** retrieves the scoring function parameters from disk
		@param mode tells which kind of parameters the function shall look for */
		void getComponents(readMode mode = STD_MODE);

		/** computes the scoring function's constant term 'b'
		@param pattern secondary structure pattern: positive samples statistic;
		@param zone protein region code;
		@param dir direction along the protein chain.
		@return constant term 'b' */
		double computeConstTerm(int pattern, int zone, int dir);

		/** computes the score
		@param zone protein region code;
		@param prof sequence window's profile;
		@param profLength length of the profile;
		@param dir direction along the protein chain.
		@return score of the amino acid sequence with the given profile */
		double getScore(int zone, double *prof, int profLength, int dir);

		/** says whether the scoring function is feasible. */
		bool isFeasible();
};

/** @brief various scoring function tools */
class scoreGear
{
	public:
		/** @brief secondary structure pattern */
		int Pattern;
		/** @brief auxiliary scorers switch */
		bool aux;
		/** @brief scoring functions (dimensions: type, id) */
		scoreMachine **sm, **auxn_sm, **auxa_sm;
		/** @brief number of scorers per pattern and secondary structure pair */
		int *partNS, *partauxnNS, *partauxaNS;
		/** @brief total number of scorers (i.e. free "super"-parameters) per pattern */
		int totalNS, totalauxnNS, totalauxaNS;
		/** @brief number of super scoring function coefficients */
		int NSSC, auxnNSSC, auxaNSSC;
		/** @brief super scoring function coefficients */
		pair <double**, double**> ssc, auxn_ssc, auxa_ssc;

		scoreGear();

		~scoreGear();

		/** sets up and retrieves the whole super-scoring functions machinery
		@param aux decides whether auxiliary scorers shall be retrieved or not. */
		void setup(bool aux = false);

		/** allocates infrastructural memory for scoreMachines */
		void setScoreTools();

		/** sets up and retrieves the scoring functions' parameters
		@param rmode read mode: identifies the scorers to be retrieved;
		@param NS partial number of scorers for the current structure pair;
		@param p0 first secondary structure pattern;
		@param p1 eventual second secondary structure pattern.
		@return a pointer to the array of <code>scoreMachines</code>. */
		scoreMachine *getScorers(readMode rmode, int NS, int p0 = -1, int p1 = -1);

		/** sets up the super-parameters machinery
		@param rmode read mode: identifies the scorers to be retrieved;
		@param N number of super parameters. */
		void getSuperCoefficients(readMode rmode, int N);

		/** retrieves the super-parameters machinery
		@param z protein region code;
		@param tag file identification tag;
		@param N number of super scoring function coefficients;
		@param rssc pair of super score coefficients for the two directions along the chain;
		@param rmode read mode: identifies the scorers to be retrieved. */
		void getSuperCoefficients(int z, string &tag, int N, pair <double**, double**> &rssc, readMode rmode);

		/** releases infrastructural scoreMachines memory */
		void eraseScoreTools();

		/** releases super-scoring function memory */
		void eraseSuperCoefficients();
};

#endif
