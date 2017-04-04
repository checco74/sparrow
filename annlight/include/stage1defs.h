//#include "DBInterface.h"
#include <list>
#include <sstream>
#include "assert.h"
//#include "Optimizer.h"
#include "time.h"
//#include "OptimizerException.h"
#include "ParserException.h"
#include "CSVLogger.h"
#include "math.h"

#ifndef _STAGE1DEFS_H_
#define _STAGE1DEFS_H_

typedef enum t_positive_class
{
  // level 0 accuracy
  POS_HELIX,
  POS_STRAND,
  POS_COIL,
  // level 1 accuracy
  POS_ALPHA_HELIX,
  POS_P_HELIX,
  POS_T_HELIX,
  POS_BETA_STRAND,
  POS_RANDOM_COIL,
  // level 2 accuracy
  POS_ALPHA_HELIX_SHORT,
  POS_ALPHA_HELIX_START,
  POS_ALPHA_HELIX_CENTER,
  POS_ALPHA_HELIX_END,
  POS_P_HELIX_SHORT,
  POS_P_HELIX_START,
  POS_P_HELIX_CENTER,
  POS_P_HELIX_END,
  POS_T_HELIX_SHORT,
  POS_T_HELIX_START,
  POS_T_HELIX_CENTER,
  POS_T_HELIX_END,
  POS_BETA_STRAND_SHORT,
  POS_BETA_STRAND_START,
  POS_BETA_STRAND_CENTER,
  POS_BETA_STRAND_END,
  POS_COIL_LONG,
  POS_COIL_SHORT,
  // level 3 accuracy (can be combined with levels 1 and 2)
  POS_ISOLATED_STRAND,
  POS_TURN,
  POS_BEND,
  // level 4 accuracy (can be combined with all level)
  POS_UNDEFINED
};

typedef enum t_negative_class
{
  // level 0 accuracy
  NEG_HELIX,
  NEG_STRAND,
  NEG_COIL,
  // level 1 accuracy
  NEG_ALPHA_HELIX,
  NEG_P_HELIX,
  NEG_T_HELIX,
  NEG_BETA_STRAND,
  NEG_RANDOM_COIL,
  // level 2 accuracy
  NEG_ALPHA_HELIX_SHORT,
  NEG_ALPHA_HELIX_START,
  NEG_ALPHA_HELIX_CENTER,
  NEG_ALPHA_HELIX_END,
  NEG_P_HELIX_SHORT,
  NEG_P_HELIX_START,
  NEG_P_HELIX_CENTER,
  NEG_P_HELIX_END,
  NEG_T_HELIX_SHORT,
  NEG_T_HELIX_START,
  NEG_T_HELIX_CENTER,
  NEG_T_HELIX_END,
  NEG_BETA_STRAND_SHORT,
  NEG_BETA_STRAND_START,
  NEG_BETA_STRAND_CENTER,
  NEG_BETA_STRAND_END,
  NEG_COIL_LONG,
  NEG_COIL_SHORT,
  // level 3 accuracy (can be combined with levels 1 and 2)
  NEG_ISOLATED_STRAND,
  NEG_TURN,
  NEG_BEND,
  // level 4 accuracy (can be combined with all levels)
  NEG_ALL,
  NEG_UNDEFINED
};

DiscriminantParameters* readParameters(unsigned dim, bool linear, bool quadratic, const std::string& filename);
double getRMSD(bool linear, bool quadratic, DiscriminantParameters* param1, DiscriminantParameters* param2);
unsigned determineLevel(t_positive_class posClass);
unsigned determineLevel(t_negative_class negClass);
bool operator==(t_positive_class pClass, t_negative_class nClass);
bool combinationAllowed(t_positive_class pClass, t_negative_class nClass);
t_positive_class toPositiveClass(std::string str);
t_negative_class toNegativeClass(std::string str);
bool isHelix(std::string& str, unsigned keyPosition);
bool isStrand(std::string& str, unsigned keyPosition);
bool isCoil(std::string& str, unsigned keyPosition);
bool isAlphaHelix(std::string& str, unsigned keyPosition);
bool isPHelix(std::string& str, unsigned keyPosition);
bool isTHelix(std::string& str, unsigned keyPosition);
bool isBetaStrand(std::string& str, unsigned keyPosition);
bool isRandomCoil(std::string& str, unsigned keyPosition);
bool isAlphaHelixShort(std::string& str, unsigned keyPosition);
bool isAlphaHelixStart(std::string& str, unsigned keyPosition);
bool isAlphaHelixCenter(std::string& str, unsigned keyPosition);
bool isAlphaHelixEnd(std::string& str, unsigned keyPosition);
bool isPHelixShort(std::string& str, unsigned keyPosition);
bool isPHelixStart(std::string& str, unsigned keyPosition);
bool isPHelixCenter(std::string& str, unsigned keyPosition);
bool isPHelixEnd(std::string& str, unsigned keyPosition);
bool isTHelixShort(std::string& str, unsigned keyPosition);
bool isTHelixStart(std::string& str, unsigned keyPosition);
bool isTHelixCenter(std::string& str, unsigned keyPosition);
bool isTHelixEnd(std::string& str, unsigned keyPosition);
bool isBetaStrandShort(std::string& str, unsigned keyPosition);
bool isBetaStrandStart(std::string& str, unsigned keyPosition);
bool isBetaStrandCenter(std::string& str, unsigned keyPosition);
bool isBetaStrandEnd(std::string& str, unsigned keyPosition);
bool isCoilLong(std::string& str, unsigned keyPosition);
bool isCoilShort(std::string& str, unsigned keyPosition);
bool isIsolatedStrand(std::string& str, unsigned keyPosition);
bool isBend(std::string& str, unsigned keyPosition);
bool isTurn(std::string& str, unsigned keyPosition);
StructureType transformToStructureType(std::string& str);
void readTranscoderData(Preprocessor* preproc, std::string filename);
bool isPositive(t_positive_class positiveClass, t_negative_class negativeClass, std::string& str, unsigned keyPosition);
bool isNegative(t_positive_class positiveClass, t_negative_class negativeClass, std::string& str, unsigned keyPosition);
void transcodeDataset(Sample* data, unsigned size, Preprocessor* preproc);
std::pair<Sample*, Sample*> getSplitting(Sample* allData, unsigned allNumber, unsigned* posNumber, unsigned* negNumber, unsigned keyPosition, t_positive_class positiveClass, t_negative_class negativeClass);
Sample* mergeSets(Sample* first, unsigned firstSize, Sample* second, unsigned secondSize, bool deleteOriginal = true);
std::list<Sample>* splitToSubsets(Sample* wholeSet, unsigned length, unsigned subsets, bool deleteOriginal = true);
Sample* mergeToLearningSet(std::list<Sample>* sets, unsigned subsets, unsigned leaveOut, unsigned* learnSize);
Sample* getValidationSet(std::list<Sample>* sets, unsigned subsets, unsigned setToTake, unsigned* validationSize);
Protein* splitProtein(Protein* fullProt, unsigned* subSetSize);
Sample* proteinsToSamples(Protein* proteinSet, unsigned proteinCount, unsigned window, unsigned keyPosition, unsigned* sampleCount);
Sample* loadData(DBInterface* db, const std::string& filename, unsigned window, unsigned keyPosition, unsigned* setSize);

#endif
