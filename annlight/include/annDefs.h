#include <string>
#include "ANN.h"
//#include "DBInterface.h"
//#include "stage1defs.h"

#ifndef _STAGE2DEFS_H_
#define _STAGE2DEFS_H_

//typedef enum t_positive_class;
//typedef enum t_negative_class;

/**
 * @brief Class to hold the data of a protein as needed be the <code>Optimizer</code>.
 */
typedef struct ProtData
{
  /** @brief The input vectors for each residue in the protein. */
  double** inputs;
  /** @brief The output vectors for each residue in the protein. */
  double** outputs;
  /** @brief The length of the protein. */
  unsigned size;
  /** @brief Determines whether a certain residue should be used in the evaluation. */
  bool* evaluate;
  /** @brief The DSSP character of the expected class used for evaluation of the program. */
  char* expectedClass;
  /** @brief The internal ID of the domain this chain belongs to. */
  unsigned domainID;
  /** @brief The sequence for this domain. */
  std::string sequence;
};

/**
 * @brief Wrapper compound to hold the <code>ProtData</code> of a dataset.
 */
typedef struct DatasetType
{
  /** @brief The data of all proteins in the dataset. */
  ProtData* data;
  /** @brief The number of proteins in the dataset. */
  unsigned size;
};

typedef enum t_assignment_mode
{
  AM_NONE,
  AM_STRICT,
  AM_LOOSE
};

const unsigned STAGE2_OUTPUT_NODES = 3;
const char STAGE2_CHAIN_SEPARATOR = '#';
const double TERMINUS_DEFAULT_VALUE = -1.0;

DatasetType* readDataset(const std::string& filename, unsigned inputs, unsigned classNum = 3, bool testMode = true, t_assignment_mode assignment = AM_NONE);
DatasetType* readDataset2(const std::string& filename, unsigned inputs);
//DatasetType* readDatasetMTX(const std::string& filename);
//DatasetType* readDatasetCoor(const std::string& filename, const ConnectionParameters* params, bool testMode = true);
std::pair<double*, double*>* parseLine(const std::string& line, unsigned inputs, bool testMode, t_assignment_mode assignment, bool* evaluation, unsigned classNum, char* expectedChar, unsigned* classCount = NULL);
double* parseLine(const std::string& line, unsigned inputs);
unsigned getSize(const DatasetType* dataset);
void freeData(DatasetType* dataset);

std::pair<int, std::string>* parseSeparatorLine(std::string& line);

double* getSubSequence(DatasetType* dataset, unsigned protID, unsigned resID, unsigned windowSize, unsigned inputs);
double* getSubSequence(DatasetType* dataset, unsigned protID, unsigned resID, unsigned windowSize, unsigned keyPosition, unsigned inputs);

//std::pair<double**, double**> getSplitting(double** allData, std::string* allCats, unsigned allNumber, unsigned* posNumber, unsigned* negNumber, t_positive_class positiveClass, t_negative_class negativeClass);

void readTest(const std::string& filename, unsigned inputs);
void dataTest(const std::string& filename, unsigned windowSize, unsigned inputs);
ActivationType toActivationType(const std::string& str);
RegularizationType toRegularizationType(const std::string& str);

unsigned getClassAssignment(std::string str, t_assignment_mode assignment, bool* evaluation);
t_assignment_mode toAssignmentMode(std::string str);

#endif
