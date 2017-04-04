#include "annDefs.h"
//#include "stage1defs.h"
#include <iostream>
#include <list>
#include <sstream>
//#include "ProfileReader.h"

typedef std::list<std::string> ResidueList;
typedef std::list<ResidueList*> ProteinList;

const std::string DATA_TEST_LOG = "datatest.csv";
const std::string READ_TEST_LOG = "readtest.csv";

DatasetType* readDataset(const std::string& filename, unsigned inputs, unsigned classNum, bool testMode, t_assignment_mode assignment)
{
  std::ifstream inFile(filename.c_str(), std::fstream::in);
  if (!inFile.is_open())
  {
    std::cerr << "ERROR: Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }

  DatasetType* dataset = new DatasetType;
  ProteinList pList;
  ResidueList* rList = new ResidueList;
  std::string line;
  while(getline(inFile, line, '\n'))
  {
    if (line.size() != 0)
    {
      if (line[0] == STAGE2_CHAIN_SEPARATOR)
      {
        if (rList->size() != 0)
        {
          pList.push_back(rList);
          rList = new ResidueList;
        }
      }
      else
      {
        rList->push_back(line);
      }
    }
  }
  if (rList->size() != 0)
  {
    pList.push_back(rList);
  }
  inFile.close();

  if (pList.size() == 0)
  {
    std::cerr << "ERROR: No protein data read from file \"" << filename << "\"!" << std::endl;
    exit(1);
  }
  dataset->data = new ProtData[pList.size()];
  dataset->size = pList.size();

  unsigned protID = 0;
  unsigned* counts = new unsigned[classNum];
  for (unsigned i = 0; i < classNum; i++) counts[i] = 0;
  for (ProteinList::iterator it = pList.begin(); it != pList.end(); it++)
  {
    ResidueList* prot = *it;
    if (prot->size() == 0)
    {
      std::cerr << "ERROR: Read protein with no residue data from file \"" << filename << "\"!" << std::endl;
      exit(1);
    }

    dataset->data[protID].inputs = new double*[prot->size()];
    if (testMode)
      dataset->data[protID].outputs = new double*[prot->size()];
    else
      dataset->data[protID].outputs = NULL;
    dataset->data[protID].evaluate = new bool[prot->size()];
    dataset->data[protID].expectedClass = new char[prot->size()];
    dataset->data[protID].size = prot->size();
    unsigned resID = 0;
    for (ResidueList::iterator it2 = prot->begin(); it2 != prot->end(); it2++)
    {
      bool evaluation = true;
      char expectedChar = '\0';
      std::pair<double*, double*>* res = parseLine(*it2, inputs, testMode, assignment, &evaluation, classNum, &expectedChar, counts);
      dataset->data[protID].inputs[resID] = res->first;
      if (testMode) dataset->data[protID].outputs[resID] = res->second;
      dataset->data[protID].evaluate[resID] = evaluation;
      dataset->data[protID].expectedClass[resID] = expectedChar;
      resID++;
      delete res;
    }
    protID++;
  }

  for (unsigned i = 0; i < classNum; i++) std::cout << "classNum[" << i << "] = " << counts[i] << std::endl;
  delete counts;
  for (ProteinList::iterator it = pList.begin(); it != pList.end(); it++) delete *it;

  return dataset;
}

DatasetType* readDataset2(const std::string& filename, unsigned inputs)
{
  std::ifstream inFile(filename.c_str(), std::fstream::in);
  if (!inFile.is_open())
  {
    std::cerr << "ERROR: Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }

  DatasetType* dataset = new DatasetType;
  ProteinList pList;
  ResidueList* rList = new ResidueList;
  std::list<std::pair<int, std::string>* > aList;
  std::string line;
  while(getline(inFile, line, '\n'))
  {
    if (line.size() != 0)
    {
      if (line[0] == STAGE2_CHAIN_SEPARATOR)
      {
        std::pair<int, std::string>* a = parseSeparatorLine(line);
        aList.push_back(a);
        if (rList->size() != 0)
        {
          pList.push_back(rList);
          rList = new ResidueList;
        }
      }
      else
      {
        rList->push_back(line);
      }
    }
  }
  if (rList->size() != 0)
  {
    pList.push_back(rList);
  }
  inFile.close();

  if (pList.size() == 0)
  {
    std::cerr << "ERROR: No protein data read from file \"" << filename << "\"!" << std::endl;
    exit(1);
  }
  if (pList.size() != aList.size())
  {
    std::cerr << "ERROR: Number of chains and number of domain IDs mismatch!" << std::endl;
    exit(1);
  }
  dataset->data = new ProtData[pList.size()];
  dataset->size = pList.size();

  unsigned protID = 0;
  std::list<std::pair<int, std::string>* >::iterator ait = aList.begin();
  for (ProteinList::iterator it = pList.begin(); it != pList.end(); it++)
  {
    ResidueList* prot = *it;
    if (prot->size() == 0)
    {
      std::cerr << "ERROR: Read protein with no residue data from file \"" << filename << "\"!" << std::endl;
      exit(1);
    }

    dataset->data[protID].inputs = new double*[prot->size()];
    dataset->data[protID].outputs = NULL;
    dataset->data[protID].evaluate = NULL;
    dataset->data[protID].size = prot->size();
    dataset->data[protID].domainID = (unsigned)(*ait)->first;
    dataset->data[protID].sequence = (*ait)->second;
    if (dataset->data[protID].sequence.size() != dataset->data[protID].size)
    {
      std::cerr << "ERROR: Sequence and domain length do not match!" << std::endl;
      exit(1);
    }
    unsigned resID = 0;
    for (ResidueList::iterator it2 = prot->begin(); it2 != prot->end(); it2++)
    {
      double* inputVec = parseLine(*it2, inputs);
      dataset->data[protID].inputs[resID] = inputVec;
      resID++;
    }
    protID++;
    ait++;
  }

  for (ProteinList::iterator it = pList.begin(); it != pList.end(); it++) delete *it;
  for (std::list<std::pair<int, std::string>* >::iterator it = aList.begin(); it != aList.end(); it++) delete *it;

  return dataset;
}

/*DatasetType* readDatasetMTX(const std::string& filename)
{
  DatasetType* dataset = new DatasetType;
  ProfileReader* pReader = new ProfileReader();

  pReader->read(filename);

  dataset->data = new ProtData[1];
  dataset->size = 1;
  dataset->data[0].size = pReader->getProteinLength();
  dataset->data[0].inputs = pReader->getDecodedProfile(true);
  dataset->data[0].outputs = NULL;

  delete pReader;
  return dataset;
}*/

/*DatasetType* readDatasetCoor(const std::string& filename, const ConnectionParameters* params, bool testMode)
{
  DatasetType* dataset = new DatasetType;

  unsigned setSize = 0;
  DBInterface* db = new DBInterface(params, 1.0, 1.0, false);
  Protein* protSet = db->loadTestset(filename, &setSize);

  std::list<Protein> chainList;
  for (unsigned i = 0; i < setSize; i++)
  {
    unsigned chainNum = 0;
    Protein* subset = splitProtein(&protSet[i], &chainNum);
    for (unsigned j = 0; j < chainNum; j++) chainList.push_back(subset[j]);
    delete [] subset;
  }

  dataset->size = chainList.size();
  dataset->data = new ProtData[chainList.size()];
  unsigned i = 0;
  for (std::list<Protein>::iterator it = chainList.begin(); it != chainList.end(); it++)
  {
    unsigned protLength = it->sequence.size();
    dataset->data[i].size = protLength;
    dataset->data[i].inputs = Transformations::getProfile(it->profile, it->sequence, true);
    if (testMode)
    {
      dataset->data[i].outputs = new double*[protLength];
      for (unsigned j = 0; j < protLength; j++)
      {
        dataset->data[i].outputs[j] = new double[STAGE2_OUTPUT_NODES];
        Category trueCat = Transformations::struct2cat(Transformations::char2struct(it->sequence[j]));
        for (unsigned k = 0; k < STAGE2_OUTPUT_NODES; k++)
        {
          if (((unsigned)trueCat) == k)
            dataset->data[i].outputs[j][k] = ANN::LITERATURE_POSITIVE_EXPECTED;
          else
            dataset->data[i].outputs[j][k] = ANN::LITERATURE_NEGATIVE_EXPECTED;
        }
      }
    }
    else
    {
      dataset->data[i].outputs = NULL;
    }
    i++;
  }

  delete [] protSet;
  delete db;
  return dataset;
}*/

void freeData(DatasetType* dataset)
{
  for (unsigned prot = 0; prot < dataset->size; prot++)
  {
    for (unsigned res = 0; res < dataset->data[prot].size; res++)
    {
      delete [] dataset->data[prot].inputs[res];
      if (dataset->data[prot].outputs != NULL) delete [] dataset->data[prot].outputs[res];
    }
    delete [] dataset->data[prot].inputs;
    if (dataset->data[prot].outputs != NULL) delete [] dataset->data[prot].outputs;
    delete [] dataset->data[prot].evaluate;
  }
  delete [] dataset->data;
  delete dataset;
}

unsigned getClassAssignment(std::string str, t_assignment_mode assignment, bool* evaluation)
{
  assert(assignment != AM_NONE);
  assert(str.size() == 1);
  *evaluation = false;
  char c = str[0];
  switch (c)
  {
    case 'H':
      return 0;
    case 'E':
      return 1;
    case 'C':
    case 'T':
    case 'S':
    case 'B':
    case '?':
      return 2;
    case 'G':
    case 'I':
      if (assignment == AM_STRICT)
      {
        *evaluation = true;
        return 2;
      }
      else
      {
        *evaluation = true;
        return 0;
      }
    default:
      std::cerr << "ERROR: Unknown SSE class: \"" << c << "\"" << std::endl;
      exit(1);
  }
}

std::pair<double*, double*>* parseLine(const std::string& line, unsigned inputs, bool testMode, t_assignment_mode assignment, bool* evaluation, unsigned classNum, char* expectedChar, unsigned* classCount)
{
  double* inputValues = new double[inputs];
  double* outputValues = NULL;
  if (testMode) outputValues = new double[classNum];

  unsigned i = 0;

  // read inputs
  unsigned inputsRead = 0;
  while (inputsRead < inputs && i < line.size())
  {
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) i++;  // skip whitespaces in line
    std::ostringstream valueStr;
    while (i < line.size() && line[i] != ' ' && line[i] != '\t')
    {
      valueStr << line[i];
      i++;
    }
    inputValues[inputsRead] = atof(valueStr.str().c_str());
    inputsRead++;
  }
  if (inputsRead < inputs)
  {
    std::cerr << "ERROR: Missing input values in line \"" << line << "\" of input file!" << std::endl;
    exit(1);
  }

  // read outputs
  if (testMode)
  {
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) i++;
    std::ostringstream valueStr;
    while (i < line.size() && line[i] != ' ' && line[i] != '\t')
    {
      valueStr << line[i];
      i++;
    }
    unsigned expectedClass;
    if (assignment == AM_NONE)
      expectedClass = atoi(valueStr.str().c_str());
    else
      expectedClass = getClassAssignment(valueStr.str(), assignment, evaluation);
    if (expectedClass >= classNum)
    {
      std::cerr << "ERROR: Expected value in line \"" << line << "\" is wrong!" << std::endl;
      exit(1);
    }
    classCount[expectedClass]++;
    //if (expectedClass == 0) std::cerr << line << std::endl;
    for (unsigned k = 0; k < classNum; k++)
    {
      if (k == expectedClass)
        outputValues[k] = ANN::LITERATURE_POSITIVE_EXPECTED;
      else
        outputValues[k] = ANN::LITERATURE_NEGATIVE_EXPECTED;
    }
    
    // read optional additional class information
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) i++;
    if (i < line.size() && line[i] != '\0' && line[i] != '\n') *expectedChar = line[i];
  }

  std::pair<double*, double*>* retPair = new std::pair<double*, double*>(inputValues, outputValues);

  return retPair;
}

double* parseLine(const std::string& line, unsigned inputs)
{
  double* inputValues = new double[inputs];
  unsigned i = 0;

  // read inputs
  unsigned inputsRead = 0;
  while (inputsRead < inputs && i < line.size())
  {
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) i++;  // skip whitespaces in line
    std::ostringstream valueStr;
    while (i < line.size() && line[i] != ' ' && line[i] != '\t')
    {
      valueStr << line[i];
      i++;
    }
    inputValues[inputsRead] = atof(valueStr.str().c_str());
    inputsRead++;
  }
  if (inputsRead < inputs)
  {
    std::cerr << "ERROR: Missing input values in line \"" << line << "\" of input file!" << std::endl;
    exit(1);
  }

  return inputValues;
}

std::pair<int, std::string>* parseSeparatorLine(std::string& line)
{
  unsigned i = 0;
  //skip separator and keyword
  while (i < line.size() && (line[i] != ' ' && line[i] != '\t')) i++;
  //skip whitespaces
  while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) i++;
  std::ostringstream vStream;
  while (i < line.size() && (line[i] != ' ' && line[i] != '\t'))
  {
    vStream << line[i];
    i++;
  }
  int ID = atoi(vStream.str().c_str());
  //skip whitespaces
  while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) i++;
  std::ostringstream sStream;
  while (i < line.size() && (line[i] != ' ' && line[i] != '\t'))
  {
    sStream << line[i];
    i++;
  }
  std::string sequence = sStream.str();
  
  std::pair<int, std::string>* p = new std::pair<int, std::string>(ID, sequence);
  
  return p;
}


void readTest(const std::string& filename, unsigned inputs)
{
  DatasetType* dataset = readDataset(filename, inputs);

  std::ofstream outFile(READ_TEST_LOG.c_str(), std::fstream::out | std::fstream::trunc);
  if (!outFile.is_open())
  {
    std::cerr << "ERROR: Could not open file \"" << READ_TEST_LOG << "\" for writing!" << std::endl;
    exit(1);
  }
  for (unsigned i = 0; i < dataset->size; i++)
  {
    outFile << "# chain " << i << std::endl;
    for (unsigned j = 0; j < dataset->data[i].size; j++)
    {
      for (unsigned k = 0; k < inputs; k++)
      {
        outFile << dataset->data[i].inputs[j][k] << '\t';
      }
      for (unsigned k = 0; k < inputs; k++)
      {
        if (dataset->data[i].outputs[j][k] == ANN::LITERATURE_POSITIVE_EXPECTED) outFile << k << std::endl;
      }
    }
  }
  outFile.close();

  freeData(dataset);
}

double* getSubSequence(DatasetType* dataset, unsigned protID, unsigned resID, unsigned windowSize, unsigned keyPosition, unsigned inputs)
{
  unsigned size = inputs * windowSize;  //size of the returned input vector
  double* retValues = new double[size];

  // fill lower vector
  unsigned vecIndex = keyPosition - 1;
  unsigned dataIndex = resID - 1;
  unsigned lower = keyPosition;
  if (keyPosition == 0) lower = 0;
  bool inBreak = (resID == 0);
  while (lower != 0)
  {
    if (inBreak)
    {
      unsigned vecInd2 = inputs * vecIndex;
      for (unsigned i = 0; i < inputs; i++)
      {
        retValues[vecInd2] = TERMINUS_DEFAULT_VALUE;
        vecInd2++;
      }
      vecIndex--;
    }
    else
    {
      unsigned vecInd2 = inputs * vecIndex;
      for (unsigned i = 0; i < inputs; i++)
      {
        retValues[vecInd2] = dataset->data[protID].inputs[dataIndex][i];
        vecInd2++;
      }
      vecIndex--;
      inBreak = (dataIndex == 0);
      dataIndex--;
    }
    lower--;
  }

  // fill upper vector
  unsigned upper = windowSize - keyPosition - 1;
  vecIndex = keyPosition + 1;
  dataIndex = resID + 1;
  inBreak = (resID == dataset->data[protID].size - 1);
  while (upper != 0)
  {
    if (inBreak)
    {
      unsigned vecInd2 = inputs * vecIndex;
      for (unsigned i = 0; i < inputs; i++)
      {
        retValues[vecInd2] = TERMINUS_DEFAULT_VALUE;
        vecInd2++;
      }
      vecIndex++;
    }
    else
    {
      unsigned vecInd2 = inputs * vecIndex;
      for (unsigned i = 0; i < inputs; i++)
      {
        retValues[vecInd2] = dataset->data[protID].inputs[dataIndex][i];
        vecInd2++;
      }
      vecIndex++;
      inBreak = (dataIndex == dataset->data[protID].size - 1);
      dataIndex++;
    }
    upper--;
  }

  // fill key position
  vecIndex = inputs * keyPosition;
  for (unsigned i = 0; i < inputs; i++)
  {
    retValues[vecIndex] = dataset->data[protID].inputs[resID][i];
    vecIndex++;
  }

  return retValues;
}

double* getSubSequence(DatasetType* dataset, unsigned protID, unsigned resID, unsigned windowSize, unsigned inputs)
{
  if (windowSize % 2 == 0)
  {
    std::cerr << "ERROR: Window size has to be odd! Current window size: " << windowSize << std::endl;
    exit(1);
  }
  unsigned size = inputs * windowSize + 2; // size of the returned input vector
  unsigned lowerBound = windowSize / 2; // lower bound up to which terminal residues have to be taken into account
  unsigned upperBound = dataset->data[protID].size - lowerBound; // upper bound starting with which termainal residues have to ba taken into account

  double* retValues = new double[size];
  bool lowerTerminus = (resID < lowerBound);
  bool upperTerminus = (resID >= upperBound) || (dataset->data[protID].size < lowerBound);

  unsigned currentIndex = 0;
  unsigned lowerFill = 0;
  unsigned effectiveWindow = windowSize;
  if (lowerTerminus)
  {
    lowerFill = lowerBound - resID;
    for (unsigned i = 0; i < lowerFill; i++)
    {
      for (unsigned v = 0; v < inputs; v++)
      {
        retValues[currentIndex] = TERMINUS_DEFAULT_VALUE;
        currentIndex++;
      }
      effectiveWindow--;
    }
  }

  unsigned startRes = 0;
  if (!lowerTerminus) startRes = resID - lowerBound;
  unsigned stopRes = dataset->data[protID].size;
  if (!upperTerminus) stopRes = resID + lowerBound + 1;

  assert (stopRes - startRes <= effectiveWindow);

  for (unsigned i = startRes; i < stopRes; i++)
  {
    for (unsigned j = 0; j < inputs; j++)
    {
      retValues[currentIndex] = dataset->data[protID].inputs[i][j];
      currentIndex++;
    }
    effectiveWindow--;
  }

  assert(upperTerminus && (effectiveWindow > 0) || !upperTerminus && (effectiveWindow == 0));
  while (effectiveWindow > 0)
  {
    for (unsigned v = 0; v < inputs; v++)
    {
      retValues[currentIndex] = TERMINUS_DEFAULT_VALUE;
      currentIndex++;
    }
    effectiveWindow--;
  }

  if (lowerTerminus)
    retValues[size - 2] = ANN::LITERATURE_POSITIVE_EXPECTED;
  else
    retValues[size - 2] = ANN::LITERATURE_NEGATIVE_EXPECTED;

  if (upperTerminus)
    retValues[size - 1] = ANN::LITERATURE_POSITIVE_EXPECTED;
  else
    retValues[size - 1] = ANN::LITERATURE_NEGATIVE_EXPECTED;

  return retValues;
}

void dataTest(const std::string& filename, unsigned windowSize, unsigned inputs)
{
  unsigned vecLength = windowSize * inputs + 2;

  DatasetType* data = readDataset(filename, inputs);

  std::ofstream outFile(DATA_TEST_LOG.c_str(), std::fstream::out | std::fstream::trunc);
  if (!outFile.is_open())
  {
    std::cerr << "ERROR: Could not open file \"" << DATA_TEST_LOG << "\" for writing!" << std::endl;
    exit(1);
  }

  for (unsigned p = 0; p < data->size; p++)
  {
    outFile << "#chain " << p << std::endl;
    for (unsigned i = 0; i < data->data[p].size; i++)
    {
      double* values = getSubSequence(data, p, i, windowSize, inputs);

      for (unsigned k = 0; k < vecLength; k++)
      {
        if (k != 0) outFile << '\t';
        outFile << values[k];
      }
      outFile << std::endl;

      delete [] values;
    }
  }

  freeData(data);

  outFile.close();
}

ActivationType toActivationType(const std::string& str)
{
  std::string changedStr = str;
  for (unsigned i = 0; i < changedStr.size(); i++)
  {
    changedStr[i] = std::tolower(changedStr[i]);
  }

  if (changedStr == "tanh")
    return ANN_TANH;
  else if (changedStr == "logistic")
    return ANN_LOGISTIC;
  else
  {
    std::cerr << "ERROR: Unknown activation function type: " << str << std::endl;
    exit(1);
  }
}

RegularizationType toRegularizationType(const std::string& str)
{
  std::string changedStr = str;
  for (unsigned i = 0; i < changedStr.size(); i++)
  {
    changedStr[i] = std::tolower(changedStr[i]);
  }

  if (changedStr == "decay")
    return W_DECAY;
  else if (changedStr == "elimination")
    return W_ELIMINATION;
  else if (changedStr == "smoother")
    return W_APPROXIMATE_SMOOTHER;
  else
  {
    std::cerr << "ERROR: Unknown activation function type: " << str << std::endl;
    exit(1);
  }
}

unsigned getSize(const DatasetType* dataset)
{
  unsigned retValue = 0;

  for (unsigned i = 0; i < dataset->size; i++)
  {
    retValue += dataset->data[i].size;
  }

  return retValue;
}

/*std::pair<double**, double**> getSplitting(double** allData, std::string* allCats, unsigned allNumber, unsigned* posNumber, unsigned* negNumber, t_positive_class positiveClass, t_negative_class negativeClass)
{
  if (positiveClass == negativeClass)
  {
    std::cerr << std::endl << "ERROR: Positive and negative class must not be the same!" << std::endl;
    exit(1);
  }
  if (positiveClass == POS_UNDEFINED)
  {
    std::cerr << std::endl << "ERROR: Positive class not defined!" << std::endl;
    exit(1);
  }
  if (negativeClass == NEG_UNDEFINED)
  {
    std::cerr << std::endl << "ERROR: Negative class not defined!" << std::endl;
    exit(1);
  }
  if (!combinationAllowed(positiveClass, negativeClass))
  {
    std::cerr << std::endl << "ERROR: The given combination of positive an negative class is not allowed! Classes have to be disjoint!" << std::endl;
    exit(1);
  }
  unsigned posCount = 0;
  unsigned negCount = 0;
  for (unsigned i = 0; i < allNumber; i++)
  {
    if (isPositive(positiveClass, negativeClass, allCats[i], 0)) posCount++;
    else if (isNegative(positiveClass, negativeClass, allCats[i], 0)) negCount++;
  }

  *posNumber = posCount;
  *negNumber = negCount;
  double** posSet = new double*[posCount];
  double** negSet = new double*[negCount];
  unsigned posIndex = 0;
  unsigned negIndex = 0;
  for (unsigned i = 0; i < allNumber; i++)
  {
    if (isPositive(positiveClass, negativeClass, allCats[i], 0))
      posSet[posIndex++] = allData[i];
    else if (isNegative(positiveClass, negativeClass, allCats[i], 0))
      negSet[negIndex++] = allData[i];
  }

  std::pair<double**, double**> myPair = std::make_pair(posSet, negSet);
  return myPair;
}*/

t_assignment_mode toAssignmentMode(std::string str)
{
  if (str == "strict")
    return AM_STRICT;
  else if (str == "loose")
    return AM_LOOSE;
  else if (str == "none")
    return AM_NONE;
  else
  {
    std::cerr << "ERROR: Unknown assignment mode: \"" << str << "\"" << std::endl;
    exit(1);
  }
}

