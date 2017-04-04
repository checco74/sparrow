#include "ANN.h"
#include "annDefs.h"
#include <iostream>
#include <sstream>
#include "ParameterParser.h"
#include "DefaultParameters.h"
#include "ANNException.h"
#include "ParserException.h"
#include "KMeansClusterer.h"
#include "math.h"

const std::string DATA_LEARN_FILE_NAME = "data.learnfile";
const std::string DATA_LEARN_FILE_DEFAULT = "";
const std::string DATA_TEST_FILE_NAME = "data.testfile";
const std::string DATA_TEST_FILE_DEFAULT = "";
const std::string ANN_MIN_ERROR_NAME = "ann.minerror";
const std::string ANN_MIN_ERROR_DEFAULT = "1e-6";
const std::string ANN_MAX_ITERATIONS_NAME = "ann.maxIt";
const std::string ANN_MAX_ITERATIONS_DEFAULT = "1000";
const std::string ANN_SHUFFLE_NAME = "ann.shuffle";
const std::string ANN_SHUFFLE_DEFAULT = "1";
const std::string ANN_SHUFFLE_SEED_NAME = "ann.shuffleSeed";
const std::string ANN_SHUFFLE_SEED_DEFAULT = "0";
const std::string ANN_RPROP_NAME = "ann.rprop";
const std::string ANN_RPROP_DEFAULT = "";
const std::string ANN_NORMALIZATION_NAME = "ann.normalization";
const std::string ANN_NORMALIZATION_DEFAULT = "";
const std::string ANN_LITERATURE_NAME = "ann.useLiteratureValues";
const std::string ANN_LITERATURE_DEFAULT = "1";
const std::string ANN_ACTIVATION_TYPE_NAME = "ann.activationType";
const std::string ANN_ACTIVATION_TYPE_DEFAULT = "tanh";
const std::string ANN_A_VALUE_NAME = "ann.aValue";
const std::string ANN_A_VALUE_DEFAULT = "1";
const std::string ANN_B_VALUE_NAME = "ann.bValue";
const std::string ANN_B_VALUE_DEFAULT = "1";
const std::string ANN_LEARNING_CONSTANT_NAME = "ann.learningConstant";
const std::string ANN_LEARNING_CONSTANT_DEFAULT = "0.1";
const std::string ANN_MOMENTUM_CONSTANT_NAME = "ann.momentumConstant";
const std::string ANN_MOMENTUM_CONSTANT_DEFAULT = "0.6";
const std::string ANN_VALIDATION_MAXIT_NAME = "ann.validationMaxIt";
const std::string ANN_VALIDATION_MAXIT_DEFAULT = "10";
const std::string ANN_ADJUST_FREQUENCY_NAME = "ann.adjustFrequency";
const std::string ANN_ADJUST_FREQUENCY_DEFAULT = "100";
const std::string ANN_ADJUST_CONSTANTS_NAME = "ann.adjustConstants";
const std::string ANN_ADJUST_CONSTANTS_DEFAULT = "0";
const std::string ANN_ADJUST_FACTOR_NAME = "ann.adjustFactor";
const std::string ANN_ADJUST_FACTOR_DEFAULT = "0.1";
const std::string ANN_ADJUST_MIN_NAME = "ann.adjustMin";
const std::string ANN_ADJUST_MIN_DEFAULT = "1e-6";
const std::string ANN_ADJUST_MAX_NAME = "ann.adjustMax";
const std::string ANN_ADJUST_MAX_DEFAULT = "50";
const std::string ANN_RHO_MINUS_NAME = "ann.rhoMinus";
const std::string ANN_RHO_MINUS_DEFAULT = "0.5";
const std::string ANN_RHO_PLUS_NAME = "ann.rhoPlus";
const std::string ANN_RHO_PLUS_DEFAULT = "1.2";
const std::string ANN_SEED_NAME = "ann.seed";
const std::string ANN_SEED_DEFAULT = "2010";
const std::string ANN_RANDOM_BOUND_NAME = "ann.rndBound";
const std::string ANN_RANDOM_BOUND_DEFAULT = "0.5";
const std::string ANN_LAMBDA_NAME = "ann.lambda";
const std::string ANN_LAMBDA_DEFAULT = "1e-10";
const std::string ANN_W_ZERO_NAME = "ann.w_zero";
const std::string ANN_W_ZERO_DEFAULT = "1e-3";
const std::string ANN_REGULARIZATION_NAME = "ann.regularization";
const std::string ANN_REGULARIZATION_DEFAULT = "decay";
const std::string ANN_POWER_NAME = "ann.power";
const std::string ANN_POWER_DEFAULT = "2";
const std::string ANN_HIDDEN_NAME = "ann.hidden";
const std::string ANN_HIDDEN_DEFAULT = "";
const std::string ANN_OUTPUT_NAME = "ann.output";
const std::string ANN_OUTPUT_DEFAULT = "3";
const std::string ANN_SCORES_PER_RES_NAME = "ann.scores";
const std::string ANN_SCORES_PER_RES_DEFAULT = "3";
const std::string ANN_WINDOW_SIZE_NAME = "ann.window";
const std::string ANN_WINDOW_SIZE_DEFAULT = "";
const std::string ANN_PERIOD_NAME = "ann.period";
const std::string ANN_PERIOD_DEFAULT = "";
const std::string LOAD_PARAMETERS_NAME = "loadParameters";
const std::string LOAD_PARAMETERS_DEFAULT = "0";
const std::string LOAD_FILE_NAME = "loadFile";
const std::string LOAD_FILE_DEFAULT = "";
const std::string TESTMODE_NAME = "testmode";
const std::string TESTMODE_DEFAULT = "0";
const std::string ASSIGNMENT_MODE_TEST_NAME = "assignmentModeTest";
const std::string ASSIGNMENT_MODE_TEST_DEFAULT = "none";
const std::string CLUSTER_USE_NAME = "clusters.use";
const std::string CLUSTER_USE_DEFAULT = "0";
const std::string CLUSTER_FILE_NAME = "clusters.filename";
const std::string CLUSTER_FILE_DEFAULT = "prototypes.prt";
const std::string EXTENDED_NAME = "useExtendedModel";
const std::string EXTENDED_DEFAULT = "";
const std::string SCORES_WRITE_NAME = "scores.write";
const std::string SCORES_WRITE_DEFAULT = "0";
const std::string SCORES_FILE_NAME = "scores.file";
const std::string SCORES_FILE_DEFAULT = "scores.csv";
const std::string SCORES_FILE_PRED_NAME = "scores.filePred";
const std::string SCORES_FILE_PRED_DEFAULT = "scores.pred.csv";

const std::string CONFIG_FILE_NAME = "stage2learn.cfg";
const std::string LOG_FILE_NAME = "stage2learn.log";

const std::string PARAMETER_PATH = "parameters";
const std::string PATH_SEPARATOR = "/";
const std::string ANN_PREFIX = "ann";
const std::string FILE_SUFFIX = ".prm";
const std::string LOG_NAME = "stage2learn.csv";

const std::string WRONG_RECALL_LOG_NAME = "recall_wrong.log";
const std::string WRONG_PREDICTION_LOG_NAME = "prediction_wrong.log";


const double TINY = 1e-20;

std::string getClassName(unsigned index, bool extended)
{
  switch (index)
  {
    case 0:
      return "alpha helix";
    case 1:
      return "beta strand";
    case 2:
      return "random coil";
    case 3:
    {
      if (extended)
        return "turn";
      else
        return "3-helix";
    }
    case 4:
    {
      if (extended)
        return "3-helix";
      else
        return "5-helix";
    }
    case 5:
      return "turn";
    case 6:
      return "bend";
    case 7:
      return "isolated beta bridge";
    default:
      std::cerr << "ERROR: Invalid class index: " << index << std::endl;
      exit(1);
  }
  return "";
}

char getClassLetter(unsigned index, bool extended)
{
  switch (index)
  {
    case 0:
      return 'H';
    case 1:
      return 'E';
    case 2:
      return 'C';
    case 3:
    {
      if (extended)
        return 'T';
      else
        return 'G';
    }
    case 4:
    {
      if (extended)
        return 'G';
      else
        return 'I';
    }
    case 5:
      return 'T';
    case 6:
      return 'S';
    case 7:
      return 'B';
    default:
      std::cerr << "ERROR: Invalid class index: " << index << std::endl;
      exit(1);
  }
  return ' ';
}

double computeCorrelationCoefficient(double **(&confusionMatrix), int catalogSize)
{
  for (int i = 0; i < catalogSize; i++)
  {
    for (int j = 0; j < catalogSize; j++)
    {
      std::cout << "confusionMatrix[" << i << "][" << j << "] = " << confusionMatrix[i][j] << std::endl;
    }
  }
  double numerator = 0.;
  double denominator = 0.;
  double sum_over_k_A = 0.;
  double sum_over_k_B = 0.;
  for(int k = 0; k < catalogSize; k++)
  {
    double sum_over_i_A = 0.;
    double sum_over_i_B = 0.;
    double sum_over_ij_A = 0.;
    double sum_over_ij_B = 0.;
    for(int i = 0; i < catalogSize; i++)
    {
      for(int j = 0; j < catalogSize; j++)
      {
        numerator += confusionMatrix[k][k] * confusionMatrix[i][j] - confusionMatrix[k][i] * confusionMatrix[j][k];
        if(j!=k)
        {
          sum_over_ij_A += confusionMatrix[j][i];
          sum_over_ij_B += confusionMatrix[i][j];
        }
      }
      sum_over_i_A += confusionMatrix[k][i];
      sum_over_i_B += confusionMatrix[i][k];
    }
    sum_over_k_A += sum_over_i_A * sum_over_ij_A;
    sum_over_k_B += sum_over_i_B * sum_over_ij_B;
  }
  denominator = sqrt(sum_over_k_A * sum_over_k_B + TINY) ;
  std::cout << "numerator = " << numerator << std::endl;
  std::cout << "denominator = " << denominator << std::endl;
  return numerator/denominator ;
}

int main(unsigned argc, char** argv)
{
  try
  {
    std::ofstream log(LOG_FILE_NAME.c_str(), std::fstream::out | std::fstream::trunc);
    if (!log.is_open())
    {
      std::cerr << "Could not open log file for writing!" << std::endl;
      exit(1);
    }
    ParameterParser* parser = new ParameterParser();
    DefaultParameters* dParam = new DefaultParameters();
    const std::list<std::pair<std::string, std::string> > defaults = dParam->getParameters();
    for (std::list<std::pair<std::string, std::string> >::const_iterator it = defaults.begin(); it != defaults.end(); it++)
    {
      parser->addParameter(it->first, it->second);
    }
    parser->addParameter(DATA_LEARN_FILE_NAME, DATA_LEARN_FILE_DEFAULT);
    parser->addParameter(DATA_TEST_FILE_NAME, DATA_TEST_FILE_DEFAULT);
    parser->addParameter(ANN_MIN_ERROR_NAME, ANN_MIN_ERROR_DEFAULT);
    parser->addParameter(ANN_MAX_ITERATIONS_NAME, ANN_MAX_ITERATIONS_DEFAULT);
    parser->addParameter(ANN_SHUFFLE_NAME, ANN_SHUFFLE_DEFAULT);
    parser->addParameter(ANN_SHUFFLE_SEED_NAME, ANN_SHUFFLE_SEED_DEFAULT);
    parser->addParameter(ANN_RPROP_NAME, ANN_RPROP_DEFAULT);
    parser->addParameter(ANN_NORMALIZATION_NAME, ANN_NORMALIZATION_DEFAULT);
    parser->addParameter(ANN_LITERATURE_NAME, ANN_LITERATURE_DEFAULT);
    parser->addParameter(ANN_ACTIVATION_TYPE_NAME, ANN_ACTIVATION_TYPE_DEFAULT);
    parser->addParameter(ANN_A_VALUE_NAME, ANN_A_VALUE_DEFAULT);
    parser->addParameter(ANN_B_VALUE_NAME, ANN_B_VALUE_DEFAULT);
    parser->addParameter(ANN_LEARNING_CONSTANT_NAME, ANN_LEARNING_CONSTANT_DEFAULT);
    parser->addParameter(ANN_MOMENTUM_CONSTANT_NAME, ANN_MOMENTUM_CONSTANT_DEFAULT);
    parser->addParameter(ANN_VALIDATION_MAXIT_NAME, ANN_VALIDATION_MAXIT_DEFAULT);
    parser->addParameter(ANN_ADJUST_FREQUENCY_NAME, ANN_ADJUST_FREQUENCY_DEFAULT);
    parser->addParameter(ANN_ADJUST_CONSTANTS_NAME, ANN_ADJUST_CONSTANTS_DEFAULT);
    parser->addParameter(ANN_ADJUST_FACTOR_NAME, ANN_ADJUST_FACTOR_DEFAULT);
    parser->addParameter(ANN_ADJUST_MIN_NAME, ANN_ADJUST_MIN_DEFAULT);
    parser->addParameter(ANN_ADJUST_MAX_NAME, ANN_ADJUST_MAX_DEFAULT);
    parser->addParameter(ANN_RHO_PLUS_NAME, ANN_RHO_PLUS_DEFAULT);
    parser->addParameter(ANN_RHO_MINUS_NAME, ANN_RHO_MINUS_DEFAULT);
    parser->addParameter(ANN_SEED_NAME, ANN_SEED_DEFAULT);
    parser->addParameter(ANN_RANDOM_BOUND_NAME, ANN_RANDOM_BOUND_DEFAULT);
    parser->addParameter(ANN_LAMBDA_NAME, ANN_LAMBDA_DEFAULT);
    parser->addParameter(ANN_W_ZERO_NAME, ANN_W_ZERO_DEFAULT);
    parser->addParameter(ANN_REGULARIZATION_NAME, ANN_REGULARIZATION_DEFAULT);
    parser->addParameter(ANN_POWER_NAME, ANN_POWER_DEFAULT);
    parser->addParameter(ANN_HIDDEN_NAME, ANN_HIDDEN_DEFAULT);
    parser->addParameter(ANN_OUTPUT_NAME, ANN_OUTPUT_DEFAULT);
    parser->addParameter(ANN_PERIOD_NAME, ANN_PERIOD_DEFAULT);
    parser->addParameter(ANN_SCORES_PER_RES_NAME, ANN_SCORES_PER_RES_DEFAULT);
    parser->addParameter(ANN_WINDOW_SIZE_NAME, ANN_WINDOW_SIZE_DEFAULT);
    parser->addParameter(LOAD_PARAMETERS_NAME, LOAD_PARAMETERS_DEFAULT);
    parser->addParameter(LOAD_FILE_NAME, LOAD_FILE_DEFAULT);
    parser->addParameter(TESTMODE_NAME, TESTMODE_DEFAULT);
    parser->addParameter(ASSIGNMENT_MODE_TEST_NAME, ASSIGNMENT_MODE_TEST_DEFAULT);
    parser->addParameter(CLUSTER_USE_NAME, CLUSTER_USE_DEFAULT);
    parser->addParameter(CLUSTER_FILE_NAME, CLUSTER_FILE_DEFAULT);
    parser->addParameter(EXTENDED_NAME, EXTENDED_DEFAULT);
    parser->addParameter(SCORES_WRITE_NAME, SCORES_WRITE_DEFAULT);
    parser->addParameter(SCORES_FILE_NAME, SCORES_FILE_DEFAULT);
    parser->addParameter(SCORES_FILE_PRED_NAME, SCORES_FILE_PRED_DEFAULT);

    std::ifstream cfgFile(CONFIG_FILE_NAME.c_str(), std::fstream::in);
    if (!cfgFile.is_open())
      parser->writeParameterFile(CONFIG_FILE_NAME);
    else
      cfgFile.close();

    parser->readParameters(CONFIG_FILE_NAME);

    double minError = atof((parser->getParameterValue(ANN_MIN_ERROR_NAME)).c_str());
    unsigned maxIterations = atoi((parser->getParameterValue(ANN_MAX_ITERATIONS_NAME)).c_str());
    bool shuffle = (strcmp((parser->getParameterValue(ANN_SHUFFLE_NAME)).c_str(), "1") == 0);
    bool useRprop = (strcmp((parser->getParameterValue(ANN_RPROP_NAME)).c_str(), "1") == 0);
    bool useNormalization = (strcmp((parser->getParameterValue(ANN_NORMALIZATION_NAME)).c_str(), "1") == 0);
    bool useLiteratureData = (strcmp((parser->getParameterValue(ANN_LITERATURE_NAME)).c_str(), "1") == 0);
    bool testmode = (strcmp((parser->getParameterValue(TESTMODE_NAME)).c_str(), "1") == 0);
    bool loadParameters = (strcmp((parser->getParameterValue(LOAD_PARAMETERS_NAME)).c_str(), "1") == 0);
    bool useExtendedModel = (strcmp((parser->getParameterValue(EXTENDED_NAME)).c_str(), "1") == 0);

    bool writeScores = (strcmp((parser->getParameterValue(SCORES_WRITE_NAME)).c_str(), "1") == 0);
    std::string scoreFileName = parser->getParameterValue(SCORES_FILE_NAME);
    std::string scoreFilePredName = parser->getParameterValue(SCORES_FILE_PRED_NAME);

    ActivationType activationType = toActivationType(parser->getParameterValue(ANN_ACTIVATION_TYPE_NAME));
    double aValue = atof((parser->getParameterValue(ANN_A_VALUE_NAME)).c_str());
    double bValue = atof((parser->getParameterValue(ANN_B_VALUE_NAME)).c_str());
    double learningConstant = atof((parser->getParameterValue(ANN_LEARNING_CONSTANT_NAME)).c_str());
    double momentumConstant = atof((parser->getParameterValue(ANN_MOMENTUM_CONSTANT_NAME)).c_str());
    unsigned validationMaxIt = atoi((parser->getParameterValue(ANN_VALIDATION_MAXIT_NAME)).c_str());
    unsigned adjustFrequency = atoi((parser->getParameterValue(ANN_ADJUST_FREQUENCY_NAME)).c_str());
    bool adjustConstants = (strcmp((parser->getParameterValue(ANN_ADJUST_CONSTANTS_NAME)).c_str(), "1") == 0);
    double adjustFactor = atof((parser->getParameterValue(ANN_ADJUST_FACTOR_NAME)).c_str());
    double adjustMin = atof((parser->getParameterValue(ANN_ADJUST_MIN_NAME)).c_str());
    double adjustMax = atof((parser->getParameterValue(ANN_ADJUST_MAX_NAME)).c_str());
    double rhoPlus = atof((parser->getParameterValue(ANN_RHO_PLUS_NAME)).c_str());
    double rhoMinus = atof((parser->getParameterValue(ANN_RHO_MINUS_NAME)).c_str());
    unsigned annSeed = atoi((parser->getParameterValue(ANN_SEED_NAME)).c_str());
    unsigned annShuffleSeed = atoi((parser->getParameterValue(ANN_SHUFFLE_SEED_NAME)).c_str());
    double annBound = atof((parser->getParameterValue(ANN_RANDOM_BOUND_NAME)).c_str());
    double lambda = atof((parser->getParameterValue(ANN_LAMBDA_NAME)).c_str());
    double w_zero = atof((parser->getParameterValue(ANN_W_ZERO_NAME)).c_str());
    RegularizationType regularization = toRegularizationType(parser->getParameterValue(ANN_REGULARIZATION_NAME));
    unsigned power = atoi((parser->getParameterValue(ANN_POWER_NAME)).c_str());
    unsigned period = atoi((parser->getParameterValue(ANN_PERIOD_NAME)).c_str());
    unsigned hiddenNodes = atoi((parser->getParameterValue(ANN_HIDDEN_NAME)).c_str());
    unsigned outputNodes = atoi((parser->getParameterValue(ANN_OUTPUT_NAME)).c_str()); //corresponds to the number of classes
    unsigned scoresPerRes = atoi((parser->getParameterValue(ANN_SCORES_PER_RES_NAME)).c_str());
    unsigned windowSize = atoi((parser->getParameterValue(ANN_WINDOW_SIZE_NAME)).c_str());

    bool useClusters = (strcmp((parser->getParameterValue(CLUSTER_USE_NAME)).c_str(), "1") == 0);
    std::string clusterFile = parser->getParameterValue(CLUSTER_FILE_NAME);

    KMeansClusterer* clusterer = NULL;
    unsigned numberOfClusters = 0;
    if (useClusters)
    {
      clusterer = new KMeansClusterer();
      clusterer->readPrototypes(clusterFile);
    }

    t_assignment_mode aMode = toAssignmentMode(parser->getParameterValue(ASSIGNMENT_MODE_TEST_NAME));

    std::string paramFileName = parser->getParameterValue(LOAD_FILE_NAME);

    unsigned inputDim = scoresPerRes * windowSize + 2;
    unsigned numOfInputs = inputDim + numberOfClusters;

    std::ofstream wrongRecallLog(WRONG_RECALL_LOG_NAME.c_str(), std::fstream::out | std::fstream::trunc);
    if (!wrongRecallLog.is_open())
    {
      std::cerr << "ERROR: Could not open file \"" << WRONG_RECALL_LOG_NAME << "\" for writing!" << std::endl;
      exit(1);
    }

    std::ofstream wrongPredictionLog(WRONG_PREDICTION_LOG_NAME.c_str(), std::fstream::out | std::fstream::trunc);
    if (!wrongPredictionLog.is_open())
    {
      std::cerr << "ERROR: Could not open file \"" << WRONG_PREDICTION_LOG_NAME << "\" for writing!" << std::endl;
      exit(1);
    }

    time_t t1 = time(NULL);
    std::cout << "Reading learn dataset...";
    std::cout.flush();

    std::string learnFile = parser->getParameterValue(DATA_LEARN_FILE_NAME);
    DatasetType* dataset = readDataset(learnFile, scoresPerRes, outputNodes, true, aMode);
    unsigned learnSize = getSize(dataset);
    double** learnScores = new double*[learnSize];
    double** learnExpected = new double*[learnSize];
    char* learnExpectedClass = new char[learnSize];
    unsigned* learnTrueStr = new unsigned[learnSize];
    unsigned* learnProtIDs = new unsigned[learnSize];
    bool* learnEvaluate = new bool[learnSize];

    unsigned currentIndex = 0;
    for (unsigned protID = 0; protID < dataset->size; protID++)
    {
      for (unsigned resID = 0; resID < dataset->data[protID].size; resID++)
      {
        double* inputs = getSubSequence(dataset, protID, resID, windowSize, scoresPerRes);
        double* distances = NULL;
        if (useClusters) distances = clusterer->getDistances(dataset->data[protID].inputs[resID]);
        learnScores[currentIndex] = new double[numOfInputs];
        for (unsigned i = inputDim; i < numOfInputs; i++) learnScores[currentIndex][i] = distances[i - inputDim];
        learnExpected[currentIndex] = new double[outputNodes];
        for (unsigned i = 0; i < outputNodes; i++)
        {
          learnExpected[currentIndex][i] = dataset->data[protID].outputs[resID][i];
          if (learnExpected[currentIndex][i] == ANN::LITERATURE_POSITIVE_EXPECTED) learnTrueStr[currentIndex] = i;
        }
        learnProtIDs[currentIndex] = protID;
        learnEvaluate[currentIndex] = dataset->data[protID].evaluate[resID];
        learnExpectedClass[currentIndex] = dataset->data[protID].expectedClass[resID];
        delete [] inputs;
        if (distances != NULL) delete [] distances;
        currentIndex++;
      }
    }


    std::cout << "done." << std::endl;

    double** testScores = NULL;
    double** testExpected = NULL;
    char* testExpectedClass = NULL;
    unsigned* testTrueStr = NULL;
    unsigned* testProtIDs = NULL;
    DatasetType* testset = NULL;
    unsigned testSize = 0;
    bool* testEvaluate = NULL;
    if (testmode)
    {
      std::string testFile = parser->getParameterValue(DATA_TEST_FILE_NAME);
      testset = readDataset(testFile, scoresPerRes, outputNodes, true, aMode);
      testSize = getSize(testset);
      testScores = new double*[testSize];
      testExpected = new double*[testSize];
      testExpectedClass = new char[testSize];
      testTrueStr = new unsigned[testSize];
      testProtIDs = new unsigned[testSize];
      testEvaluate = new bool[testSize];

      unsigned currentIndex2 = 0;
      for (unsigned protID = 0; protID < testset->size; protID++)
      {
        for (unsigned resID = 0; resID < testset->data[protID].size; resID++)
        {
          double* inputs = getSubSequence(testset, protID, resID, windowSize, scoresPerRes);
          double* distances = NULL;
          if (useClusters) distances = clusterer->getDistances(testset->data[protID].inputs[resID]);
          testScores[currentIndex2] = new double[numOfInputs];
          for (unsigned i = 0; i < inputDim; i++) testScores[currentIndex2][i] = inputs[i];
          for (unsigned i = inputDim; i < numOfInputs; i++) testScores[currentIndex2][i] = distances[i - inputDim];
          //testScores[currentIndex2] = getSubSequence(testset, protID, resID, windowSize, scoresPerRes);
          testExpected[currentIndex2] = new double[outputNodes];
          for (unsigned i = 0; i < outputNodes; i++)
          {
            testExpected[currentIndex2][i] = testset->data[protID].outputs[resID][i];
            if (testExpected[currentIndex2][i] == ANN::LITERATURE_POSITIVE_EXPECTED) testTrueStr[currentIndex2] = i;
          }
          testProtIDs[currentIndex2] = protID;
          testEvaluate[currentIndex2] = testset->data[protID].evaluate[resID];
          testExpectedClass[currentIndex2] = testset->data[protID].expectedClass[resID];
          delete [] inputs;
          if (distances != NULL) delete [] distances;
          currentIndex2++;
        }
      }
    }

    time_t t2 = time(NULL);
    std::cout << "Setting up ANN...";
    std::cout.flush();
    // set up learning
    ANN* myAnn = new ANN(numOfInputs, hiddenNodes, outputNodes, annSeed, annBound);
    if (useLiteratureData)
      myAnn->useLiteratureValues();
    else
    {
      myAnn->setFuncType(activationType);
      myAnn->setA(aValue);
      myAnn->setB(bValue);
    }
    myAnn->setLearningConstant(learningConstant);
    myAnn->setMomentumConstant(momentumConstant);
    myAnn->setValidationMaxIt(validationMaxIt);
    myAnn->setShuffleSet(shuffle);
    myAnn->setAdjustFrequency(adjustFrequency);
    myAnn->setAdjustConstants(adjustConstants);
    myAnn->setAdjustFactor(adjustFactor);
    myAnn->setAdjustMin(adjustMin);
    myAnn->setAdjustMax(adjustMax);
    myAnn->setRhoPlus(rhoPlus);
    myAnn->setRhoMinus(rhoMinus);
    myAnn->setLambda(lambda);
    myAnn->setWZero(w_zero);
    myAnn->setRegularization(regularization);
    myAnn->setPower(power);
#ifdef ANN_DEBUG_MODE
    myAnn->setPeriod(period);
#endif

    if (!loadParameters) myAnn->initLearnData(learnScores, learnExpected, learnSize, false);
    if (!loadParameters && useNormalization) myAnn->normalizeData();
    if (loadParameters) myAnn->readParameters(paramFileName);

    std::ostringstream paramOutFileName;
    paramOutFileName << /*PARAMETER_PATH << PATH_SEPARATOR <<*/ ANN_PREFIX << "_" << numOfInputs << "_" << hiddenNodes << "_" << outputNodes;
    paramOutFileName << FILE_SUFFIX;
    std::string outfileName = paramOutFileName.str();

    std::cout << "done." << std::endl;
    time_t t3 = time(NULL);
    if (!loadParameters)
    {
      std::cout << "Iterations = " << maxIterations << std::endl;
      std::cout << "Learning...";
      std::cout.flush();
      if (useRprop)
        myAnn->learn2(maxIterations, minError, annShuffleSeed);
      else
        myAnn->learn(maxIterations, minError, annShuffleSeed);

      myAnn->writeParameters(outfileName);
      //myAnn->readParameters(outfileName);

      std::cout << "done." << std::endl;
    }
    time_t t4 = time(NULL);

#ifdef ANN_DEBUG_MODE
    if (period != 0)
    {
      unsigned iterations = 0;
      std::ofstream logFile(LOG_NAME.c_str(), std::fstream::out | std::fstream::trunc);
      if (!logFile.is_open())
      {
        std::cerr << "Could not open file \"" << LOG_NAME << "\" for writing!" << std::endl;
        exit(1);
      }
      logFile << "iteration\trecall";
      if (testmode)
        logFile << "\tprediction";
      logFile << std::endl;
      std::cout << "Creating iteration log...";
      std::cout.flush();
      while (iterations < maxIterations)
      {
        std::ostringstream filename;
        filename << ANN::WRITE_PREFIX_STRING << iterations << ANN::WRITE_SUFFIX_STRING;
        std::ifstream check(filename.str().c_str(), std::fstream::in | std::fstream::binary);
        if (!check.is_open())
          break;
        else
          check.close();
       std::string filenameStr = filename.str();
       myAnn->readParameters(filenameStr);

       double correctRec = 0;
       for (unsigned i = 0; i < learnSize; i++)
       {
         const double* outputActivation = myAnn->getActivation(learnScores[i]);
         double max = outputActivation[0];
         unsigned maxStr = 0;
         for (unsigned k = 1; k < outputNodes; k++)
         {
           if (outputActivation[k] >= max)
           {
             max = outputActivation[k];
             maxStr = k;
           }
         }
         if (maxStr == learnTrueStr[i]) correctRec++;
       }
       double perfRec = 100.0 * correctRec / ((double)learnSize);

       double perfPred = 0;
       if (testmode)
       {
         double correctPred = 0;
         for (unsigned i = 0; i < testSize; i++)
         {
           const double* outputActivation = myAnn->getActivation(testScores[i]);
           double max = outputActivation[0];
           unsigned maxStr = 0;
           for (unsigned k = 1; k < outputNodes; k++)
           {
             if (outputActivation[k] >= max)
             {
               max = outputActivation[k];
               maxStr = k;
             }
           }
           if (maxStr == testTrueStr[i]) correctPred++;
         }
         perfPred = 100.0 * correctPred / ((double)testSize);
       }

       logFile << iterations << '\t' << perfRec;
       if (testmode)
         logFile << '\t' << perfPred;
       logFile << std::endl;
       iterations += period;
     }

     if (iterations >= maxIterations)
     {
       std::ostringstream filename;
       filename << ANN::WRITE_PREFIX_STRING << (maxIterations) << ANN::WRITE_SUFFIX_STRING;
       std::ifstream check(filename.str().c_str(), std::fstream::in | std::fstream::binary);
       if (check.is_open())
       {
         check.close();
         std::string filenameStr = filename.str();
         myAnn->readParameters(filenameStr);

         double correctRec = 0;
         for (unsigned i = 0; i < learnSize; i++)
         {
           const double* outputActivation = myAnn->getActivation(learnScores[i]);
           double max = outputActivation[0];
           unsigned maxStr = 0;
           for (unsigned k = 1; k < outputNodes; k++)
           {
             if (outputActivation[k] >= max)
             {
               max = outputActivation[k];
               maxStr = k;
             }
           }
           if (maxStr == learnTrueStr[i]) correctRec++;
         }
         double perfRec = 100.0 * correctRec / ((double)learnSize);

         double perfPred = 0;
         if (testmode)
         {
           double correctPred = 0;
           for (unsigned i = 0; i < testSize; i++)
           {
             const double* outputActivation = myAnn->getActivation(testScores[i]);
             double max = outputActivation[0];
             unsigned maxStr = 0;
             for (unsigned k = 1; k < outputNodes; k++)
             {
               if (outputActivation[k] >= max)
               {
                 max = outputActivation[k];
                 maxStr = k;
               }
             }
             if (maxStr == testTrueStr[i]) correctPred++;
           }
           perfPred = 100.0 * correctPred / ((double)testSize);
         }

         logFile << iterations << '\t' << perfRec;
         if (testmode)
           logFile << '\t' << perfPred;
         logFile << std::endl;
         iterations += period;
       }
    }
    std::cout << "done." << std::endl;
    logFile.close();
  }
  if (!loadParameters)
    myAnn->readParameters(outfileName);
  else
    myAnn->readParameters(paramFileName);
#endif
  std::cout << "Performing recall test..."<< std::endl;
  double correct = 0;
  double correctProt = 0;
  unsigned protSize = 0;
  unsigned learnSamples = 0;
  unsigned lastRecallID = learnProtIDs[0];
  double* correctPerClass = new double[outputNodes];
  double* totalPerClass = new double[outputNodes];
  double* correctPerClassProt = new double[outputNodes];
  double* totalPerClassProt = new double[outputNodes];
  double** globalRecallConfusionMatrix = new double*[outputNodes];
  double** localRecallConfusionMatrix = new double*[outputNodes];
  std::ofstream scoreFile;
  if (writeScores)
  {
    scoreFile.open(scoreFileName.c_str(), std::fstream::out | std::fstream::trunc);
    if (!scoreFile.is_open())
    {
      std::cerr << "Could not open file \"" << scoreFileName << "\" for writing!" << std::endl;
      exit(1);
    }
    scoreFile << STAGE2_CHAIN_SEPARATOR << std::endl;
  }
  std::ofstream rawRecallFile("recall.raw", std::fstream::out | std::fstream::trunc);
  if (!scoreFile.is_open())
  {
    std::cerr << "Could not open raw recall data file for writing!" << std::endl;
    exit(1);
  }
  for (unsigned i = 0; i < outputNodes; i++)
  {
    correctPerClass[i] = 0;
    totalPerClass[i] = 0;
    correctPerClassProt[i] = 0;
    totalPerClassProt[i] = 0;
    globalRecallConfusionMatrix[i] = new double[outputNodes];
    localRecallConfusionMatrix[i] = new double[outputNodes];
    for (unsigned j = 0; j < outputNodes; j++)
    {
      globalRecallConfusionMatrix[i][j] = 0.0;
      localRecallConfusionMatrix[i][j] = 0.0;
    }
  }
  std::string structureStr = "";
  for (unsigned i = 0; i < learnSize; i++)
  {
    if (!learnEvaluate[i]) continue;
    if (lastRecallID != learnProtIDs[i])
    {
      double perfRecallProt = 100.0 * correctProt / ((double)protSize);
      double* perClassPerfProt = new double[outputNodes];
      for (unsigned c = 0; c < outputNodes; c++)
        perClassPerfProt[c] = 100.0 * correctPerClassProt[c] / totalPerClassProt[c];
      double lgmcc = computeCorrelationCoefficient(localRecallConfusionMatrix, outputNodes);
      std::cout << "=> Recall for protein with ID " << lastRecallID << ":" << std::endl;
      for (unsigned c = 0; c < outputNodes; c++)
        std::cout << getClassName(c, useExtendedModel) << " recall: " << perClassPerfProt[c] << "% (" << correctPerClassProt[c] << "/" << totalPerClassProt[c] << ")" << std::endl;
      std::cout << "Recall: " << perfRecallProt << "% (" << correctProt << "/" << protSize << ")" << std::endl;
      std::cout << "Generalized MCC: " << lgmcc << std::endl;
      std::cout << "Structure:" << std::endl << "str >> " << structureStr << std::endl;
      std::cout << "============================================================================================================" << std::endl;
      for (unsigned k = 0; k < outputNodes; k++)
      {
        correctPerClassProt[k] = 0.0;
        totalPerClassProt[k] = 0.0;
        for (unsigned j = 0; j < outputNodes; j++)
          localRecallConfusionMatrix[k][j] = 0.0;
      }
      protSize = 0;
      correctProt = 0;
      lastRecallID = learnProtIDs[i];
      structureStr = "";
      if (writeScores) scoreFile << STAGE2_CHAIN_SEPARATOR << std::endl;
      delete [] perClassPerfProt;
    }
    const double* outputActivation = myAnn->getActivation(learnScores[i]);
    double max = outputActivation[0];
    if (writeScores) scoreFile << outputActivation[0] << ' ';
    unsigned maxStr = 0;
    for (unsigned k = 1; k < outputNodes; k++)
    {
      if (writeScores) scoreFile << outputActivation[k] << ' ';
      if (outputActivation[k] >= max)
      {
        max = outputActivation[k];
        maxStr = k;
      }
    }
    if (writeScores) scoreFile << learnTrueStr[i] << std::endl;
    rawRecallFile << (maxStr == learnTrueStr[i]) << std::endl;
    if (maxStr == learnTrueStr[i])
    {
      correct++;
      correctProt++;
      correctPerClass[learnTrueStr[i]]++;
      correctPerClassProt[learnTrueStr[i]]++;
    }
    else
    {
      wrongRecallLog << learnProtIDs[i] << " " << protSize << " " << maxStr;
      if (learnExpectedClass[i] != '\0') wrongRecallLog << " " << learnExpectedClass[i];
      wrongRecallLog << std::endl;
    }
    structureStr += getClassLetter(maxStr, useExtendedModel);
    protSize++;
    totalPerClass[learnTrueStr[i]]++;
    totalPerClassProt[learnTrueStr[i]]++;
    globalRecallConfusionMatrix[maxStr][learnTrueStr[i]]++;
    localRecallConfusionMatrix[maxStr][learnTrueStr[i]]++;
    learnSamples++;
  }
  double perfRecallProt = 100.0 * correctProt / ((double)protSize);
  double* perClassPerfProt = new double[outputNodes];
  for (unsigned c = 0; c < outputNodes; c++)
    perClassPerfProt[c] = 100.0 * correctPerClassProt[c] / totalPerClassProt[c];
  double lgmcc = computeCorrelationCoefficient(localRecallConfusionMatrix, outputNodes);
  std::cout << "=> Recall for protein with ID " << lastRecallID << ":" << std::endl;
  for (unsigned c = 0; c < outputNodes; c++)
    std::cout << getClassName(c, useExtendedModel) << " recall: " << perClassPerfProt[c] << "% (" << correctPerClassProt[c] << "/" << totalPerClassProt[c] << ")" << std::endl;
  std::cout << "Recall: " << perfRecallProt << "% (" << correctProt << "/" << protSize << ")" << std::endl;
  std::cout << "Generalized MCC: " << lgmcc << std::endl;
  std::cout << "Structure:" << std::endl << "str >> " << structureStr << std::endl;
  std::cout << "============================================================================================================" << std::endl;
  delete [] perClassPerfProt;
  double perf = 100.0 * correct / ((double)learnSamples);
  double* recallPerf = new double[outputNodes];
  for (unsigned i = 0; i < outputNodes; i++)
    recallPerf[i] = 100.0 * correctPerClass[i] / totalPerClass[i];
  double gmcc = computeCorrelationCoefficient(globalRecallConfusionMatrix, outputNodes);
  std::cout << "=> Overall recall:" << std::endl;
  for(unsigned i = 0; i < outputNodes; i++)
    std::cout << getClassName(i, useExtendedModel) << " recall: " << recallPerf[i] << "% (" << correctPerClass[i] << "/" << totalPerClass[i] << ")" << std::endl;
  std::cout << "Recall: " << perf << "% (" << correct << "/" << learnSamples << ")" << std::endl;
  std::cout << "Generalized MCC: " << gmcc << std::endl;
  if (writeScores) scoreFile.close();
  rawRecallFile.close();

  time_t t5 = time(NULL);

  if (testmode)
  {
    structureStr = "";
    std::cout << "Performing prediction test..." << std::endl;
    double correctPred = 0.0;
    double correctPredProt = 0.0;
    unsigned testSamples = 0;
    protSize = 0;
    unsigned lastPredID = testProtIDs[0];
    double** globalPredConfusionMatrix = new double*[outputNodes];
    double** localPredConfusionMatrix = new double*[outputNodes];
    std::ofstream scoreFilePred;
    if (writeScores)
    {
      scoreFilePred.open(scoreFilePredName.c_str(), std::fstream::out | std::fstream::trunc);
      if (!scoreFilePred.is_open())
      {
        std::cerr << "Could not open file \"" << scoreFilePredName << "\" for writing!" << std::endl;
        exit(1);
      }
      scoreFilePred << STAGE2_CHAIN_SEPARATOR << std::endl;
    }
    std::ofstream rawPredictionFile("prediction.raw", std::fstream::out | std::fstream::trunc);
    if (!rawPredictionFile.is_open())
    {
      std::cerr << "Could not open raw prediction file for writing!" << std::endl;
      exit(1);
    }
    for (unsigned i = 0; i < outputNodes; i++)
    {
      correctPerClass[i] = 0;
      totalPerClass[i] = 0;
      correctPerClassProt[i] = 0;
      totalPerClassProt[i] = 0;
      globalPredConfusionMatrix[i] = new double[outputNodes];
      localPredConfusionMatrix[i] = new double[outputNodes];
      for (unsigned j = 0; j < outputNodes; j++)
      {
        globalPredConfusionMatrix[i][j] = 0.0;
        localPredConfusionMatrix[i][j] = 0.0;
      }
    }
    for (unsigned i = 0; i < testSize; i++)
    {
      if (!testEvaluate[i]) continue;
      if (lastPredID != testProtIDs[i])
      {
        double perfPredProt = 100.0 * correctPredProt / ((double)protSize);
        double* perClassPredPerfProt = new double[outputNodes];
        for (unsigned c = 0; c < outputNodes; c++)
          perClassPredPerfProt[c] = 100.0 * correctPerClassProt[c] / totalPerClassProt[c];
        double lgmcc = computeCorrelationCoefficient(localPredConfusionMatrix, outputNodes);
        std::cout << "=> Prediction for protein with ID " << lastPredID << ":" << std::endl;
        for (unsigned c = 0; c < outputNodes; c++)
          std::cout << getClassName(c, useExtendedModel) << " prediction: " << perClassPredPerfProt[c] << "% (" << correctPerClassProt[c] << "/" << totalPerClassProt[c] << ")" << std::endl;
        std::cout << "Prediction: " << perfPredProt << "% (" << correctPredProt << "/" << protSize << ")" << std::endl;
        std::cout << "Generalized MCC: " << lgmcc << std::endl;
        std::cout << "Structure:" << std::endl << "str >> " << structureStr << std::endl;
        std::cout << "============================================================================================================" << std::endl;
        for (unsigned k = 0; k < outputNodes; k++)
        {
          correctPerClassProt[k] = 0;
          totalPerClassProt[k] = 0;
          for (unsigned j = 0; j < outputNodes; j++)
            localPredConfusionMatrix[k][j] = 0.0;
        }
        protSize = 0;
        correctPredProt = 0;
        structureStr = "";
        lastPredID = testProtIDs[i];
        delete [] perClassPredPerfProt;
        if (writeScores) scoreFilePred << STAGE2_CHAIN_SEPARATOR << std::endl;
      }
      const double* outputActivation = myAnn->getActivation(testScores[i]);
      double max = outputActivation[0];
      if (writeScores) scoreFilePred << outputActivation[0] << ' ';
      unsigned maxStr = 0;
      for (unsigned k = 1; k < outputNodes; k++)
      {
        if (writeScores) scoreFilePred << outputActivation[k] << ' ';
        if (outputActivation[k] >= max)
        {
          max = outputActivation[k];
          maxStr = k;
        }
      }
      if (writeScores) scoreFilePred << testTrueStr[i] << std::endl;
      rawPredictionFile << (maxStr == testTrueStr[i]) << std::endl;
      if (maxStr == testTrueStr[i])
      {
        correctPred++;
        correctPredProt++;
        correctPerClass[testTrueStr[i]]++;
        correctPerClassProt[testTrueStr[i]]++;
      }
      else
      {
        wrongPredictionLog << testProtIDs[i] << " " << protSize << " " << maxStr;
        if (testExpectedClass[i] != '\0') wrongPredictionLog << " " << testExpectedClass[i];
        wrongPredictionLog << std::endl;
      }
      structureStr += getClassLetter(maxStr,useExtendedModel);
      protSize++;
      totalPerClass[testTrueStr[i]]++;
      totalPerClassProt[testTrueStr[i]]++;
      globalPredConfusionMatrix[maxStr][testTrueStr[i]]++;
      localPredConfusionMatrix[maxStr][testTrueStr[i]]++;
      testSamples++;
    }
    double perfPredProt = 100.0 * correctPredProt / ((double)protSize);
    double* perClassPredPerfProt = new double[outputNodes];
    for (unsigned c = 0; c < outputNodes; c++)
      perClassPredPerfProt[c] = 100.0 * correctPerClassProt[c] / totalPerClassProt[c];
    double lgmcc = computeCorrelationCoefficient(localPredConfusionMatrix, outputNodes);
    std::cout << "=> Prediction for protein with ID " << lastPredID << ":" << std::endl;
    for (unsigned c = 0; c < outputNodes; c++)
      std::cout << getClassName(c,useExtendedModel) << " prediction: " << perClassPredPerfProt[c] << "% (" << correctPerClassProt[c] << "/" << totalPerClassProt[c] << ")" << std::endl;
    std::cout << "Prediction: " << perfPredProt << "% (" << correctPredProt << "/" << protSize << ")" << std::endl;
    std::cout << "Generalized MCC: " << lgmcc << std::endl;
    std::cout << "Structure:" << std::endl << "str >> " << structureStr << std::endl;
    std::cout << "============================================================================================================" << std::endl;
    //delete [] perClassPredPerfProt;
    double perfPred = 100.0 * correctPred / ((double)testSamples);
    double* perClassPredPerf = new double[outputNodes];
    for (unsigned c = 0; c < outputNodes; c++)
      perClassPredPerf[c] = 100.0 * correctPerClass[c] / totalPerClass[c];
    gmcc = computeCorrelationCoefficient(globalPredConfusionMatrix, outputNodes);
    std::cout << "=> Overall prediction:" << std::endl;
    for (unsigned c = 0; c < outputNodes; c++)
      std::cout << getClassName(c, useExtendedModel) << " prediction: " << perClassPredPerf[c] << "% (" << correctPerClass[c] << "/" << totalPerClass[c] << ")" << std::endl;
    std::cout << "Prediction: " << perfPred << "% (" << correctPred << "/" << testSamples << ")" << std::endl;
    std::cout << "Generalized MCC: " << gmcc << std::endl;
    for (unsigned i = 0; i < outputNodes; i++)
    {
      delete [] globalPredConfusionMatrix[i];
      delete [] localPredConfusionMatrix[i];
    }
    delete [] globalPredConfusionMatrix;
    delete [] localPredConfusionMatrix;
    delete [] perClassPredPerfProt;
    if (writeScores) scoreFilePred.close();
    rawPredictionFile.close();
  }

  std::cout << "Freeing memory...";
  std::cout.flush();
  // free allocated memory
  delete [] correctPerClass;
  delete [] totalPerClass;
  delete [] correctPerClassProt;
  delete [] totalPerClassProt;

  for (unsigned i = 0; i < learnSize; i++)
  {
    delete [] learnScores[i];
    delete [] learnExpected[i];
  }
  delete [] learnScores;
  delete [] learnTrueStr;
  delete [] learnExpected;
  delete [] learnProtIDs;
  delete [] learnEvaluate;
  for (unsigned i = 0; i < outputNodes; i++)
  {
    delete [] globalRecallConfusionMatrix[i];
    delete [] localRecallConfusionMatrix[i];
  }
  delete [] globalRecallConfusionMatrix;
  delete [] localRecallConfusionMatrix;
  delete [] recallPerf;

  if (testmode)
  {
    for (unsigned i = 0; i < testSize; i++)
    {
      delete [] testScores[i];
      delete [] testExpected[i];
    }
    delete [] testScores;
    delete [] testTrueStr;
    delete [] testExpected;
    delete [] testProtIDs;
    delete [] testEvaluate;

    freeData(testset);
  }

  delete dParam;
  delete parser;
  delete myAnn;
  if (clusterer != NULL) delete clusterer;
  freeData(dataset);
  log.close();
  wrongRecallLog.close();
  wrongPredictionLog.close();
  std::cout << "done." << std::endl;
  time_t t6 = time(NULL);

  std::cout << "Data loading time: " << (t2 - t1) << " seconds" << std::endl;
  std::cout << "Initialization time: " << (t3 - t2) << " seconds" << std::endl;
  std::cout << "Learning time: " << (t4 - t3) << " seconds" << std::endl;
  std::cout << "Recall time: " << (t5 - t4) << " seconds" << std::endl;
  std::cout << "Memory freeing time: " << (t6 - t5) << " seconds" << std::endl;
  }
  catch (ParserException e)
  {
    std::cout << std::endl << "ParserException: " << e.what() << std::endl;
    exit(1);
  }
  catch (ANNException e)
  {
    std::cout << std::endl << "ANNException: " << e.what() << std::endl;
    exit(1);
  }
}
