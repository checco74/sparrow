#include "ANN.h"
#include "time.h"
#include <iostream>
#include <string>
#include <list>
#include <sstream>

const unsigned INPUT_NODES = 3;
const unsigned OUTPUT_NODES = 3;

const std::string CHAIN_SEPARATOR = "###################";

void readDataset(double** input, double** output, std::string filename)
{
  std::string line = "";
  std::ifstream inFile(filename.c_str(), std::fstream::in);
  if (!inFile.is_open())
  {
    std::cerr << "ERROR: Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }

  unsigned i = 0;
  unsigned lineNum = 1;
  while (getline(inFile, line, '\n'))
  {
    if (line != CHAIN_SEPARATOR)
    {
      unsigned currentLineIndex = 0;
      for (unsigned k = 0; k < INPUT_NODES; k++)
      {
        std::ostringstream valStr;
        while (currentLineIndex < line.size() && line[currentLineIndex] != ' ' && line[currentLineIndex] != '\t')
        {
          valStr << line[currentLineIndex];
          currentLineIndex++;
        }
        input[i][k] = atof(valStr.str().c_str());
        while (currentLineIndex < line.size() && (line[currentLineIndex] == ' ' || line[currentLineIndex] == '\t')) currentLineIndex++;
  
        if (currentLineIndex == line.size())
        {
          if (k == INPUT_NODES - 1)
          {
            std::cerr << "ERROR: Missing expected value in file \"" << filename << "\" line " << lineNum << std::endl;
            exit(1);
          }
          else
          {
            std::cerr << "ERROR: Missing value for input entry " << (k+1) <<  " in file \"" << filename << "\" line " << lineNum << std::endl;
            exit(1);
          }
        }
      }


      // now the character at current index should be an character different from a whitespace character
      std::string ending = line.substr(currentLineIndex, line.size() - currentLineIndex);
      unsigned expectedClass = atoi(ending.c_str());
      if (expectedClass >= OUTPUT_NODES)
      {
        std::cerr << "ERROR: Expected value in line " << lineNum << " in file \"" << filename << "\"" << std::endl;
        exit(1);
      }
      for (unsigned k = 0; k < OUTPUT_NODES; k++)
      {
        if (k == expectedClass)
          output[i][k] = ANN::LITERATURE_POSITIVE_EXPECTED;
        else
          output[i][k] = ANN::LITERATURE_NEGATIVE_EXPECTED;
      }

      i++;
    }
    lineNum++;
  }

  inFile.close();
}

void freeData(double** input, double** output, unsigned size)
{
  for (unsigned i = 0; i < size; i++)
  {
    delete [] input[i];
    delete [] output[i];
  }
  delete [] input;
  delete [] output;
}

void countInputLines(unsigned* size, std::string filename)
{
  std::string line = "";
  std::ifstream inFile(filename.c_str(), std::fstream::in);
  if (!inFile.is_open())
  {
    std::cerr << "ERROR: Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }

  unsigned setSize = 0;
  while (getline(inFile, line, '\n'))
  {
    if (line != CHAIN_SEPARATOR) setSize++;
  }
  *size = setSize;

  inFile.close();
}

const char LEARN_FILE_OPTION = 'l';
const char VALIDATION_FILE_OPTION = 'v';
const char PREDICTION_FILE_OPTION = 'p';
const char SHUFFLE_OPTION = 's';
const char RANDOM_OPTION = 'r';
const char MAXIT_OPTION = 'm';
const char MIN_ERROR_OPTION = 'e';
const char TEST_OPTION = 't';
const char HIDDEN_OPTION = 'h';
const char NORMALIZATION_OPTION = 'n';

std::string learnFilename = "learnData.in";
std::string validationFilename = "";
std::string predictionFilename = "";
bool shuffle = false;
unsigned randomSeed = 0;
unsigned hidden = 0;
unsigned maxIt = 500;
double minError = 1e-5;
bool recallTest = false;
bool normalize = false;

void printUsage(char* progName)
{
  std::cout << "Usage:" << progName << " [<options>] <number of hidden nodes>" << std::endl << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "\t-" << LEARN_FILE_OPTION << " <filename>\t\tTells the program to use the specified file to load the training set." << std::endl;
  std::cout << "\t-" << VALIDATION_FILE_OPTION << " <filename>\t\tTells the program to use the specified file to load the validation set." << std::endl;
  std::cout << "\t-" << PREDICTION_FILE_OPTION << " <filename>\t\tTells the program to use the specified file to load the prediction set." << std::endl;
  std::cout << "\t-" << SHUFFLE_OPTION << " \t\t\tTells the program to shuffle the order of the learning set after each iteration." << std::endl;
  std::cout << "\t-" << RANDOM_OPTION << " <unsigned integer>\tTells the program to use the specified number as seed when randomly initializing the synaptic weights of the ANN" << std::endl << "\t\t\t\t(0 means no seed)." << std::endl;
  std::cout << "\t-" << MAXIT_OPTION << " <unsigned integer>\tTells the program to use the specified number as maximum number of iterations for the learning algorithm." << std::endl;
  std::cout << "\t-" << MIN_ERROR_OPTION << " <real number>\tTells the program to stop the learning algorithm if the error change falls below the specified number." << std::endl;
  std::cout << "\t-" << HIDDEN_OPTION << " <unsigned integer>\tTells the program to stop use the specified number as number of hidden nodes for the ANN." << std::endl;
  std::cout << "\t-" << TEST_OPTION << " \t\t\tTells the program to perform a recall test after learning." << std::endl;
  std::cout << "\t-" << NORMALIZATION_OPTION << " \t\t\tTells the program to normalize all input data fed into the network." << std::endl << std::endl;
}

void readCommandLine(unsigned argc, char** argv)
{
  unsigned i = 1;
  while (i < argc)
  {
    if (argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
        case LEARN_FILE_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing learn file name for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          learnFilename = argv[i];
          break;
        }
        case VALIDATION_FILE_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing validation file name for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          validationFilename = argv[i];
          break;
        }
        case PREDICTION_FILE_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing prediction file name for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          predictionFilename = argv[i];
          break;
        }
        case MAXIT_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing maximum number of iterations for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          maxIt = atoi(argv[i]);
          break;
        }
        case MIN_ERROR_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing minimum error change for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          minError = atof(argv[i]);
          break;
        }
        case HIDDEN_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing number of hidden nodes for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          hidden = atoi(argv[i]);
          break;
        }
        case SHUFFLE_OPTION:
          shuffle = true;
          break;
        case RANDOM_OPTION:
        {
          if (i + 1 >= argc)
          {
            std::cerr << "ERROR: Missing seed number for option " << argv[i] << "!" << std::endl;
            printUsage(argv[0]);
            exit(1);
          }
          i++;
          randomSeed = atoi(argv[i]);
          break;
        }
        case TEST_OPTION:
          recallTest = true;
          break;
        case NORMALIZATION_OPTION:
          normalize = true;
          break;
        default:
        {
          std::cerr << "ERROR: Unknown option: " << argv[i] << std::endl;
          printUsage(argv[0]);
          exit(1);
        }
      }
    }
    else
    {
      if (hidden == 0) hidden = atoi(argv[i]);
    }
    i++;
  }
  if (hidden == 0)
  {
    std::cerr << "ERROR: Number of hidden nodes has not been set!" << std::endl;
    exit(1);
  }
}

void prediction(double** input, double** output, unsigned samples, ANN* ann, bool recall)
{
  if (recall)
    std::cout << "Performing recall test...";
  else
    std::cout << "Performing prediction test...";
  std::cout.flush();
  double correct = 0;
  double total = 0;

  time_t t1 = time(NULL);

  for (unsigned i = 0; i < samples; i++)
  {
    double max = -2;
    unsigned expected = 0;
    for (unsigned j = 0; j < OUTPUT_NODES; j++)
    {
      if (output[i][j] > max)
      {
        max = output[i][j];
        expected = j;
      }
    }

    if (!recall && normalize) ann->normalize(input[i]);
    const double* activation = ann->getActivation(input[i]);

    max = -2;
    unsigned predicted = 0;
    for (unsigned j = 0; j < OUTPUT_NODES; j++)
    {
      if (activation[j] > max)
      {
        max = activation[j];
        predicted = j;
      }
    }

    total++;
    if (expected == predicted) correct++;
  }

  time_t t2 = time(NULL);

  std::cout << "done!" << std::endl;
  if (recall)
  {
    std::cout << "Recalled: " << correct << "/" << total << " (" << (100.0*correct/total) << "%)" << std::endl;
    std::cout << "Recall time: " << (t2 - t1) << " seconds" << std::endl;
  }
  else
  {
    std::cout << "Predicted: " << correct << "/" << total << " (" << (100.0*correct/total) << "%)" << std::endl;
    std::cout << "Prediction time: " << (t2 - t1) << " seconds" << std::endl;
  }
}

int main(unsigned argc, char** argv)
{
  if (argc == 1)
  {
    printUsage(argv[0]);
    exit(1);
  }
  readCommandLine(argc, argv);
  std::cout << "Initializing ANN...";
  std::cout.flush();
  ANN* myANN = new ANN(INPUT_NODES, hidden, OUTPUT_NODES, randomSeed);
  std::cout << "done!" << std::endl;

  unsigned samples = 0;
  double** learnData = NULL;
  double** expected = NULL;
  unsigned validationSamples = 0;
  double** validationData = NULL;
  double** validationExpected = NULL;
  unsigned predictionSamples = 0;
  double** predictionData = NULL;
  double** predictionExpected = NULL;

  std::cout << "Reading learning data...";
  std::cout.flush();
  countInputLines(&samples, learnFilename);
  learnData = new double*[samples];
  expected = new double*[samples];
  for (unsigned i = 0; i < samples; i++)
  {
    learnData[i] = new double[INPUT_NODES];
    expected[i] = new double[OUTPUT_NODES];
  }
  readDataset(learnData, expected, learnFilename);
//   std::ofstream tFile("test.out", std::fstream::out | std::fstream::trunc);
//   for (unsigned i = 0; i < samples; i++)
//   {
//     for (unsigned j = 0; j < INPUT_NODES; j++)
//     {
//       tFile << learnData[i][j] << "\t";
//     }
//     for (unsigned j = 0; j < INPUT_NODES; j++)
//     {
//       if (expected[i][j] == ANN::LITERATURE_POSITIVE_EXPECTED)
//         tFile << j << std::endl;
//     }
//   }
//   tFile.close();
  std::cout << samples << " samples read!" << std::endl;
  if (validationFilename != "")
  {
    std::cout << "Reading validation data...";
    std::cout.flush();
    countInputLines(&validationSamples, validationFilename);
    validationData = new double*[validationSamples];
    validationExpected = new double*[validationSamples];
    for (unsigned i = 0; i < validationSamples; i++)
    {
      validationData[i] = new double[INPUT_NODES];
      validationExpected[i] = new double[OUTPUT_NODES];
    }
    readDataset(validationData, validationExpected, validationFilename);
    std::cout << validationSamples << " samples read!" << std::endl;
  }
  if (predictionFilename != "")
  {
    std::cout << "Reading prediction data...";
    std::cout.flush();
    countInputLines(&predictionSamples, predictionFilename);
    predictionData = new double*[predictionSamples];
    predictionExpected = new double*[predictionSamples];
    for (unsigned i = 0; i < predictionSamples; i++)
    {
      predictionData[i] = new double[INPUT_NODES];
      predictionExpected[i] = new double[OUTPUT_NODES];
    }
    readDataset(predictionData, predictionExpected, predictionFilename);
    std::cout << predictionSamples << " samples read!" << std::endl;
  }

  time_t t1 = time(NULL);
  myANN->useLiteratureValues();
  myANN->setShuffleSet(shuffle);
  myANN->setAdjustConstants(true);
  myANN->setAdjustFrequency(50);
  std::cout << "Initializing learning data...";
  std::cout.flush();
  myANN->initLearnData(learnData, expected, samples, false);
  std::cout << "done!" << std::endl;
  if (validationSamples != 0)
  {
    std::cout << "Initializing validation data...";
    std::cout.flush();
    myANN->initValidationData(validationData, validationExpected, validationSamples, false);
    std::cout << "done!" << std::endl;
  }
  if (normalize)
  {
    std::cout << "Normalizing data...";
    std::cout.flush();
    myANN->normalizeData();
    std::cout << "done!" << std::endl;
  }
  std::cout << "Learning...";
  std::cout.flush();
  //myANN->learn(maxIt, minError);
  myANN->learn2(maxIt, minError);
  std::cout << "done!" << std::endl;
  time_t t2 = time(NULL);

  std::cout << "Learning time: " << (t2-t1) << " seconds" << std::endl;

  if (recallTest) prediction(learnData, expected, samples, myANN, true);
  if (predictionSamples != 0) prediction(predictionData, predictionExpected, predictionSamples, myANN, false);

  std::cout << "Freeing allocated memory...";
  std::cout.flush();
  delete myANN;
  freeData(learnData, expected, samples);
  if (validationSamples != 0) freeData(validationData, validationExpected, validationSamples);
  if (predictionSamples != 0) freeData(predictionData, predictionExpected, predictionSamples);
  std::cout << "done!" << std::endl;
}
