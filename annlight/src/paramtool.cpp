#include <iostream>
#include <fstream>
#include <sstream>
#include "DiscriminantParameters.h"
#include "PredictionFile.h"
#include "math.h"

const char DIM_LIN_CHAR = 'l';
const char DIM_QUAD_CHAR = 'q';
const char DIM_STAGE_CHAR = 's';

const char PROC_MEANS_CHAR = 'm';
const char PROC_LEN_CHAR = 'l';
const char PROC_DIFF_CHAR = 'd';
const char PROC_DIST_CHAR = 'D';
const char PROC_RMSD_CHAR = 'r';
const char PROC_VALUES_CHAR = 'v';
const char PROC_COMP_CHAR = 'c';

const std::string DIFF_FILE_NAME = "difference.dat";
const std::string OUTPUT_FILE_EXTENSION = ".dat";

char* progName;

unsigned dim = 0;
bool linear = false;
bool quadratic = false;
bool stage2 = false;

bool means = false;
bool lengths = false;
bool rmsd = false;
bool distances = false;
bool differences = false;
bool values = false;
bool comp = false;

std::string getBoolStr(bool b)
{
  if (b)
    return "yes";
  else
    return "no";
}

void printUsage()
{
  std::cout << "Usage: ./" << progName << " --<dimension options> -<processing options> <filename1> [<filename2>]" << std::endl << std::endl;
  std::cout << "<filename1> and <filename2>(optional) are the names of the files where the parameter sets to process are stored (in binary)" << std::endl << std::endl;
  std::cout << "Dimension options (if using more than one option do not add spaces in between)" << std::endl;
  std::cout << "\t" << DIM_LIN_CHAR << "\tDetermines whether parameters contain linear terms" << std::endl;
  std::cout << "\t" << DIM_QUAD_CHAR << "\tDetermines whether parameters contain quadratic terms" << std::endl;
  std::cout << "\t" << DIM_STAGE_CHAR << "\tDetermines whether parameters are stage 2 parameters" << std::endl;
  std::cout << "Immediatly after the dimension options (no whitespace) should follow an ineteger number representing the window length for the read parameters" << std::endl << std::endl;
  std::cout << "Processing parameters (if using more than one option do not add spaces in between)" << std::endl;
  std::cout << "\t" << PROC_MEANS_CHAR << "\tDetermines whether to display the mean parameter values for the given parameter sets" << std::endl;
  std::cout << "\t" << PROC_LEN_CHAR << "\tDetermines whether to display the length (squareroot of sum of squares) for the given parameter sets" << std::endl;
  std::cout << "\t" << PROC_DIST_CHAR << "\tDetermines whether to display the distance between the two parameter sets" << std::endl;
  std::cout << "\t" << PROC_RMSD_CHAR << "\tDetermines whether to display the RMSD of the two parameter sets (or RMSD with respect to zero vector if only one set is given)" << std::endl;
  std::cout << "\t" << PROC_COMP_CHAR << "\tDetermines whether to display a similarity measure" << std::endl;
  std::cout << "\t" << PROC_DIFF_CHAR << "\tDetermines whether to write the differences between the entries of both sets into \"" << DIFF_FILE_NAME << "\"" << std::endl;
  std::cout << "\t" << PROC_VALUES_CHAR << "\tDetermines whether to write the parameter values in human-readable form into files" << std::endl;
}

void parseProcessingParameters(char* str)
{
  unsigned i = 0;
  while (str[i] != '\0')
  {
    if (str[i] != '-')
    {
      switch (str[i])
      {
        case PROC_MEANS_CHAR:
          means = true;
          break;
        case PROC_LEN_CHAR:
          lengths = true;
          break;
        case PROC_DIFF_CHAR:
          differences = true;
          break;
        case PROC_DIST_CHAR:
          distances = true;
          break;
        case PROC_RMSD_CHAR:
          rmsd = true;
          break;
        case PROC_VALUES_CHAR:
          values = true;
          break;
        case PROC_COMP_CHAR:
          comp = true;
          break;
        default:
          std::cout << "ERROR: Unknown character '" << str[i] << "' in option string \"" << str << "\"!" << std::endl;
          printUsage();
          exit(1);
      }
    }
    i++;
  }

  std::cout << "Processing parameters:" << std::endl;
  std::cout << "\tdisplay mean values: " << getBoolStr(means) << std::endl;
  std::cout << "\tdisplay lengths: " << getBoolStr(lengths) << std::endl;
  std::cout << "\tdisplay distances: " << getBoolStr(distances) << std::endl;
  std::cout << "\tdisplay RMSD: " << getBoolStr(rmsd) << std::endl;
  std::cout << "\tdisplay comparison factor: " << getBoolStr(comp) << std::endl;
  std::cout << "\twrite values: " << getBoolStr(values) << std::endl;
  std::cout << "\twrite differences: " << getBoolStr(differences) << std::endl << std::endl;
}

void parseDimensionParameters(char* str)
{
  unsigned i = 0;
  while (str[i] != '\0')
  {
    if (str[i] != '-')
    {
      if (str[i] == DIM_LIN_CHAR)
        linear = true;
      else if (str[i] == DIM_QUAD_CHAR)
        quadratic = true;
      else if (str[i] == DIM_STAGE_CHAR)
        stage2 = true;
      else if (str[i] >= '0' && str[i] <= '9')
      {
        dim = atoi(&(str[i]));
        break;
      }
      else
      {
        std::cout << "ERROR: Unknown character '" << str[i] << "' in option string \"" << str << "\"!" << std::endl;
        printUsage();
        exit(1);
      }
    }
    i++;
  }

  std::cout << "Dimension parameters: " << std::endl;
  std::cout << "\tdimension: " << dim << std::endl;
  std::cout << "\tlinear terms: " << getBoolStr(linear) << std::endl;
  std::cout << "\tquadratic terms: " << getBoolStr(quadratic) << std::endl;
  std::cout << "\tstage 2: " << getBoolStr(stage2) << std::endl << std::endl;
}

DiscriminantParameters* loadParameters(const std::string& filename)
{
  if (dim == 0)
  {
    std::cout << "ERROR: No input sequence window length specified!" << std::endl;
    printUsage();
    exit(1);
  }
  if (filename != "")
  {
    std::ifstream paramFile(filename.c_str(), std::fstream::in | std::fstream::binary);
    if (!paramFile.is_open())
    {
      std::cout << "Parameter file \"" << filename << "\" could not be opened for reading!" << std::endl;
      exit(1);
    }
    else
    {
      std::cout << "Reading parameters from file \"" << filename << "\"..." << std::endl;
    }

    unsigned vecDim = 20 * dim;
    if (stage2)
      vecDim = PredictionFileEntry::NUMBER_OF_SCORES;
    DiscriminantParameters* param = new DiscriminantParameters(vecDim);

    if (linear)
    {
      for (unsigned i = 0; i < param->getVecDim(); i++)
      {
        double value;
        paramFile.read((char*)(&value), sizeof(double));
        param->set(i, value);
      }
    }

    if (quadratic)
    {
      for (unsigned i = 0; i < param->getMatDim(); i++)
      {
        for (unsigned j = 0; j <= i; j++)
        {
          double value;
          paramFile.read((char*)(&value), sizeof(double));
          if (i != j)
            param->set(i, j, (value / 2.0));
          else
            param->set(i, j, value);
        }
      }
    }

    double val;
    paramFile.read((char*)(&val), sizeof(double));
    param->set(val);

    paramFile.close();

    return param;
  }
  else
    return NULL;
}

void displayMeans(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  double meanVal1 = 0;
  double meanVal2 = 0;
  double paramNum = 0;
  unsigned vecDim = 0;
  unsigned matDim = 0;
  if (param1 != NULL)
  {
    vecDim = param1->getVecDim();
    matDim = param1->getMatDim();
  }
  else if (param2 != NULL)
  {
    vecDim = param2->getVecDim();
    matDim = param2->getMatDim();
  }

  if (linear)
  {
    for (unsigned i = 0; i < vecDim; i++)
    {
      if (param1 != NULL) meanVal1 += param1->get(i);
      if (param2 != NULL) meanVal2 += param2->get(i);
      paramNum += 1;
    }
  }

  if (quadratic)
  {
    for (unsigned i = 0; i < param1->getMatDim(); i++)
    {
      for (unsigned j = i; j <= i; j++)
      {
        if (param1 != NULL) meanVal1 += param1->get(i, j);
        if (param2 != NULL) meanVal2 += param2->get(i, j);
        paramNum += 1;
      }
    }
  }

  if (param1 != NULL) meanVal1 += param1->get();
  if (param2 != NULL) meanVal2 += param2->get();
  paramNum += 1;

  meanVal1 /= paramNum;
  meanVal2 /= paramNum;

  std::cout << "Mean values:" << std::endl;
  if (param1 != NULL) std::cout << "\t" << filename1 << ": " << meanVal1 << std::endl;
  if (param2 != NULL) std::cout << "\t" << filename2 << ": " << meanVal2 << std::endl;
  std::cout << std::endl;
}

void displayLengths(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  double length1 = 0;
  double length2 = 0;
  unsigned vecDim = 0;
  unsigned matDim = 0;
  if (param1 != NULL)
  {
    vecDim = param1->getVecDim();
    matDim = param1->getMatDim();
  }
  else if (param2 != NULL)
  {
    vecDim = param2->getVecDim();
    matDim = param2->getMatDim();
  }

  if (linear)
  {
    for (unsigned i = 0; i < vecDim; i++)
    {
      if (param1 != NULL) length1 += pow(param1->get(i), 2);
      if (param2 != NULL) length2 += pow(param2->get(i), 2);
    }
  }

  if (quadratic)
  {
    for (unsigned i = 0; i < param1->getMatDim(); i++)
    {
      for (unsigned j = i; j <= i; j++)
      {
        if (param1 != NULL) length1 += pow(param1->get(i, j), 2);
        if (param2 != NULL) length2 += pow(param2->get(i, j), 2);
      }
    }
  }

  if (param1 != NULL) length1 += pow(param1->get(), 2);
  if (param2 != NULL) length2 += pow(param2->get(), 2);

  if (param1 != NULL) length1 = sqrt(length1);
  if (param2 != NULL) length2 = sqrt(length2);

  std::cout << "Lengths:" << std::endl;
  if (param1 != NULL) std::cout << "\t" << filename1 << ": " << length1 << std::endl;
  if (param2 != NULL) std::cout << "\t" << filename2 << ": " << length2 << std::endl;
  std::cout << std::endl;
}

void displayRMSDs(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  DiscriminantParameters* first = NULL;
  DiscriminantParameters* second = NULL;

  if (param1 != NULL && param2 != NULL)
  {
    first = param1;
    second = param2;
  }
  else
  {
    if (param1 != NULL)
      first = param1;
    if (param2 != NULL)
      first = param2;
    second = new DiscriminantParameters(first->getVecDim());
  }

  double sum = 0;
  double paramNum = 0;
  double vecDim = first->getVecDim();
  double matDim = first->getMatDim();

  if (linear)
  {
    for (unsigned i = 0; i < vecDim; i++)
    {
      double delta = first->get(i) - second->get(i);
      sum += (delta * delta);
      paramNum += 1;
    }
  }

  if (quadratic)
  {
    for (unsigned i = 0; i < matDim; i++)
    {
      for (unsigned j = i; j <= i; j++)
      {
        double delta = first->get(i, j) - second->get(i, j);
        sum += (delta * delta);
        paramNum += 1;
      }
    }
  }

  double delta = first->get() - second->get();
  sum += (delta * delta);
  paramNum += 1;

  std::cout << "RMSD: " << sqrt(sum/paramNum) << std::endl << std::endl;

  if (param1 == NULL || param2 == NULL) delete second;
}

void displayDistances(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  if (param1 != NULL && param2 != NULL)
  {
    double sum = 0;
    double vecDim = param1->getVecDim();
    double matDim = param1->getMatDim();

    if (linear)
    {
      for (unsigned i = 0; i < vecDim; i++)
      {
        double delta = param1->get(i) - param2->get(i);
        sum += (delta * delta);
      }
    }

    if (quadratic)
    {
      for (unsigned i = 0; i < matDim; i++)
      {
        for (unsigned j = i; j <= i; j++)
        {
          double delta = param1->get(i, j) - param2->get(i, j);
          sum += (delta * delta);
        }
      }
    }

    double delta = param1->get() - param2->get();
    sum += (delta * delta);

    std::cout << "Distance: " << sqrt(sum) << std::endl << std::endl;
  }
}

void displayComparison(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  if (param1 != NULL && param2 != NULL)
    std::cout << "Similarity factor: " << DiscriminantParameters::compareParameters(param1, param2) << std::endl << std::endl;
}

void writeDifferences(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  if (param1 != NULL && param2 != NULL)
  {
    DiscriminantParameters* difference = new DiscriminantParameters(param1->getVecDim());
    *difference = (*param1) + ((*param2)*(-1));

    std::ofstream diffFile(DIFF_FILE_NAME.c_str(), std::fstream::out | std::fstream::trunc);
    if (!diffFile.is_open())
    {
      std::cout << "Difference file \"" << DIFF_FILE_NAME << "\" could not be opened for writing!";
      exit(1);
    }

    unsigned vecDim = difference->getVecDim();
    diffFile << "LINEAR PARAMETERS:" << std::endl;
    if (linear)
    {
      for (unsigned i = 0; i < vecDim; i++)
      {
        diffFile << difference->get(i) << std::endl;
      }
    }

    unsigned matDim = difference->getMatDim();
    diffFile << "QUADRATIC PARAMETERS:" << std::endl;
    if (quadratic)
    {
      for (unsigned i = 0; i < matDim; i++)
      {
        for (unsigned j = 0; j <= i; j++)
        {
          diffFile << difference->get(i, j) << std::endl;
        }
      }
    }

    diffFile << "CONSTANT PARAMETER:" << std::endl;
    diffFile <<  difference->get() << std::endl;

    std::cout << "Differences written into \"" << DIFF_FILE_NAME << "\"..." << std::endl;
    diffFile.close();
    delete difference;
  }
}

std::string getFileName(std::string path)
{
  int slashPos = path.size() - 1;
  while (slashPos >= 0 && path[slashPos] != '/') slashPos--;
  int dotPos = path.size() - 1;
  while (dotPos > slashPos && path[dotPos] != '.') dotPos--;

  std::string protname;
  if (slashPos == dotPos)
    protname = path.substr(slashPos + 1, path.size() - (slashPos + 1));
  else
  {
    protname = path.substr(slashPos + 1, dotPos - (slashPos + 1));
  }

  std::ostringstream filename;
  filename << protname << OUTPUT_FILE_EXTENSION;
  return filename.str();
}

void writeValues(std::string filename1, DiscriminantParameters* param1, std::string filename2, DiscriminantParameters* param2)
{
  std::ofstream paramFile1;
  std::ofstream paramFile2;
  if (param1 != NULL)
  {
    std::string fName1 = getFileName(filename1);
    paramFile1.open(fName1.c_str(), std::fstream::out | std::fstream::trunc);
    if (!paramFile1.is_open())
    {
      std::cout << "Parameter file \"" << fName1 << "\" could not be opened for writing!";
      exit(1);
    }
  }
  if (param2 != NULL)
  {
    std::string fName2 = getFileName(filename2);
    paramFile2.open(fName2.c_str(), std::fstream::out | std::fstream::trunc);
    if (!paramFile2.is_open())
    {
      std::cout << "Parameter file \"" << fName2 << "\" could not be opened for writing!";
      exit(1);
    }
  }

  unsigned vecDim = 0;
  unsigned matDim = 0;
  if (param1 != NULL)
  {
    vecDim = param1->getVecDim();
    matDim = param1->getMatDim();
  }
  else if (param2 != NULL)
  {
    vecDim = param2->getVecDim();
    matDim = param2->getMatDim();
  }
  if (param1 != NULL) paramFile1 << "LINEAR PARAMETERS:" << std::endl;
  if (param2 != NULL) paramFile2 << "LINEAR PARAMETERS:" << std::endl;
  if (linear)
  {
    for (unsigned i = 0; i < vecDim; i++)
    {
      if (param1 != NULL) paramFile1 << param1->get(i) << std::endl;
      if (param2 != NULL) paramFile2 << param2->get(i) << std::endl;
    }
  }

  if (param1 != NULL) paramFile1 << "QUADRATIC PARAMETERS:" << std::endl;
  if (param2 != NULL) paramFile2 << "QUADRATIC PARAMETERS:" << std::endl;
  if (quadratic)
  {
    for (unsigned i = 0; i < matDim; i++)
    {
      for (unsigned j = 0; j <= i; j++)
      {
        if (param1 != NULL) paramFile1 << param1->get(i, j) << std::endl;
        if (param2 != NULL) paramFile2 << param2->get(i, j) << std::endl;
      }
    }
  }

  if (param1 != NULL)
  {
    paramFile1 << "CONSTANT PARAMETER:" << std::endl;
    paramFile1 <<  param1->get() << std::endl;
    paramFile1.close();
    std::cout << "Written parameter values from \"" << filename1 << "\" to file \"" << getFileName(filename1) << "\"..." << std::endl;
  }
  if (param2 != NULL)
  {
    paramFile2 << "CONSTANT PARAMETER:" << std::endl;
    paramFile2 <<  param2->get() << std::endl;
    paramFile2.close();
    std::cout << "Written parameter values from \"" << filename2 << "\" to file \"" << getFileName(filename2) << "\"..." << std::endl;
  }
}

int main(int argc, char** argv)
{
  progName = argv[0];

  DiscriminantParameters* param1 = NULL;
  DiscriminantParameters* param2 = NULL;
  std::string file1 = "";
  std::string file2 = "";

  if (argc > 5)
  {
    std::cout << "ERROR: Too many parameters!" << std::endl;
    printUsage();
    exit(1);
  }

  unsigned paramCount = 0;
  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] != '-')
    {
      if (paramCount == 0)
      {
        file1 = argv[i];
        paramCount++;
      }
      else if (paramCount == 1)
      {
        file2 = argv[i];
        paramCount++;
      }
    }
    else
    {
      if(argv[i][1] != '-')
        parseProcessingParameters(argv[i]);
      else
        parseDimensionParameters(argv[i]);
    }
  }

  param1 = loadParameters(file1);
  param2 = loadParameters(file2);

  if (param1 == NULL && param2 == NULL)
  {
    std::cout << "No parameters read!" << std::endl;
  }
  else
  {
    std::cout << std::endl;
    std::cout << "Number of parameters: ";
    if (param1 != NULL)
      std::cout << param1->getNumberOfParameters(linear, quadratic) << std::endl;
    else if (param2 != NULL)
      std::cout << param2->getNumberOfParameters(linear, quadratic) << std::endl;
    std::cout << std::endl;
    if (means)
      displayMeans(file1, param1, file2, param2);
    if (lengths)
      displayLengths(file1, param1, file2, param2);
    if (rmsd)
      displayRMSDs(file1, param1, file2, param2);
    if (distances)
      displayDistances(file1, param1, file2, param2);
    if (comp)
      displayComparison(file1, param1, file2, param2);
    if (differences)
      writeDifferences(file1, param1, file2, param2);
    if (values)
      writeValues(file1, param1, file2, param2);
  }

  delete param1;
  delete param2;
}
