#include "SVDSolver.h"
#include "SVDException.h"
#include "ParameterParser.h"
#include "ParserException.h"
#include "assert.h"
#include <fstream>
#include "math.h"
#include <sstream>
#include "time.h"

const std::string OUTPUT_FILE_NAME = "outputFile";
const std::string OUTPUT_FILE_DEFAULT = "";
const std::string SOLVING_METHOD_NAME = "method";
const std::string SOLVING_METHOD_DEFAULT = "";
const std::string THRESHOLD_NAME = "threshold";
const std::string THRESHOLD_DEFAULT = "";
const std::string DIMENSION_NAME = "dimension";
const std::string DIMENSION_DEFAULT = "";
const std::string WRITE_EIGENVALUES_NAME = "writeEigenvalues";
const std::string WRITE_EIGENVALUES_DEFAULT = "";
const std::string READ_EIGENVALUES_NAME = "readEigenvalues";
const std::string READ_EIGENVALUES_DEFAULT = "";
const std::string EIGENVALUE_OUTPUT_FILE_NAME = "eigenvalueFile";
const std::string EIGENVALUE_OUTPUT_FILE_DEFAULT = "";
const std::string WRITE_EIGENVECTORS_NAME = "writeEigenvectors";
const std::string WRITE_EIGENVECTORS_DEFAULT = "";
const std::string READ_EIGENVECTORS_NAME = "readEigenvectors";
const std::string READ_EIGENVECTORS_DEFAULT = "";
const std::string EIGENVECTOR_OUTPUT_FILE_NAME = "eigenvectorFile";
const std::string EIGENVECTOR_OUTPUT_FILE_DEFAULT = "";
const std::string WRITE_OFFDIAG_NAME = "writeOffdiagonal";
const std::string WRITE_OFFDIAG_DEFAULT = "";
const std::string READ_OFFDIAG_NAME = "readOffdiagonal";
const std::string READ_OFFDIAG_DEFAULT = "";
const std::string OFFDIAG_FILE_NAME = "offdiagonalFile";
const std::string OFFDIAG_FILE_DEFAULT = "";
const std::string WINDOW_NAME = "window";
const std::string WINDOW_DEFAULT = "";
const std::string KEYPOS_NAME = "keypos";
const std::string KEYPOS_DEFAULT = "";
const std::string LAMBDALIN_NAME = "lambdaLin";
const std::string LAMBDALIN_DEFAULT = "";
const std::string LAMBDAQUAD_NAME = "lambdaQuad";
const std::string LAMBDAQUAD_DEFAULT = "";
const std::string LAMBDA_DELTA_SQUARE_NAME = "lambdaDeltaSquare";
const std::string LAMBDA_DELTA_SQUARE_DEFAULT = "";
const std::string LAMBDA_QUAD_MAX_NAME = "lambdaQuadMax";
const std::string LAMBDA_QUAD_MAX_DEFAULT = "";
const std::string PATTERN_NAME = "posClass";
const std::string PATTERN_DEFAULT = "";
const std::string ANTIPATTERN_NAME = "negClass";
const std::string ANTIPATTERN_DEFAULT = "";

const std::string CONFIG_FILE_NAME = "solvePCA.cfg";

const std::string TAKE_BEST_STR = "take_best";
const std::string TAKE_ABOVE_THRESHOLD_STR = "take_above_threshold";
const std::string LEAVE_OUT_WORST_STR = "leave_out_worst";
const std::string LEAVE_OUT_BELOW_THRESHOLD_STR = "leave_out_below_threshold";
const std::string TAKE_PERCENTAGE_STR = "take_percentage";

const std::string WEIGHTS_KEYWORD = "#weight";
const std::string MATRIX_KEYWORD = "#matrix";
const std::string VECTOR_KEYWORD = "#vector";

const int PATTERN_NUM = 3;
const double TINY = 1e-20;
const unsigned NUM_OF_AMINO_ACIDS = 20;

/**
 * @brief A wrapper to hold some additional data for the linear equation system.
 */
typedef struct t_data_parameters
{
  /** @brief The dimension of the linear equation system meaning the number of compunents of the solution vector. */
  unsigned dim;
  /** @brief An array holding the matrix file names for the three secondary structure types. */
  std::string* matFileNames;
  /** @brief An array holding the vector file names for the three secondary structure types. */
  std::string* vecFileNames;
  /** @brief An array holding the weights for the three secondary structure types. */
  double* weights;
  /** @brief The window size for which the linear equation system should be solved. */
  unsigned window;
  /** @brief The key position for which the linear equation system should be solved. */
  unsigned keyPos;
  /** @brief The lambda factor for the linear terms. */
  double lambdaLin;
  /** @brief The lambda factor for the quadratic terms. */
  double lambdaQuad;
  /** @brief The width of the distribution of the quadratic lambda terms. */
  double lambdaDeltaSquare;
  /** @brief The maximum value of the distribution of the quadratic lambda terms. */
  double lambdaQuadMax;
  /** @brief The number of linear terms. */
  unsigned Ntot;
};

int getResiduePosition(uint residueNumber, uint index, uint length)
{
  uint jcomb = 0 ;
  uint profLength = NUM_OF_AMINO_ACIDS*length ;
  int position = -1 ;
  for(uint jm = 0; jm<profLength; jm++)
  {
    jcomb = jm ;
    for(uint jpmp = 0; jpmp<profLength; jpmp++)
    {
      jcomb += jpmp ;
      if(jpmp>=jm && jcomb==index)
      {
        switch(residueNumber)
        {
          case 0:
            position = (int)(jm/NUM_OF_AMINO_ACIDS);
            break;
          case 1:
            position = (int)(jpmp/NUM_OF_AMINO_ACIDS);
            break;
          default :
            std::cout << "\n>invalid residue number: " << residueNumber << std::endl ;
            exit(0) ;
        }
      }
    }
  }
  return position ;
}

double* readMatrix(std::string& filename, unsigned dim)
{
  // TODO Check correctness!
  std::ifstream inFile(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  double* mtx = new double[(dim * (dim + 1))/2];

  for(unsigned k = 0; k < dim; k++)
  {
    for(unsigned h = k; h < dim; h++)
    {
      double value;
      inFile.read((char*)&value, sizeof(value));
      mtx[dim * h  - ((h + 1) * h)/2 + k] = value;
    }
  }
  
  inFile.close();
  return mtx;
}

double* readVector(std::string& filename, unsigned dim)
{
  // TODO Check correctness!
  std::ifstream inFile(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  
  double* vec = new double[dim];
  
  for (unsigned i = 0; i < dim; i++)
  {
    double value;
    inFile.read((char*)&value, sizeof(value));
    vec[i] = value;
  }
  
  return vec;
}

double **setLambdaStrength(t_data_parameters* param, unsigned length, unsigned centre)
{
  double **vector = NULL ;
  vector = new double*[length] ;
  for(int k = 0; k < (int)length; k++)
  {
    vector[k] = NULL ;
    vector[k] = new double[length] ;
    for(int h = 0; h < (int)length; h++)
    {
      vector[k][h] = param->lambdaQuad * exp((double)(k-h)*(double)(k-h)/param->lambdaDeltaSquare) ;
      if(vector[k][h] > param->lambdaQuadMax) vector[k][h] = param->lambdaQuadMax ;
    }
  }

  return vector ;
}

void setupLES(t_data_parameters* param, int pattern, int antipattern, double *coefficients, double* righthand, double* averageX)
{
  uint i, j; // indexes for the two dimensions of the coefficient matrix
  uint ic = 0; // to compute matrix positions in the symmetric storage mode
  uint ij = 0; // matrix position ij in the symmetric storage mode
  double **Qstat; // samples statistic vectors
  double entry; // temporary entry from the matrices Nw
  double *norm = new double[PATTERN_NUM];
  double **actual_lambdaquad = NULL;
  unsigned int *N = new unsigned int[PATTERN_NUM];
  FILE **fp = new FILE*[PATTERN_NUM];
  actual_lambdaquad = setLambdaStrength(param, param->window, param->keyPos);
  Qstat = new double*[PATTERN_NUM];
  for(int p = 0; p < PATTERN_NUM; p++)
  {
    N[p] = 0 ;
    Qstat[p] = readVector(param->vecFileNames[p], param->dim);
    fp[p] = fopen(param->matFileNames[p].c_str(), "rb");
    if(fp[p] == NULL)
    {
      std::cout << "\n>error opening pattern " << p << " matrix file." << std::endl ;
      exit(1);
    }
    else
    {
      fseek(fp[p], -sizeof(unsigned int), SEEK_END);
      if(fread(&N[p], sizeof(unsigned int), 1, fp[p])!=1)
      {
        std::cout << "\nstats count failed." << std::endl;
        exit(1);
      }
      else rewind(fp[p]);
    }
    if(antipattern != -1 && p != antipattern && p != pattern) param->weights[p] = 0. ;
    norm[p] = param->weights[p]/((double)N[p] + TINY) ;
  }
  unsigned matDim = (param->dim * (param->dim + 1))/2;
  for (unsigned k = 0; k < matDim; k++) coefficients[k] = 0.0;
  for (j = 0; j < param->dim; j++)
  {
    for(int p = 0; p < PATTERN_NUM; p++) averageX[j] += norm[p]*Qstat[p][j] ;
    righthand[j] = param->weights[pattern]*(Qstat[pattern][j]/((double)N[pattern]+TINY) - averageX[j]);
  }
  for (j = 0; j < param->dim; j++)
  {
    ic += j ;
    ij = ic + j ;
    for (i = j; i < param->dim; i++)
    {
      for(int p = 0; p < PATTERN_NUM; p++)
      {
        if(fread(&entry, sizeof(double), 1, fp[p])!=1)
        {
          std::cout << "\nerror: couldn't read statistics matrix entries. Aborting...." << std::endl;
          std::cout << "j = " << j << "; i = " << i << "; p = " << p << std::endl;
          exit(1);
        }
        coefficients[ij] += norm[p]*(double)entry ;
      }
      coefficients[ij] -= averageX[i] * averageX[j] ;
      if(i == j)
      {
        double lambda = param->lambdaLin ;
        if(i >= param->Ntot)
        {
          int posA = getResiduePosition(0, i - param->Ntot, param->window) ;
          int posB = getResiduePosition(1, i - param->Ntot, param->window) ;
          if(posA < 0 || posB < 0)
          {
            std::cerr << "something's wrong with the indexing." << std::endl;
            exit(1);
          }
          else lambda = actual_lambdaquad[posA][posB] ;
        }
        coefficients[ij] += lambda ;
      }
      ij += i + 1 ;
    }
  }
  
  delete [] norm;
  delete [] N;
  for (int i = 0; i < PATTERN_NUM; i++)
  {
    delete [] Qstat[i];
    fclose(fp[i]);
  }
  delete [] fp;
  delete [] Qstat;
  for (unsigned i = 0; i < param->window; i++)
    delete [] actual_lambdaquad[i];
  delete [] actual_lambdaquad;
}

double calcConstTerm(double weight, unsigned dim, double* solution, double* xav)
{
  double sProd = 0;
  
  for (unsigned i = 0; i < dim; i++)
    sProd += solution[i] * xav[i];
  
  return (weight - sProd);
}

void writeSolution(SVDSolver* solver, std::string filename, double* xav, double weight)
{
  // TODO Check correctness!
  std::ofstream outFile(filename.c_str(), std::fstream::out | std::fstream::binary);
  if (!outFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for writing!" << std::endl;
    exit(1);
  }
  
  std::cout << "Writing solution to file \"" << filename << "\"..." << std::endl;
  
  double* sol = solver->getSolution();
  for (unsigned i = 0; i < solver->getDim(); i++)
    outFile.write(((char*)&sol[i]), sizeof(sol[i]));
  
  double constantTerm = calcConstTerm(weight, solver->getDim(), sol, xav);
  outFile.write(((char*)&constantTerm), sizeof(constantTerm));
  
  outFile.close();
}

void writeEigenvectors(SVDSolver* solver, std::string filename)
{
  std::ofstream outFile(filename.c_str(), std::fstream::out | std::fstream::binary);
  if (!outFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for writing!" << std::endl;
    exit(1);
  }
  
  std::cout << "Writing eigenvectors to file \"" << filename << "\"..." << std::endl;
  
  for (unsigned ev = 0; ev < solver->getDim(); ev++)
  {
    for (unsigned c = 0; c < solver->getDim(); c++)
    {
      double value = solver->getEigenvectorValue(c, ev, false);
      outFile.write(((char*)&value), sizeof(value));
    }
  }
  
  outFile.close();
}

void writeEigenvalues(SVDSolver* solver, std::string filename, std::string filename2)
{
  std::ofstream outFileBin(filename.c_str(), std::fstream::out | std::fstream::binary);
  if (!outFileBin.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for writing!" << std::endl;
    exit(1);
  }
  std::ofstream outFilePlain(filename2.c_str(), std::fstream::out | std::fstream::trunc);
  if (!outFilePlain.is_open())
  {
    std::cerr << "Could not open file \"" << filename2 << "\" for writing!" << std::endl;
    exit(1);
  }
  
  std::cout << "Writing eigenvalues to file \"" << filename << "\"..." << std::endl;
  
  for (unsigned ev = 0; ev < solver->getDim(); ev++)
  {
    double value = solver->getEigenvalue(ev, false);
    outFileBin.write(((char*)&value), sizeof(value));
  }
  
  std::cout << "Writing eigenvalues to file \"" << filename2 << "\"..." << std::endl;
  
  for (unsigned ev = 0; ev < solver->getDim(); ev++)
    outFilePlain << solver->getEigenvalue(ev, false) << std::endl; 
  
  outFileBin.close();
  outFilePlain.close();
}

void writeOffdiagonal(SVDSolver* solver, std::string filename)
{
  std::ofstream outFile(filename.c_str(), std::fstream::out | std::fstream::binary);
  if (!outFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for writing!" << std::endl;
    exit(1);
  }
  
  std::cout << "Writing offdiagonal elements to file \"" << filename << "\"..." << std::endl;
  
  for (unsigned i = 0; i < solver->getDim(); i++)
  {
    double value = solver->getOffDiagonalElement(i, false);
    outFile.write(((char*)&value), sizeof(value));
  }
  
  outFile.close();
}

double* readEigenvalues(SVDSolver* solver, std::string filename)
{
  std::ifstream inFile(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  
  std::cout << "Reading eigenvalues from file \"" << filename << "\"..." << std::endl;
  
  double* evVec = new double[solver->getDim()];
  for (unsigned ev = 0; ev < solver->getDim(); ev++)
  {
    double value;
    inFile.read(((char*)&value), sizeof(value));
    evVec[ev] = value;
  }
  
  inFile.close();
  
  return evVec;
}

double* readEigenvectors(SVDSolver* solver, std::string filename)
{
  std::ifstream inFile(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  
  std::cout << "Reading eigenvectors from file \"" << filename << "\"..." << std::endl;
  
  double* evVec = new double[(solver->getDim() * (solver->getDim() + 1))/2];
  for (unsigned ev = 0; ev < solver->getDim(); ev++)
  {
    for (unsigned c = 0; c < solver->getDim(); c++)
    {
      double value;
      inFile.read(((char*)&value), sizeof(value));
      evVec[c + ev * solver->getDim()] = value;
    }
  }
  
  inFile.close();
  
  return evVec;
}

double* readOffdiagonal(SVDSolver* solver, std::string filename)
{
  std::ifstream inFile(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::cerr << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  
  std::cout << "Reading offdiagonal elements from file \"" << filename << "\"..." << std::endl;
  
  double* od = new double[solver->getDim()];
  for (unsigned i = 0; i < solver->getDim(); i++)
  {
    double value;
    inFile.read(((char*)&value), sizeof(value));
    od[i] = value;
  }
  
  inFile.close();
  
  return od;
}

void initialiseParser(ParameterParser* parser)
{
  parser->addParameter(OUTPUT_FILE_NAME, OUTPUT_FILE_DEFAULT);
  parser->addParameter(SOLVING_METHOD_NAME, SOLVING_METHOD_DEFAULT);
  parser->addParameter(THRESHOLD_NAME, THRESHOLD_DEFAULT);
  parser->addParameter(DIMENSION_NAME, DIMENSION_DEFAULT);
  parser->addParameter(WRITE_EIGENVALUES_NAME, WRITE_EIGENVALUES_DEFAULT);
  parser->addParameter(READ_EIGENVALUES_NAME, READ_EIGENVALUES_DEFAULT);
  parser->addParameter(EIGENVALUE_OUTPUT_FILE_NAME, EIGENVALUE_OUTPUT_FILE_DEFAULT);
  parser->addParameter(WRITE_EIGENVECTORS_NAME, WRITE_EIGENVECTORS_DEFAULT);
  parser->addParameter(READ_EIGENVECTORS_NAME, READ_EIGENVECTORS_DEFAULT);
  parser->addParameter(EIGENVECTOR_OUTPUT_FILE_NAME, EIGENVECTOR_OUTPUT_FILE_DEFAULT);
  parser->addParameter(WRITE_OFFDIAG_NAME, WRITE_OFFDIAG_DEFAULT);
  parser->addParameter(READ_OFFDIAG_NAME, READ_OFFDIAG_DEFAULT);
  parser->addParameter(OFFDIAG_FILE_NAME, OFFDIAG_FILE_DEFAULT);
  parser->addParameter(WINDOW_NAME, WINDOW_DEFAULT);
  parser->addParameter(KEYPOS_NAME, KEYPOS_DEFAULT);
  parser->addParameter(LAMBDALIN_NAME, LAMBDALIN_DEFAULT);
  parser->addParameter(LAMBDAQUAD_NAME, LAMBDAQUAD_DEFAULT);
  parser->addParameter(LAMBDA_DELTA_SQUARE_NAME, LAMBDA_DELTA_SQUARE_DEFAULT);
  parser->addParameter(LAMBDA_QUAD_MAX_NAME, LAMBDA_QUAD_MAX_DEFAULT);
  parser->addParameter(PATTERN_NAME, PATTERN_DEFAULT);
  parser->addParameter(ANTIPATTERN_NAME, ANTIPATTERN_DEFAULT);
}

t_solving_type toSolvingType(std::string & str)
{
  std::string lowerStr = "";
  for (unsigned i = 0; i < str.length(); i++)
    lowerStr += std::tolower(str[i]);
  assert (lowerStr.length() == str.length());
  
  if (lowerStr == TAKE_BEST_STR)
    return TAKE_BEST;
  else if (lowerStr == TAKE_ABOVE_THRESHOLD_STR)
    return TAKE_ABOVE_THRESHOLD;
  else if (lowerStr == LEAVE_OUT_WORST_STR)
    return LEAVE_OUT_WORST;
  else if (lowerStr == LEAVE_OUT_BELOW_THRESHOLD_STR)
    return LEAVE_OUT_BELOW_THRESHOLD;
  else if (lowerStr == TAKE_PERCENTAGE_STR)
    return TAKE_PERCENTAGE;
  else
  {
    std::cerr << "Unknown solvin method string: " << str << std::endl;
    exit(1);
  }
}

void getWeights(t_data_parameters* param, std::string& line)
{
  unsigned i = 0;
  while (line[i] == ' ' && i < line.size()) i++; // skip leading whitespaces
  while (line[i] != ' ' && i < line.size()) i++; // skip keyword
  int str = 0;
  while (i < line.size())
  {
    while (line[i] == ' ' && i < line.size()) i++;
    unsigned start = i;
    while (line[i] != ' ' && i < line.size()) i++;
    unsigned length = i - start;
    if (length != 0 && i <= line.size())
    {
      param->weights[str] = atof(line.substr(start, length).c_str());
      str++;
    }
  }
  if (str < PATTERN_NUM)
  {
    std::cerr << PATTERN_NUM - str << " weights are missing!" << std::endl;
    exit(1);
  }
  
  double sum = 0;
  for (int i = 0; i < PATTERN_NUM; i++)
    sum += param->weights[i];
  for (int i = 0; i < PATTERN_NUM; i++)
    param->weights[i] /= sum;
}

void getMatrix(t_data_parameters* param, std::string& line)
{
  unsigned i = 0;
  while (line[i] == ' ' && i < line.size()) i++; // skip leading whitespaces
  while (line[i] != ' ' && i < line.size()) i++; // skip keyword
  while (line[i] == ' ' && i < line.size()) i++;
  unsigned start = i;
  while (line[i] != ' ' && i < line.size()) i++; // read structure index
  unsigned length = i - start;
  if (length == 0 || i == line.size())
  {
    std::cerr << "Missing pattern entry in matrix file description!" << std::endl;
    exit(1);
  }
  unsigned str = atoi(line.substr(start, length).c_str());
  if (param->matFileNames[str] != "")
  {
    std::cerr << "Matrix file name for pattern " << str << " has already been read!" << std::endl;
    exit(1);
  }
  str = atoi(line.substr(start, length).c_str());
  while (line[i] == ' ' && i < line.size()) i++;
  start = i;
  while (line[i] != ' ' && i < line.size()) i++; // read file name
  length = i - start;
  if (length == 0)
  {
    std::cerr << "Missing file name in matrix file description!" << std::endl;
    exit(1);
  }
  param->matFileNames[str] = line.substr(start, length);
}

void getVector(t_data_parameters* param, std::string& line)
{
  unsigned i = 0;
  while (line[i] == ' ' && i < line.size()) i++; // skip leading whitespaces
  while (line[i] != ' ' && i < line.size()) i++; // skip keyword
  while (line[i] == ' ' && i < line.size()) i++;
  unsigned start = i;
  while (line[i] != ' ' && i < line.size()) i++; // read structure index
  unsigned length = i - start;
  if (length == 0 || i == line.size())
  {
    std::cerr << "Missing pattern entry in matrix file description!" << std::endl;
    exit(1);
  }
  unsigned str = atoi(line.substr(start, length).c_str());
  if (param->vecFileNames[str] != "")
  {
    std::cerr << "Vector file name for pattern " << str << " has already been read!" << std::endl;
    exit(1);
  }
  str = atoi(line.substr(start, length).c_str());
  while (line[i] == ' ' && i < line.size()) i++;
  start = i;
  while (line[i] != ' ' && i < line.size()) i++; // read file name
  length = i - start;
  if (length == 0)
  {
    std::cerr << "Missing file name in vector file description!" << std::endl;
    exit(1);
  }
  param->vecFileNames[str] = line.substr(start, length);
}

void readParameters(t_data_parameters* param, const std::string& filename)
{
  param->matFileNames = new std::string[PATTERN_NUM];
  param->vecFileNames = new std::string[PATTERN_NUM];
  param->weights = new double[PATTERN_NUM];
  for (int i = 0; i < PATTERN_NUM; i++)
  {
    param->matFileNames[i] = "";
    param->vecFileNames[i] = "";
    param->weights[i] = 0;
  }
  
  std::ifstream cfgFile(filename.c_str(), std::fstream::in);
  std::string line;
  while (getline(cfgFile, line, '\n'))
  {
    unsigned i = 0;
    while (line[i] == ' ' && i < line.size()) i++;
    std::string keyword = line.substr(i, WEIGHTS_KEYWORD.size());
    if (keyword == WEIGHTS_KEYWORD) getWeights(param, line);
    if (keyword == MATRIX_KEYWORD) getMatrix(param, line);
    if (keyword == VECTOR_KEYWORD) getVector(param, line);
  }
  
  bool errors = false;
  for (int i = 0; i < PATTERN_NUM; i++)
  {
    if (param->matFileNames[i] == "")
    {
      std::cerr << "No matrix file for pattern " << i << " specified!" << std::endl;
      errors = true;
    }
    if (param->vecFileNames[i] == "")
    {
      std::cerr << "No vector file for pattern " << i << " specified!" << std::endl;
      errors = true;
    }
  }  
}

void deleteParameters(t_data_parameters* param)
{
  delete [] param->matFileNames;
  delete [] param->vecFileNames;
  delete [] param->weights;
  delete param;
}

void writeDataParameters(t_data_parameters* param)
{
  std::cout << "dim = " << param->dim << std::endl;
  std::cout << "window = " << param->window << std::endl;
  std::cout << "keyPos = " << param->keyPos << std::endl;
  std::cout << "lambdaLin = " << param->lambdaLin << std::endl;
  std::cout << "lambdaQuad = " << param->lambdaQuad << std::endl;
  std::cout << "lambdaDeltaSquare = " << param->lambdaDeltaSquare << std::endl;
  std::cout << "lambdaQuadMax = " << param->lambdaQuadMax << std::endl;
  std::cout << "Ntot = " << param->Ntot << std::endl;
  std::cout << "matFileNames:" << std::endl;
  for (int i = 0; i < PATTERN_NUM; i++)
    std::cout << '\t' << i << ":\t" << param->matFileNames[i] << std::endl;
  std::cout << "vecFileNames:" << std::endl;
  for (int i = 0; i < PATTERN_NUM; i++)
    std::cout << '\t' << i << ":\t" << param->vecFileNames[i] << std::endl;
  std::cout << "weights:" << std::endl;
  for (int i = 0; i < PATTERN_NUM; i++)
    std::cout << '\t' << i << ":\t" << param->weights[i] << std::endl;
}

std::string timeString(unsigned totalSeconds)
{
  unsigned seconds = totalSeconds % 60;
  unsigned temp = totalSeconds / 60; // time in minutes
  unsigned minutes = temp % 60;
  temp = temp / 60; // time in hours
  unsigned hours = temp % 24;
  unsigned days = temp / 24;
  std::ostringstream timeStr;
  if(days < 10) timeStr << "0";
  timeStr << days << ":";
  if(hours < 10) timeStr << "0";
  timeStr << hours << ":";
  if(minutes < 10) timeStr << "0";
  timeStr << minutes << ":";
  if(seconds < 10) timeStr << "0";
  timeStr << seconds;
  return timeStr.str();
}

int main(unsigned argc, char** argv)
{
  try
  {
    ParameterParser* parser = new ParameterParser();
  
    initialiseParser(parser);
  
    std::ifstream cfgFile(CONFIG_FILE_NAME.c_str(), std::fstream::in);
    if (!cfgFile.is_open())
      parser->writeParameterFile(CONFIG_FILE_NAME);
    else
      cfgFile.close();
  
    parser->readParameters(CONFIG_FILE_NAME);
  
    unsigned dimension = atoi((parser->getParameterValue(DIMENSION_NAME)).c_str());
    std::string solvingMethodStr = parser->getParameterValue(SOLVING_METHOD_NAME);
    t_solving_type method = toSolvingType(solvingMethodStr);
    double threshold = atof((parser->getParameterValue(THRESHOLD_NAME)).c_str());
    bool writeEigenvals = (strcmp((parser->getParameterValue(WRITE_EIGENVALUES_NAME)).c_str(), "1") == 0);
    bool readEigenvals = (strcmp((parser->getParameterValue(READ_EIGENVALUES_NAME)).c_str(), "1") == 0);
    std::string eigenvalueFile = parser->getParameterValue(EIGENVALUE_OUTPUT_FILE_NAME);
    bool writeEigenvecs = (strcmp((parser->getParameterValue(WRITE_EIGENVECTORS_NAME)).c_str(), "1") == 0);
    bool readEigenvecs = (strcmp((parser->getParameterValue(READ_EIGENVECTORS_NAME)).c_str(), "1") == 0);
    std::string eigenvectorFile = parser->getParameterValue(EIGENVECTOR_OUTPUT_FILE_NAME);
    bool writeOffdiag = (strcmp((parser->getParameterValue(WRITE_OFFDIAG_NAME)).c_str(), "1") == 0);
    bool readOffdiag = (strcmp((parser->getParameterValue(READ_OFFDIAG_NAME)).c_str(), "1") == 0);
    std::string offdiagFile = parser->getParameterValue(OFFDIAG_FILE_NAME);
    std::string outputFile = parser->getParameterValue(OUTPUT_FILE_NAME);
    unsigned window = atoi((parser->getParameterValue(WINDOW_NAME)).c_str());
    unsigned keyPos = atoi((parser->getParameterValue(KEYPOS_NAME)).c_str());
    double lambdaLin = atof((parser->getParameterValue(LAMBDALIN_NAME)).c_str());
    double lambdaQuad = atof((parser->getParameterValue(LAMBDAQUAD_NAME)).c_str());
    double lambdaDeltaSquare = atof((parser->getParameterValue(LAMBDA_DELTA_SQUARE_NAME)).c_str());
    double lambdaQuadMax = atof((parser->getParameterValue(LAMBDA_QUAD_MAX_NAME)).c_str());
    int pattern = atoi((parser->getParameterValue(PATTERN_NAME)).c_str());
    int antipattern = atoi((parser->getParameterValue(ANTIPATTERN_NAME)).c_str());
    
    double* matrix = new double[(dimension * (dimension + 1))/2];
    double* vector = new double[dimension];
    double* xav = new double[dimension];
    
    t_data_parameters* param = new t_data_parameters;
    param->dim = dimension;
    param->window = window;
    param->keyPos = keyPos;
    param->lambdaLin = lambdaLin;
    param->lambdaQuad = lambdaQuad;
    param->lambdaDeltaSquare = lambdaDeltaSquare;
    param->lambdaQuadMax = lambdaQuadMax;
    param->Ntot = NUM_OF_AMINO_ACIDS * window;
    readParameters(param, CONFIG_FILE_NAME);
    
    writeDataParameters(param);

    time_t t0 = time(NULL);
    std::cout << "Setting up linear equation system..." << std::endl;
    setupLES(param, pattern, antipattern, matrix, vector, xav);
    time_t t1 = time(NULL);
    std::cout << "Initialising SVDSolver..." << std::endl;
    SVDSolver* solver = new SVDSolver(matrix, vector, dimension);
    if (readEigenvals) solver->setEigenvalues(readEigenvalues(solver, eigenvalueFile));
    if (readEigenvecs) solver->setEigenvectors(readEigenvectors(solver, eigenvectorFile));
    if (readOffdiag) solver->setOffDiagonal(readOffdiagonal(solver, offdiagFile));
    time_t t2 = time(NULL);
    
    std::cout << "Solving linear equation system..." << std::endl;
    solver->solve(method, threshold);
    time_t t3 = time(NULL);
    
    std::cout << "Writing data..." << std::endl;
    writeSolution(solver, outputFile, xav, param->weights[pattern]);
    if (writeEigenvecs) writeEigenvectors(solver, eigenvectorFile);
    if (writeEigenvals)
    { 
      std::ostringstream plainFileName;
      plainFileName << eigenvalueFile << ".plain";
      writeEigenvalues(solver, eigenvalueFile, plainFileName.str().c_str());
    }
    if (writeOffdiag) writeOffdiagonal(solver, offdiagFile);
    time_t t4 = time(NULL);
    
    std::cout << "Set up time: " << timeString(t1 - t0) << std::endl;
    std::cout << "Initialisation time: " << timeString(t2 - t1) << std::endl;
    std::cout << "Solving time: " << timeString(t3 - t2) << std::endl;
    std::cout << "Writing time: " << timeString(t4 - t3) << std::endl;
    
    delete parser;
    delete solver;
    delete [] xav;
    deleteParameters(param);
  }
  catch (ParserException e)
  {
    std::cerr << "ParserException: " << e.what() << std::endl;
    exit(1);
  }
  catch (SVDException e)
  {
    std::cerr << "SVDException: " << e.what() << std::endl;
    exit(1);
  }
}
