#include "ANN.h"
#include "ANNException.h"
#include <sstream>
#include <list>
#include "math.h"
//#include "assert.h"
#include <stdlib.h>

#ifdef ANN_DEBUG_MODE
const std::string DEBUG_FILE_NAME = "ann.debug";
const std::string DEBUG_DATA_FILE_NAME = "ann.data";
const std::string ANN::WRITE_PREFIX_STRING = "annlearn";
const std::string ANN::WRITE_SUFFIX_STRING = ".prm";
#endif

/**
 * @brief Standard constructor for the <code>ANN</code> class.
 *
 * @param inputNodes The number of nodes in the input layer
 * @param hiddenNodes The number of nodes in the hidden layer
 * @param outputNodes The number of nodes in the output layer
 * @param seed The random seed for the random number generator used to randomly initialize the parameters of the ANN
 * @param rndBound The maximum absolute value for the randomly initialized parameters (parameters vary between <code>-rndBound</code> and <code>rndBound</code>)
 */
ANN::ANN(unsigned inputNodes, unsigned hiddenNodes, unsigned outputNodes, unsigned seed, double rndBound)
{
  this->inputNodes = inputNodes;
  this->hiddenNodes = hiddenNodes;
  this->outputNodes = outputNodes;

#ifdef ANN_DEBUG_MODE
  this->debugFile.open(DEBUG_FILE_NAME.c_str(), std::fstream::out | std::fstream::trunc);
  if (!this->debugFile.is_open())
  {
    std::ostringstream errMsg;
    errMsg << "Could not open file \"" << DEBUG_FILE_NAME << "\" for writing!";
    throw ANNException(errMsg.str().c_str());
  }
  this->debugFile << "Creating new ANN(" << this->inputNodes << ", " << this->hiddenNodes << ", " << this->outputNodes << ", " << seed << ", " << rndBound << ")" << std::endl;

  this->dataFile.open(DEBUG_DATA_FILE_NAME.c_str(), std::fstream::out | std::fstream::trunc);
  if (!this->dataFile.is_open())
  {
    std::ostringstream errMsg;
    errMsg << "Could not open file \"" << DEBUG_DATA_FILE_NAME << "\" for writing!";
    throw ANNException(errMsg.str().c_str());
  }
  this->dataFile << "iterations\terror\tlast error value\terror value change" << std::endl;
  this->period = 0;
#endif

  this->inputActivation = new double[this->inputNodes + 1];
  this->hiddenActivation = new double[this->hiddenNodes + 1];
  this->outputActivation = new double[this->outputNodes];

  this->hiddenParameters = new double*[this->hiddenNodes];
  for (unsigned i = 0; i < this->hiddenNodes; i++)
    this->hiddenParameters[i] = new double[this->inputNodes + 1];

  this->outputParameters = new double*[this->outputNodes];
  for (unsigned i = 0; i < this->outputNodes; i++)
    this->outputParameters[i] = new double[this->hiddenNodes + 1];

  for (unsigned i = 0; i < this->inputNodes; i++)
    this->inputActivation[i] = 0;
  this->inputActivation[this->inputNodes] = 1.0;
  this->hiddenActivation[this->hiddenNodes] = 1.0;

  if (seed == 0)
    srand48(time(NULL));
  else
    srand48(seed);
  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    this->hiddenActivation[i] = 0;
    for (unsigned j = 0; j <= this->inputNodes; j++)
    {
      this->hiddenParameters[i][j] = 2.0 * rndBound * drand48() - rndBound;
      while (this->hiddenParameters[i][j] == 0) this->hiddenParameters[i][j] = 2.0 * rndBound * drand48() - rndBound;
      assert(this->hiddenParameters[i][j] != 0);
    }
  }

  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    this->outputActivation[i] = 0;
    for (unsigned j = 0; j <= this->hiddenNodes; j++)
    {
      this->outputParameters[i][j] = 2.0 * rndBound * drand48() - rndBound;
      while (this->outputParameters[i][j] == 0) this->outputParameters[i][j] = 2.0 * rndBound * drand48() - rndBound;
      assert(this->outputParameters != 0);
    }
  }

  this->funcType = ANN::DEFAULT_FUNCTION_TYPE;
  this->learningConstant = ANN::DEFAULT_LEARNING_CONSTANT;
  this->momentumConstant = ANN::DEFAULT_MOMENTUM_CONSTANT;

  this->learnData = NULL;
  this->expectedOutput = NULL;
  this->samples = 0;
  this->learnDataCopied = false;

  this->validating = false;
  this->validationIt = 0;
  this->validationData = NULL;
  this->validationExpected = NULL;
  this->validationSamples = 0;
  this->validationCopied = false;
  this->validationMaxIt = ANN::DEFAULT_VALIDATION_MAX_IT;
  this->bestHiddenParameters = NULL;
  this->bestOutputParameters = NULL;
  this->bestErrorValue = 1e100;

  this->maxIt = ANN::DEFAULT_MAX_IT;
  this->minError = ANN::DEFAULT_MIN_ERROR;

  this->a = ANN::DEFAULT_A;
  this->b = ANN::DEFAULT_B;

  this->shuffleSet = ANN::DEFAULT_SHUFFLE_SET;
  this->adjustFrequency = ANN::DEFAULT_ADJUST_FREQUENCY;
  this->adjustConstants = ANN::DEFAULT_ADJUST_CONSTANTS;
  this->adjustFactor = ANN::DEFAULT_ADJUST_FACTOR;
  this->adjustMin = ANN::DEFAULT_ADJUST_MIN;
  this->adjustMax = ANN::DEFAULT_ADJUST_MAX;
  this->rhoPlus = ANN::DEFAULT_RHO_PLUS;
  this->rhoMinus = ANN::DEFAULT_RHO_MINUS;

  this->w_zero = ANN::DEFAULT_W_ZERO;
  this->regularization = ANN::DEFAULT_REGULARIZATION;
  this->lambda = ANN::DEFAULT_LAMBDA;
  this->power = ANN::DEFAULT_POWER;
}

/**
 * @brief Standard destructor for the <code>ANN</code> class.
 */
ANN::~ANN()
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Deleting ANN..." << std::endl;
#endif
  delete [] this->inputActivation;
  delete [] this->hiddenActivation;
  delete [] this->outputActivation;

  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    delete [] this->hiddenParameters[i];
    if (this->bestHiddenParameters != NULL) delete [] this->bestHiddenParameters[i];
  }
  delete [] this->hiddenParameters;
  if (this->bestHiddenParameters != NULL) delete [] this->bestHiddenParameters;

  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    delete [] this->outputParameters[i];
    if (this->bestOutputParameters != NULL) delete [] this->bestOutputParameters[i];
  }
  delete [] this->outputParameters;
  if (this->bestOutputParameters != NULL) delete [] this->bestOutputParameters;

  if (this->learnDataCopied && this->learnData != NULL) this->freeDataset(this->learnData, this->expectedOutput, this->samples);
  if (this->validationCopied && this->learnData != NULL) this->freeDataset(this->validationData, this->validationExpected, this->validationSamples);

#ifdef ANN_DEBUG_MODE
  this->debugFile.close();
  this->dataFile.close();
#endif
}

/**
 * @brief Function to initialize the learning dataset for this ANN.
 *
 * @param learnData The input values of the learning dataset
 * @param expectedOutput The expected outputs for the learning dataset
 * @param samples The number of samples within the dataset
 * @param makeCopy Determines whether to make a copy of the dataset or use the given pointers directly
 */
void ANN::initLearnData(double** learnData, double** expectedOutput, unsigned samples, bool makeCopy)
{
  assert(learnData != NULL);
  assert(expectedOutput != NULL);
  assert(samples != 0);

#ifdef ANN_DEBUG_MODE
  this->debugFile << "Initializing learn dataset (samples = " << samples << ", makeCopy = " << makeCopy << ")" << std::endl;
#endif

  if (this->learnDataCopied && this->learnData != NULL) this->freeDataset(this->learnData, this->expectedOutput, this->samples);
  this->samples = samples;
  this->learnDataCopied = makeCopy;

  if (makeCopy)
  {
    this->learnData = new double*[this->samples];
    this->expectedOutput = new double*[this->samples];
    for (unsigned i = 0; i < this->samples; i++)
    {
      this->expectedOutput[i] = new double[this->outputNodes];
      this->learnData[i] = new double[this->inputNodes];
      for (unsigned j = 0; j < this->inputNodes; j++)
        this->learnData[i][j] = learnData[i][j];
      for (unsigned j = 0; j < this->outputNodes; j++)
        this->expectedOutput[i][j] = expectedOutput[i][j];
    }
  }
  else
  {
    this->learnData = learnData;
    this->expectedOutput = expectedOutput;
  }
}

/**
 * @brief Function to initialize the validation dataset for this ANN.
 *
 * @param validationData The input values of the validation dataset
 * @param validationExpected The expected outputs for the validation dataset
 * @param validationSamples The number of samples within the dataset
 * @param makeCopy Determines whether to make a copy of the dataset or use the given pointers directly
 */
void ANN::initValidationData(double** validationData, double** validationExpected, unsigned validationSamples, bool makeCopy)
{
  assert(validationData != NULL);
  assert(validationExpected != NULL);
  assert(validationSamples != 0);

#ifdef ANN_DEBUG_MODE
  this->debugFile << "Initializing validation dataset (samples = " << validationSamples << ", makeCopy = " << makeCopy << ")" << std::endl;
#endif

  if (this->validationCopied && this->validationData != NULL) this->freeDataset(this->validationData, this->validationExpected, this->validationSamples);
  this->validationSamples = validationSamples;
  this->validationCopied = makeCopy;

  if (makeCopy)
  {
    this->validationData = new double*[this->validationSamples];
    this->validationExpected = new double*[this->validationSamples];
    for (unsigned i = 0; i < this->validationSamples; i++)
    {
      this->validationExpected[i] = new double[this->outputNodes];
      this->validationData[i] = new double[this->inputNodes];
      for (unsigned j = 0; j < this->inputNodes; j++)
        this->validationData[i][j] = validationData[i][j];
      for (unsigned j = 0; j < this->outputNodes; j++)
        this->validationExpected[i][j] = validationExpected[i][j];
    }
  }
  else
  {
    this->validationData = validationData;
    this->validationExpected = validationExpected;
  }

  this->validating = true;

  if (this->bestHiddenParameters == NULL)
  {
    this->bestHiddenParameters = new double*[this->hiddenNodes];
    for (unsigned i = 0; i < this->hiddenNodes; i++)
    {
      this->bestHiddenParameters[i] = new double[this->inputNodes + 1];
      for (unsigned j = 0; j <= this->inputNodes; j++)
          this->bestHiddenParameters[i][j] = 0;
    }
  }

  if (this->bestOutputParameters == NULL)
  {
    this->bestOutputParameters = new double*[this->outputNodes];
    for (unsigned i = 0; i < this->outputNodes; i++)
    {
      this->bestOutputParameters[i] = new double[this->hiddenNodes + 1];
      for (unsigned j = 0; j <= this->hiddenNodes; j++)
      {
        this->bestOutputParameters[i][j] = 0;
      }
    }
  }
}

/**
 * The learning algorithm is the simplest version of the on-line backpropagation algorithm. It needs the learning dataset (<code>learnData</code> and 
 * <code>expectedOutput</code>) to be initialized.
 * @brief Function running the backpropagation learning algorithm.
 *
 * @param maxIterations The maximum number of iterations for the learning algorithm (if 0 the algorithm is not bounded by any number of iterations)
 * @param minError The lower bound of the change in error causing the algorithm to stop if the current error gets below this value
 * @param seed Seed used when randomly shuffling the data set (default value 0 means that no see is used)
 */
void ANN::learn(unsigned maxIterations, double minError, unsigned seed)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Starting learning procedure (maxIt = " << maxIterations << ", minError = " << minError << ")" << std::endl;
  this->debugFile << this->getStats() << std::endl;
#endif
  this->maxIt = maxIterations;
  this->minError = minError;

  // initialize arrays for the deltas necessary for the backpropagation algorithm
  double* hiddenDeltas = new double[this->hiddenNodes];
  double* outputDeltas = new double[this->outputNodes];

  // initialize structure to store the previous parameter increments (necessary for the momentum term)
  double** previousHiddenParameters = NULL;
  double** previousOutputParameters = NULL;
  if (this->momentumConstant != 0)
  {
    previousHiddenParameters = new double*[this->hiddenNodes];
    for (unsigned i = 0; i < this->hiddenNodes; i++)
    {
      previousHiddenParameters[i] = new double[this->inputNodes + 1];
      for (unsigned j = 0; j <= this->inputNodes; j++)
        previousHiddenParameters[i][j] = 0;
    }

    previousOutputParameters = new double*[this->outputNodes];
    for (unsigned i = 0; i < this->outputNodes; i++)
    {
      previousOutputParameters[i] = new double[this->hiddenNodes + 1];
      for (unsigned j = 0; j <= this->hiddenNodes; j++)
        previousOutputParameters[i][j] = 0;
    }
  }

  bool done = false;
  double lastError = 0;
  unsigned iterations = 0;
  unsigned* setIndices = NULL;
  unsigned* list1 = NULL;
  unsigned* list2 = NULL;
  unsigned list1Index = 0;
  unsigned list2Index = 0;
  if (this->shuffleSet)
  {
    setIndices = new unsigned[this->samples];
    list1 = new unsigned[this->samples];
    list2 = new unsigned[this->samples];
    for (unsigned i = 0; i < this->samples; i++) setIndices[i] = i;
    if (seed == 0)
      srand48(time(NULL));
    else
      srand48(seed);
  }
  while (!done)
  {
#ifdef ANN_DEBUG_MODE
    this->debugFile << "ITERATION " << iterations << std::endl;
#endif
    if (this->shuffleSet)
    {
      list1Index = 0;
      list2Index = 0;
    }
    //begin epoch
    for (unsigned sCount = 0; sCount < this->samples; sCount++)
    {
      unsigned s = sCount;
      if (this->shuffleSet)
        s = setIndices[sCount];
      // forward propagation
      this->updateActivation(this->learnData[s]);

      // backward propagation
      // first compute output deltas
      for (unsigned i = 0; i < this->outputNodes; i++)
      {
        double d = this->expectedOutput[s][i]; // The desired output of the i-th output node
        double o = this->outputActivation[i]; // The output of the i-th output node
        // calculation follows the formula delta[i] = error[i] * activationFunction'(sum[i]), where sum[i] = sum(w[i][j] * x[j]) (x is output of hidden layer, error[i] = d -o)
        switch (this->funcType)
        {
          case ANN_LOGISTIC:
            //activationFunction'(sum[i]) = a * y[i] * (1 - y[i]) where y[i] = output of i-th output neuron
            outputDeltas[i] = this->a * (d - o) * o * (1.0 - o);
            break;
          case ANN_TANH:
            //activationFunction'(sum[i]) = b/a * (a - y[i]) * (a + y[i]) where y[i] = output of i-th output neuron
            outputDeltas[i] = this->b/this->a * (d - o) * (this->a - o) * (this->a + o);
            break;
          default:
            assert(false);
        }
      }

      // then compute hidden deltas
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        // compute the backpropagated local gradient g[i] = sum{k}(outputDeltas[k] * outputParameters[k][i])
        double g = 0;
        for (unsigned k = 0; k < this->outputNodes; k++)
          g += outputDeltas[k] * this->outputParameters[k][i];

        double y = this->hiddenActivation[i]; // output of the i-th hidden neuron

        // calculation follow the formula delta[i] = g * activationFunction'(sum[i]), where sum[i] = sum(w[i][j] * x[j]) (x is output of input layer)
        switch (this->funcType)
        {
          case ANN_LOGISTIC:
            //activationFunction'(sum[i]) = a * y[i] * (1 - y[i]) where y[i] = output of i-th output neuron
            hiddenDeltas[i] = this->a * y * (1 - y) * g;
            break;
          case ANN_TANH:
            //activationFunction'(sum[i]) = b/a * (a - y[i]) * (a + y[i]) where y[i] = output of i-th output neuron
            hiddenDeltas[i] = this->b/this->a * (a - y) * (a + y) * g;
            break;
          default:
            assert(false);
        }
      }

      // adjust synaptic weights (parameters) for the output layer
      for (unsigned i = 0; i < this->outputNodes; i++)
      {
        double factor = this->learningConstant * outputDeltas[i];
        for (unsigned j = 0; j <= this->hiddenNodes; j++)
        {
          double deltaW = factor * this->hiddenActivation[j];

          // add the derivative of the complexity term
          if (this->lambda != 0)
            deltaW += this->learningConstant * this->getPenalityDerivative(true, i, j);

          // add the momentum term
          if (this->momentumConstant != 0)
          {
            deltaW += this->momentumConstant * previousOutputParameters[i][j];
            previousOutputParameters[i][j] = deltaW;
          }
          this->outputParameters[i][j] += deltaW;
        }
      }

      // adjust synaptic weights (parameters) for the hidden layer
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        double factor = this->learningConstant * hiddenDeltas[i];
        for (unsigned j = 0; j <= this->inputNodes; j++)
        {
          double deltaW = factor * this->inputActivation[j];

          // add the derivative of the complexity term
          if (this->lambda != 0)
            deltaW += this->learningConstant * this->getPenalityDerivative(false, i, j);

          // add the momentum term
          if (this->momentumConstant != 0)
          {
            deltaW += this->momentumConstant * previousHiddenParameters[i][j];
            previousHiddenParameters[i][j] = deltaW;
          }
          this->hiddenParameters[i][j] += deltaW;
        }
      }

      if (this->shuffleSet)
      {
        double rnd = drand48();
        if (rnd > 0.5)
          list1[list1Index++] = s;
        else
          list2[list2Index++] = s;
      }
    }

    // check if already done
#ifdef ANN_DEBUG_MODE
    double error = this->getError();
    double errorChange = lastError - error;
    if (errorChange < 0)
      errorChange = -errorChange;
#ifdef ANN_PROCENTUAL_ERROR
    if (lastError != 0) errorChange /= lastError;
#endif
#endif
    if (!this->validating)
    {
#ifndef ANN_DEBUG_MODE
      double error = 0;
#endif
      if (this->minError != 0.0)
      {
#ifndef ANN_DEBUG_MODE
        error = this->getError();
        double errorChange = lastError - error;
        if (errorChange < 0)
          errorChange = -errorChange;
#ifdef ANN_PROCENTUAL_ERROR
        if (lastError != 0) errorChange /= lastError;
#endif
        lastError = error;
#endif
        done = (errorChange <= this->minError);
#ifdef ANN_DEBUG_MODE
        this->debugFile << "done = (errorChange <= this->minError) = " << done << std::endl;
#endif
      }
      if (this->maxIt != 0)
      {
        if (this->minError != 0.0)
        {
          done = (done || (iterations >= this->maxIt - 1));
#ifdef ANN_DEBUG_MODE
          this->debugFile << "done = (done || (iterations >= this->maxIt - 1)) = " << done << std::endl;
#endif
        }
        else
        {
          done = (iterations >= this->maxIt - 1);
#ifdef ANN_DEBUG_MODE
          this->debugFile << "done = (iterations >= this->maxIt - 1) = " << done << std::endl;
#endif
        }
      }
    }
    else
    {
      this->updateValidation();
      if (this->maxIt != 0)
        done = (iterations >= this->maxIt);
      done = done || (this->validationIt >= this->validationMaxIt);
    }

#ifdef ANN_DEBUG_MODE
    this->debugFile << iterations << ": error = " << error << "; lastError = " << lastError << "; change = " << errorChange << "; minError = " << this->minError << std::endl;
    this->dataFile << iterations << "\t" << error << "\t" << lastError << "\t" << errorChange << std::endl;
    if (this->period != 0 && iterations % this->period == 0)
    {
      std::ostringstream filename;
      filename << WRITE_PREFIX_STRING << iterations << WRITE_SUFFIX_STRING;
      std::string filenameStr = filename.str();
      this->writeParameters(filenameStr);
    }
    lastError = error;
#endif

    // randomize order of training examples by merging list1 and list2
    if (this->shuffleSet)
    {
      assert(list1Index + list2Index == this->samples);
      unsigned counter = 0;
      unsigned i1 = 0;
      unsigned i2 = 0;
      while (i1 < list1Index || i2 < list2Index)
      {
        if (i1 < list1Index)
          setIndices[counter++] = list1[i1++];
        if (i2 < list2Index)
          setIndices[counter++] = list2[i2++];
      }
    }

    if (this->adjustConstants) this->adjustLearnConstant(iterations);
    iterations++;
    // TODO optimal brain surgeon?
  }

#ifdef ANN_DEBUG_MODE
  if (this->period != 0)
  {
    std::ostringstream filename;
    filename << WRITE_PREFIX_STRING << iterations << WRITE_SUFFIX_STRING;
    std::string filenameStr = filename.str();
    this->writeParameters(filenameStr);
  }
#endif

  if (this->validating)
  {
    for (unsigned i = 0; i < this->hiddenNodes; i++)
    {
      for (unsigned j = 0; j <= this->inputNodes; j++)
      {
        this->hiddenParameters[i][j] = this->bestHiddenParameters[i][j];
      }
    }
    for (unsigned i = 0; i < this->outputNodes; i++)
    {
      for (unsigned j = 0; j <= this->hiddenNodes; j++)
      {
        this->outputParameters[i][j] = this->bestOutputParameters[i][j];
      }
    }
  }

  delete [] hiddenDeltas;
  delete [] outputDeltas;
  if (this->shuffleSet)
  {
    delete [] setIndices;
    delete [] list1;
    delete [] list2;
  }

  if (this->momentumConstant != 0)
  {
    for (unsigned i = 0; i < this->hiddenNodes; i++)
      delete [] previousHiddenParameters[i];
    delete [] previousHiddenParameters;

    for (unsigned i = 0; i < this->outputNodes; i++)
      delete [] previousOutputParameters[i];
    delete [] previousOutputParameters;
  }
}

/**
 * The learning algorithm is a version of the on-line backpropagation algorithm which uses different learning constants for different weights. It needs 
 * the learning dataset (<code>learnData</code> and <code>expectedOutput</code>) to be initialized.
 * @brief Function running the backpropagation learning algorithm.
 *
 * @param maxIterations The maximum number of iterations for the learning algorithm (if 0 the algorithm is not bounded by any number of iterations)
 * @param minError The lower bound of the change in error causing the algorithm to stop if the current error gets below this value
 * @param seed Seed used when randomly shuffling the data set (default value 0 means that no see is used)
 */
void ANN::learn2(unsigned maxIterations, double minError, unsigned seed)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Starting learning procedure (maxIt = " << maxIterations << ", minError = " << minError << ")" << std::endl;
  this->debugFile << this->getStats() << std::endl;
#endif
  this->maxIt = maxIterations;
  this->minError = minError;

  // initialize arrays for the deltas necessary for the backpropagation algorithm
  double* hiddenDeltas = new double[this->hiddenNodes];
  double* outputDeltas = new double[this->outputNodes];

  // initialize structure to store the previous parameter increments (necessary for the momentum term)
  double** previousHiddenParameters = NULL;
  double** previousOutputParameters = NULL;
  if (this->momentumConstant != 0)
  {
    previousHiddenParameters = new double*[this->hiddenNodes];
    for (unsigned i = 0; i < this->hiddenNodes; i++)
    {
      previousHiddenParameters[i] = new double[this->inputNodes + 1];
      for (unsigned j = 0; j <= this->inputNodes; j++)
        previousHiddenParameters[i][j] = 0;
    }

    previousOutputParameters = new double*[this->outputNodes];
    for (unsigned i = 0; i < this->outputNodes; i++)
    {
      previousOutputParameters[i] = new double[this->hiddenNodes + 1];
      for (unsigned j = 0; j <= this->hiddenNodes; j++)
        previousOutputParameters[i][j] = 0;
    }
  }

  // initialize data necessary to adjust learning constants for each parameter independently
  double** hiddenLearnConstants = new double*[this->hiddenNodes];
  double** previousHiddenSign = new double*[this->hiddenNodes];
  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    hiddenLearnConstants[i] = new double[this->inputNodes + 1];
    previousHiddenSign[i] = new double[this->inputNodes + 1];
    for (unsigned j = 0; j <= this->inputNodes; j++)
    {
      hiddenLearnConstants[i][j] = this->learningConstant;
      previousHiddenSign[i][j] = 0;
    }
  }

  double** outputLearnConstants = new double*[this->outputNodes];
  double** previousOutputSign = new double*[this->outputNodes];
  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    outputLearnConstants[i] = new double[this->hiddenNodes + 1];
    previousOutputSign[i] = new double[this->hiddenNodes + 1];
    for (unsigned j = 0; j <= this->hiddenNodes; j++)
    {
      outputLearnConstants[i][j] = this->learningConstant;
      previousOutputSign[i][j] = 0;
    }

  }

  bool done = false;
  double lastError = 0;
  unsigned iterations = 0;
  unsigned* setIndices = NULL;
  unsigned* list1 = NULL;
  unsigned* list2 = NULL;
  unsigned list1Index = 0;
  unsigned list2Index = 0;
  if (this->shuffleSet)
  {
    setIndices = new unsigned[this->samples];
    list1 = new unsigned[this->samples];
    list2 = new unsigned[this->samples];
  }
  for (unsigned i = 0; this->shuffleSet && i < this->samples; i++) setIndices[i] = i;
  if (this->shuffleSet) 
  {
    if (seed == 0)
      srand48(time(NULL));
    else
      srand48(seed);
  }
  while (!done)
  {
#ifdef ANN_DEBUG_MODE
    this->debugFile << "ITERATION " << iterations << std::endl;
#endif
    if (this->shuffleSet)
    {
      list1Index = 0;
      list2Index = 0;
    }
    //begin epoch
    for (unsigned sCount = 0; sCount < this->samples; sCount++)
    {
      unsigned s = sCount;
      if (this->shuffleSet)
        s = setIndices[sCount];
      // forward propagation
      this->updateActivation(this->learnData[s]);

      // backward propagation
      // first compute output deltas
      for (unsigned i = 0; i < this->outputNodes; i++)
      {
        double d = this->expectedOutput[s][i]; // The desired output of the i-th output node
        double o = this->outputActivation[i]; // The output of the i-th output node
        // calculation follows the formula delta[i] = error[i] * activationFunction'(sum[i]), where sum[i] = sum(w[i][j] * x[j]) (x is output of hidden layer, error[i] = d -o)
        switch (this->funcType)
        {
          case ANN_LOGISTIC:
            //activationFunction'(sum[i]) = a * y[i] * (1 - y[i]) where y[i] = output of i-th output neuron
            outputDeltas[i] = this->a * (d - o) * o * (1.0 - o);
            break;
          case ANN_TANH:
            //activationFunction'(sum[i]) = b/a * (a - y[i]) * (a + y[i]) where y[i] = output of i-th output neuron
            outputDeltas[i] = this->b/this->a * (d - o) * (this->a - o) * (this->a + o);
            break;
          default:
            assert(false);
        }
      }

      // then compute hidden deltas
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        // compute the backpropagated local gradient g[i] = sum{k}(outputDeltas[k] * outputParameters[k][i])
        double g = 0;
        for (unsigned k = 0; k < this->outputNodes; k++)
          g += outputDeltas[k] * this->outputParameters[k][i];

        double y = this->hiddenActivation[i]; // output of the i-th hidden neuron

        // calculation follow the formula delta[i] = g * activationFunction'(sum[i]), where sum[i] = sum(w[i][j] * x[j]) (x is output of input layer)
        switch (this->funcType)
        {
          case ANN_LOGISTIC:
            //activationFunction'(sum[i]) = a * y[i] * (1 - y[i]) where y[i] = output of i-th output neuron
            hiddenDeltas[i] = this->a * y * (1 - y) * g;
            break;
          case ANN_TANH:
            //activationFunction'(sum[i]) = b/a * (a - y[i]) * (a + y[i]) where y[i] = output of i-th output neuron
            hiddenDeltas[i] = this->b/this->a * (a - y) * (a + y) * g;
            break;
          default:
            assert(false);
        }
      }

      // adjust synaptic weights (parameters) for the output layer
      for (unsigned i = 0; i < this->outputNodes; i++)
      {
        for (unsigned j = 0; j <= this->hiddenNodes; j++)
        {
          double deltaW = outputDeltas[i] * this->hiddenActivation[j];
          double sign = 1.0;
          //if (deltaW < 0) sign = -1.0;
          if (outputDeltas[i] < 0) sign = -1.0;
          double product = sign * previousOutputSign[i][j];
          if (product > 0){
            outputLearnConstants[i][j] = std::min(outputLearnConstants[i][j] * this->rhoPlus, this->adjustMax);
            previousOutputSign[i][j] = sign;
          }
          else if (product < 0){
            outputLearnConstants[i][j] = std::max(outputLearnConstants[i][j] * this->rhoMinus, this->adjustMin);
            previousOutputSign[i][j] = sign;
          }
          else{
            previousOutputSign[i][j] = sign;
          }
          deltaW = outputLearnConstants[i][j] * deltaW;
          if (this->momentumConstant != 0)
          {
            deltaW += this->momentumConstant * previousOutputParameters[i][j];
            previousOutputParameters[i][j] = deltaW;
          }
          this->outputParameters[i][j] += deltaW;
        }
      }

      // adjust synaptic weights (parameters) for the hidden layer
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        for (unsigned j = 0; j <= this->inputNodes; j++)
        {
          double deltaW = hiddenDeltas[i] * this->inputActivation[j];
          double sign = 1.0;
          //if (deltaW < 0) sign = -1.0;
          if (hiddenDeltas[i] < 0) sign = -1.0;
          double product = sign * previousHiddenSign[i][j];
          if (product > 0){
            hiddenLearnConstants[i][j] = std::min(hiddenLearnConstants[i][j] * this->rhoPlus, this->adjustMax);
            previousHiddenSign[i][j] = sign;
          }
          else if (product < 0){
            hiddenLearnConstants[i][j] = std::max(hiddenLearnConstants[i][j] * this->rhoMinus, this->adjustMin);
            previousHiddenSign[i][j] = sign;
          }
          else{
            previousHiddenSign[i][j] = sign;
          }
          deltaW = hiddenLearnConstants[i][j] * deltaW;
          if (this->momentumConstant != 0)
          {
            deltaW += this->momentumConstant * previousHiddenParameters[i][j];
            previousHiddenParameters[i][j] = deltaW;
          }
          this->hiddenParameters[i][j] += deltaW;
        }
      }

      if (this->shuffleSet)
      {
        double rnd = drand48();
        if (rnd > 0.5)
          list1[list1Index++] = s;
        else
          list2[list2Index++] = s;
      }
    }

    // check if already done
#ifdef ANN_DEBUG_MODE
    double error = this->getError();
    double errorChange = lastError - error;
    if (errorChange < 0)
      errorChange = -errorChange;
#ifdef ANN_PROCENTUAL_ERROR
    if (lastError != 0) errorChange /= lastError;
#endif
#endif
    if (!this->validating)
    {
#ifndef ANN_DEBUG_MODE
      double error = 0;
#endif
      if (this->minError != 0.0)
      {
#ifndef ANN_DEBUG_MODE
        error = this->getError();
        double errorChange = lastError - error;
        if (errorChange < 0)
          errorChange = -errorChange;
#ifdef ANN_PROCENTUAL_ERROR
        if (lastError != 0) errorChange /= lastError;
#endif
        lastError = error;
#endif
        done = (errorChange <= this->minError);
#ifdef ANN_DEBUG_MODE
        this->debugFile << "done = (errorChange <= this->minError) = " << done << std::endl;
#endif
      }
      if (this->maxIt != 0)
      {
        if (this->minError != 0.0)
        {
          done = (done || (iterations >= this->maxIt - 1));
#ifdef ANN_DEBUG_MODE
          this->debugFile << "done = (done || (iterations >= this->maxIt - 1)) = " << done << std::endl;
#endif
        }
        else
        {
          done = (iterations >= this->maxIt - 1);
#ifdef ANN_DEBUG_MODE
          this->debugFile << "done = (iterations >= this->maxIt - 1) = " << done << std::endl;
#endif
        }
      }
    }
    else
    {
      this->updateValidation();
      if (this->maxIt != 0)
        done = (iterations >= this->maxIt);
      done = done || (this->validationIt >= this->validationMaxIt);
    }

#ifdef ANN_DEBUG_MODE
    this->debugFile << iterations << ": error = " << error << "; lastError = " << lastError << "; change = " << errorChange << "; minError = " << this->minError << std::endl;
    this->dataFile << iterations << "\t" << error << "\t" << lastError << "\t" << errorChange << std::endl;
    if (this->period != 0 && iterations % this->period == 0)
    {
      std::ostringstream filename;
      filename << WRITE_PREFIX_STRING << iterations << WRITE_SUFFIX_STRING;
      std::string filenameStr = filename.str();
      this->writeParameters(filenameStr);
    }
    lastError = error;
#endif

    // randomize order of training examples by merging list1 and list2
    if (this->shuffleSet)
    {
      assert(list1Index + list2Index == this->samples);
      unsigned counter = 0;
      unsigned i1 = 0;
      unsigned i2 = 0;
      while (i1 < list1Index || i2 < list2Index)
      {
        if (i1 < list1Index)
          setIndices[counter++] = list1[i1++];
        if (i2 < list2Index)
          setIndices[counter++] = list2[i2++];
      }
    }

    if (this->adjustConstants) this->adjustLearnConstant(iterations);
    iterations++;
    // TODO optimal brain surgeon?
  }

#ifdef ANN_DEBUG_MODE
  if (this->period != 0)
  {
    std::ostringstream filename;
    filename << WRITE_PREFIX_STRING << iterations << WRITE_SUFFIX_STRING;
    std::string filenameStr = filename.str();
    this->writeParameters(filenameStr);
  }
#endif

  if (this->validating)
  {
    for (unsigned i = 0; i < this->hiddenNodes; i++)
    {
      for (unsigned j = 0; j <= this->inputNodes; j++)
      {
        this->hiddenParameters[i][j] = this->bestHiddenParameters[i][j];
      }
    }
    for (unsigned i = 0; i < this->outputNodes; i++)
    {
      for (unsigned j = 0; j <= this->hiddenNodes; j++)
      {
        this->outputParameters[i][j] = this->bestOutputParameters[i][j];
      }
    }
  }

  delete [] hiddenDeltas;
  delete [] outputDeltas;
  if (this->shuffleSet)
  {
    delete [] setIndices;
    delete [] list1;
    delete [] list2;
  }

  if (this->momentumConstant != 0)
  {
    for (unsigned i = 0; i < this->hiddenNodes; i++)
      delete [] previousHiddenParameters[i];
    delete [] previousHiddenParameters;

    for (unsigned i = 0; i < this->outputNodes; i++)
      delete [] previousOutputParameters[i];
    delete [] previousOutputParameters;
  }

  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    delete [] previousHiddenSign[i];
    delete [] hiddenLearnConstants[i];
  }
  delete [] previousHiddenSign;
  delete [] hiddenLearnConstants;

  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    delete [] previousOutputSign[i];
    delete [] outputLearnConstants[i];
  }
  delete [] previousOutputSign;
  delete [] outputLearnConstants;
}

/**
 * @brief Writes all necessary parameters of the current ANN into a file.
 *
 * @param filename The name of the file where to store the parameters
 */
void ANN::writeParameters(std::string& filename)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Writing parameters to file \"" << filename << "\"..." << std::endl;
#endif
  std::ofstream outfile(filename.c_str(), std::fstream::out | std::fstream::trunc | std::fstream::binary);

  if(!outfile.is_open())
  {
    std::ostringstream errMsg;
    errMsg << "Could not open file \"" << filename << "\" for writing!" ;
    throw ANNException(errMsg.str().c_str());
  }

  unsigned uiValue = this->inputNodes;
  outfile.write((char*)(&uiValue), sizeof(unsigned));
  uiValue = this->hiddenNodes;
  outfile.write((char*)(&uiValue), sizeof(unsigned));
  uiValue = this->outputNodes;
  outfile.write((char*)(&uiValue), sizeof(unsigned));
  ActivationType fType = this->funcType;
  outfile.write((char*)(&fType), sizeof(ActivationType));

  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    for (unsigned j = 0; j <= this->inputNodes; j++)
    {
      double value = this->hiddenParameters[i][j];
      outfile.write((char*)(&value), sizeof(double));
    }
  }

  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    for (unsigned j = 0; j <= this->hiddenNodes; j++)
    {
      double value = this->outputParameters[i][j];
      outfile.write((char*)(&value), sizeof(double));
    }
  }

  outfile.close();
}

/**
 * @brief Reads all necessary parameters of the ANN from a file.
 *
 * @param filename The name of the file where the parameters are stored
 */
void ANN::readParameters(std::string& filename)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Reading parameters from file \"" << filename << "\"..." << std::endl;
#endif
  std::ifstream infile(filename.c_str(), std::fstream::in | std::fstream::binary);

  if(!infile.is_open())
  {
    std::ostringstream errMsg;
    errMsg << "Could not open file \"" << filename << "\" for reading!" ;
    throw ANNException(errMsg.str().c_str());
  }

  unsigned inRead;
  infile.read((char*)(&inRead), sizeof(unsigned));

  unsigned hnRead;
  infile.read((char*)(&hnRead), sizeof(unsigned));

  unsigned onRead;
  infile.read((char*)(&onRead), sizeof(unsigned));

  ActivationType fType;
  infile.read((char*)(&fType), sizeof(ActivationType));
  this->funcType = fType;

  if (inRead != this->inputNodes || hnRead != this->hiddenNodes || onRead != this->outputNodes)
  {
    delete [] this->inputActivation;
    delete [] this->hiddenActivation;
    delete [] this->outputActivation;
  
    for (unsigned i = 0; i < this->hiddenNodes; i++)
      delete [] this->hiddenParameters[i];
    delete [] this->hiddenParameters;

    for (unsigned i = 0; i < this->outputNodes; i++)
      delete [] this->outputParameters[i];
    delete [] this->outputParameters;

    this->inputNodes = inRead;
    this->hiddenNodes = hnRead;
    this->outputNodes = onRead;
    
    this->inputActivation = new double[this->inputNodes + 1];
    this->hiddenActivation = new double[this->hiddenNodes + 1];
    this->outputActivation = new double[this->outputNodes];
    
    for (unsigned i = 0; i < this->inputNodes; i++)
      this->inputActivation[i] = 0;
    this->inputActivation[this->inputNodes] = 1.0;
    this->hiddenActivation[this->hiddenNodes] = 1.0;

    this->hiddenParameters = new double*[this->hiddenNodes];
    for (unsigned i = 0; i < this->hiddenNodes; i++)
      this->hiddenParameters[i] = new double[this->inputNodes + 1];

    this->outputParameters = new double*[this->outputNodes];
    for (unsigned i = 0; i < this->outputNodes; i++)
      this->outputParameters[i] = new double[this->hiddenNodes + 1];
  }

  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    for (unsigned j = 0; j <= this->inputNodes; j++)
    {
      double value;
      infile.read((char*)(&value), sizeof(double));
      this->hiddenParameters[i][j] = value;
    }
  }

  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    for (unsigned j = 0; j <= this->hiddenNodes; j++)
    {
      double value;
      infile.read((char*)(&value), sizeof(double));
      this->outputParameters[i][j] = value;
    }
  }

  infile.close();
}

/**
 * @brief Calculates and returns the output activation of the current ANN.
 *
 * @param input The input vector for which the activation should be computed
 * @return A pointer to the <code>outputActivation</code> member variable
 */
const double* ANN::getActivation(double* input)
{
// #ifdef ANN_DEBUG_MODE
//   this->debugFile << "Getting activation..." << std::endl;
// #endif
  this->updateActivation(input);

  return this->outputActivation;
}

#ifndef ANN_INLINES
/**
 * @brief Sets the type of the neurons' activation functions.
 *
 * @param type The new value of the <code>funcType</code> member variable
 */
void ANN::setFuncType(ActivationType type)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting funcType to \"";
  switch (type)
  {
    case ANN_LOGISTIC:
      this->debugFile << "LOGISTIC\"" << std::endl;
      break;
    case ANN_TANH:
      this->debugFile << "TANH\"" << std::endl;
      break;
  }
#endif
  this->funcType = type;
}

/**
 * Recommended values of the learning constant are positive values between 0.01 and 0.9.
 * @brief Sets a new value for the learning constant.
 *
 * @param value The new value for the <code>learningConstant</code> member variable
 */
void ANN::setLearningConstant(double value)
{
  assert(value >= 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting learningConstant to " << value << std::endl;
#endif
  this->learningConstant = value;
}

/**
 * Recommended values of the momentum constant are positive values between 0.6 and 0.9.
 * @brief Sets a new value for the momentum constant.
 *
 * @param value The new value for the <code>momentumConstant</code> member variable
 */
void ANN::setMomentumConstant(double value)
{
  assert(value >= 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting momentumConstant to " << value << std::endl;
#endif
  this->momentumConstant = value;
}

/**
 * @brief Sets the maximum number of iteration for the validation mechanism.
 *
 * @param validationMaxIt The new value for the <code>validationMaxIt</code> member variable
 */
void ANN::setValidationMaxIt(unsigned validationMaxIt)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting validationMaxIt to " << validationMaxIt << std::endl;
#endif
  this->validationMaxIt = validationMaxIt;
}

/**
 * @brief Sets the maximum number of iterations after which the learning algorithm should stop.
 *
 * @param maxIt The new value for the <code>maxIt</code> member variable
 */
void ANN::setMaxIt(unsigned maxIt)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting maxIt to " << maxIt << std::endl;
#endif
  this->maxIt = maxIt;
}

/**
 * @brief Sets the <code>a</code> parameter of the neurons' activation function.
 *
 * @param value The new value for the <code>a</code> member variable
 */
void ANN::setA(double value)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting a to " << value << std::endl;
#endif
  this->a = value;
}

/**
 * @brief Sets the <code>b</code> parameter of the neurons' activation function.
 *
 * @param value The new value for the <code>b</code> member variable
 */
void ANN::setB(double value)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting b to " << value << std::endl;
#endif
  this->b = value;
}

/**
 * @brief Activates or deactivates the shuffling of the learning dataset during the run of the learning algorithm.
 *
 * @param shuffleSet The new value of the <code>shuffleSet</code> member variable
 */
void ANN::setShuffleSet(bool shuffleSet)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting shuffleSet to " << shuffleSet << std::endl;
#endif
  this->shuffleSet = shuffleSet;
}

/**
 * @brief Sets the lower threshold for the error used to stop the learning algorithm.
 *
 * @param minError The new value for the <code>minError</code> member variable
 */
void ANN::setMinError(double minError)
{
  assert(minError >= 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting minError to " << minError << std::endl;
#endif
  this->minError = minError;
}

/**
 * @brief Sets the frequency according to which the learning constant should be adjusted.
 *
 * @param adjustFrequency The new value for the <code>adjustFrequency</code> member variable
 */
void ANN::setAdjustFrequency(unsigned adjustFrequency)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting adjustFrequency to " << adjustFrequency << std::endl;
#endif
  this->adjustFrequency = adjustFrequency;
}

/**
 * @brief Activates or deactivates the readjustment of the learning constant.
 *
 * @param adjustConstants The new value for the <code>adjustConstants</code> member variable
 */
void ANN::setAdjustConstants(bool adjustConstants)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting adjustConstants to " << adjustConstants << std::endl;
#endif
  this->adjustConstants = adjustConstants;
}

/**
 * @brief Sets the factor to use for readjusting the learning constant.
 *
 * @param adjustFactor The new value of the <code>adjustFactor</code> member variable
 */
void ANN::setAdjustFactor(double adjustFactor)
{
  assert(adjustFactor >= 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting adjustFactor to " << adjustFactor << std::endl;
#endif
  this->adjustFactor = adjustFactor;
}

/**
 * @brief Sets the minimum value the learning constant can not fall below when readjusting it.
 *
 * @param adjustMin The new value for the <code>adjustMin</code> member variable
 */
void ANN::setAdjustMin(double adjustMin)
{
  assert(adjustMin >= 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting adjustMin to " << adjustMin << std::endl;
#endif
  this->adjustMin = adjustMin;
}

/**
 * @brief Sets the maximum value the learning constant can not rise above when readjusting it.
 *
 * @param adjustMax The new value for the <code>adjustMax</code> member variable
 */
void ANN::setAdjustMax(double adjustMax)
{
  assert(adjustMax >= 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting adjustMax to " << adjustMax << std::endl;
#endif
  this->adjustMax = adjustMax;
}

/**
 * @brief Sets the multiplicative factor used to increase the learning constant when learning using <code>ANN::learn2()</code>.
 *
 * @param rhoPlus The new value for the <code>rhoPlus</code> member variable
 */
void ANN::setRhoPlus(double rhoPlus)
{
  assert(rhoPlus > 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting rhoPlus to " << rhoPlus << std::endl;
#endif
  this->rhoPlus = rhoPlus;
}

/**
 * @brief Sets the multiplicative factor used to decrease the learning constant when learning using <code>ANN::learn2()</code>.
 *
 * @param rhoMinus The new value for the <code>rhoMinus</code> member variable
 */
void ANN::setRhoMinus(double rhoMinus)
{
  assert(rhoMinus > 0);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting rhoMinus to " << rhoMinus << std::endl;
#endif
  this->rhoMinus = rhoMinus;
}

/**
 * @brief Sets the factor representing the importance of the regularization term.
 *
 * @param lambda The new value for the <code>lambda</code> member variable
 */
void ANN::setLambda(double lambda)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting lambda to " << lambda << std::endl;
#endif
  this->lambda = lambda;
}

/**
 * @brief Sets the value of the w0 parameter used for the weight elimination term.
 *
 * @param w_zero The new value for the <code>w_zero</code> member variable
 */
void ANN::setWZero(double w_zero)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting w_zero to " << w_zero << std::endl;
#endif
  this->w_zero = w_zero;
}

/**
 * @brief Sets the value of the type of the regularization type used.
 *
 * @param regularization The new value for the <code>regularization</code> member variable
 */
void ANN::setRegularization(RegularizationType regularization)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting regularization to " << regularization << std::endl;
#endif
  this->regularization = regularization;
}

/**
 * @brief Sets the power parameter used for the approximate smoother regularization term.
 *
 * @param power The new value for the <code>power</code> member variable
 */
void ANN::setPower(unsigned power)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Setting power to " << power << std::endl;
#endif
  this->power = power;
}

/**
 * If not debugging this function does nothing.
 * @brief Sets the number of iterations after which an intermediate parameter set should be written out.
 *
 * @param period The new value for the <code>period</code> member variable
 */
void ANN::setPeriod(unsigned period)
{
#ifdef ANN_DEBUG_MODE
  this->period = period;
#endif
}

/**
 * @brief Returns the number of input nodes.
 *
 * @return The value of the <code>inputNodes</code> member variable
 */
unsigned ANN::getInputNodes() const
{
  return this->inputNodes;
}

/**
 * @brief Returns the number of hidden nodes.
 *
 * @return The value of the <code>hiddenNodes</code> member variable
 */
unsigned ANN::getHiddenNodes() const
{
  return this->hiddenNodes;
}

/**
 * @brief Returns the number output nodes.
 *
 * @return The value of the <code>outputNodes</code> member variable
 */
unsigned ANN::getOutputNodes() const
{
  return this->outputNodes;
}

/**
 * @brief Returns the type of the neurons' activation function.
 *
 * @return The value of the <code>funcType</code> member variable
 */
ActivationType ANN::getFuncType() const
{
  return this->funcType;
}

/**
 * @brief Returns the value of the learning constant.
 *
 * @return The value of the <code>learningConstant</code> member variable
 */
double ANN::getLearningConstant() const
{
  return this->learningConstant;
}

/**
 * @brief Returns the value of the momentum constant.
 *
 * @return The value of the <code>momentumConstant</code> member variable
 */
double ANN::getMomentumConstant() const
{
  return this->momentumConstant;
}

/**
 * @brief Returns the maximum number of successing iterations for which the error on the validation set may not improve before the learning algorithms stops.
 * @return The value of the <code>validationMaxIt</code> member variable
 */
unsigned ANN::getValidationMaxIt() const
{
  return this->validationMaxIt;
}

/**
 * @brief Returns the maximum number of iterations after which the learning algorithm should stop.
 *
 * @return The value of the <code>maxIt</code> member variable
 */
unsigned ANN::getMaxIt() const
{
  return this->maxIt;
}

/**
 * @brief Returns the lower threshold for the error used to stop the learning algorithm.
 *
 * @return The value of the <code>minError</code> member variable
 */
double ANN::getMinError() const
{
  return this->minError;
}

/**
 * @brief Returns the value of the <code>a</code> parameter of the neurons' activation function.
 *
 * @return The value of the <code>a</code> member variable
 */
double ANN::getA() const
{
  return this->a;
}

/**
 * @brief Returns the value of the <code>b</code> paramter of the neurons' activation function.
 *
 * @return The value of the <code>b</code> member variable
 */
double ANN::getB() const
{
  return this->b;
}

/**
 * @brief Returns a value telling whether the learning set is being shuffled after an iteration of the learning algortihm.
 *
 * @return The value of the <code>shuffleSet</code> member variable
 */
bool ANN::getShuffleSet() const
{
  return this->shuffleSet;
}

/**
 * @brief Returns the frequency according to which the learining constant is readjusted.
 *
 * @return The value of the <code>adjustFrequency</code> member variable
 */
unsigned ANN::getAdjustFrequency() const
{
  return this->adjustFrequency;
}

/**
 * @brief Returns a value telling whether the learning constant is being readjusted durin the learning procedure.
 *
 * @return The value if the <code>adjustConstants</code> member variable
 */
bool ANN::getAdjustConstants() const
{
  return this->adjustConstants;
}

/**
 * @brief Returns the multiplicative factor used to readjust the learning constant.
 *
 * @return The value of the <code>adjustFactor</code> member variable
 */
double ANN::getAdjustFactor() const
{
  return this->adjustFactor;
}

/**
 * @brief Returns the minimum possible value for the learning constant when adjusting it.
 *
 * @return The value of the <code>adjustMin</code> member variable
 */
double ANN::getAdjustMin() const
{
  return this->adjustMin;
}

/**
 * @brief Returns the maximum possible value for the learning constant when adjusting it.
 *
 * @return The value of the <code>adjustMax</code> member variable
 */
double ANN::getAdjustMax() const
{
  return this->adjustMax;
}

/**
 * @brief Returns the multiplicative factor used to increase the learning constant when learning using <code>ANN::learn2()</code>.
 *
 * @return The value of the <code>rhoPlus</code> member variable
 */
double ANN::getRhoPlus() const
{
  return this->rhoPlus;
}

/**
 * @brief Returns the multiplicative factor used to decrease the learning constant when learning using <code>ANN::learn2()</code>.
 *
 * @return The value of the <code>rhoMinus</code> member variable
 */
double ANN::getRhoMinus() const
{
  return this->rhoMinus;
}

/**
 * @brief Returns the factor representing the importance of the regularization term.
 *
 * @return The value of the <code>lambda</code> member variable
 */
double ANN::getLambda() const
{
  return this->lambda;
}

/**
 * @brief Returns the value of the w0 parameter used for the weight elimination term.
 *
 * @return The value of the <code>w_zero</code> member variable
 */
double ANN::getWZero() const
{
  return this->w_zero;
}

/**
 * @brief Returns the value of the type of the regularization type used.
 *
 * @return The value of the <code>regularization</code> member variable
 */
RegularizationType ANN::getRegularization() const
{
  return this->regularization;
}

/**
 * @brief Returns the power parameter used for the approximate smoother regularization term.
 *
 * @return The value of the <code>power</code> member variable
 */
unsigned ANN::getPower() const
{
  return this->power;
}

/**
 * @brief Returns the number of iterations after which an intermediate parameter set should be written out.
 *
 * @return The value of the <code>period</code> member variable
 */
unsigned ANN::getPeriod() const
{
  return this->period;
}

#endif

/**
 * @brief The activation function for the ANN's neurons.
 *
 * @param x The input value for the activation function
 * @return The value of the activation function given <code>x</code>
 */
double ANN::activationFunction(double x)
{
  switch (this->funcType)
  {
    case ANN_TANH: return (this->a * tanh(this->b * x));
    case ANN_LOGISTIC: return 1.0/(1.0+exp(-this->a * x));
    default: 
      assert(false);
      return x;
  }
}

/**
 * @brief Returns the classification error for one of the embedded datasets.
 *
 * @param validation If <code>true</code> the error is computed for the validation dataset otherwise for the learning dataset
 * @returns The mean squared error for the given dataset
 */
double ANN::getError(bool validation)
{
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Calculating error (validation = " << validation << ")" << std::endl;
#endif
  if (!validation)
  {
    assert (this->learnData != NULL);
    assert (this->expectedOutput != NULL);
    assert (this->samples != 0);
  }
  else
  {
    assert (this->validationData != NULL);
    assert (this->validationExpected != NULL);
    assert (this->validationSamples != 0);
  }

  double error = 0;
  unsigned N = this->samples;
  double** expected = this->expectedOutput;
  if (validation)
  {
    N = this->validationSamples;
    expected = this->validationExpected;
  }
  for (unsigned i = 0; i < N; i++)
  {
    if (validation)
      this->updateActivation(this->validationData[i]);
    else
      this->updateActivation(this->learnData[i]);

    for (unsigned j = 0; j < this->outputNodes; j++)
    {
      double e = expected[i][j] - this->outputActivation[j];
      error += (e * e);
    }
  }

  if (!validation && this->lambda != 0)
    error += this->getPenality();

  return error/(2.0*(double)N);
}

/**
 * @brief Updates the activation for all neurons in the net with respect to a given input vector.
 *
 * @param iVec The input vector to use
 */
void ANN::updateActivation(double* iVec)
{
  // TODO different activation functions for different layers (either totally differnt functions or just differently parametrized functions)
  for (unsigned i = 0; i < this->inputNodes; i++)
    this->inputActivation[i] = iVec[i];

  for (unsigned i = 0; i < this->hiddenNodes; i++)
  {
    double sum = 0;
    for (unsigned j = 0; j <= this->inputNodes; j++)
      sum += this->hiddenParameters[i][j] * this->inputActivation[j];
    this->hiddenActivation[i] = this->activationFunction(sum);
  }

  for (unsigned i = 0; i < this->outputNodes; i++)
  {
    double sum = 0;
    for (unsigned j = 0; j <= this->hiddenNodes; j++)
      sum += this->outputParameters[i][j] * this->hiddenActivation[j];
    this->outputActivation[i] = this->activationFunction(sum);
  }
}

/**
 * @brief Frees the memory allocated for the given dataset.
 *
 * @param inputData The part of the dataset containing the input data for the set
 * @param outputData The part of the dataset containing the output data for the set
 * @param setSize The size of the dataset
 */
void ANN::freeDataset(double** inputData, double** outputData, unsigned setSize)
{
  assert(inputData != NULL);
  assert(outputData != NULL);

#ifdef ANN_DEBUG_MODE
  this->debugFile << "Freeing dataset (size = " << setSize << ")" << std::endl;
#endif

  for (unsigned i = 0; i < setSize; i++)
  {
    delete [] inputData[i];
    inputData[i] = NULL;
  }
  delete [] inputData;
  inputData = NULL;

  for (unsigned i = 0; i < setSize; i++)
  {
    delete [] outputData[i];
    outputData[i] = NULL;
  }
  delete [] outputData;
  outputData = NULL;
}

/**
 * @brief Updates all internal data used for validation learning.
 */
void ANN::updateValidation()
{
  assert(this->validating);
#ifdef ANN_DEBUG_MODE
  this->debugFile << "Updating validation data..." << std::endl;
#endif
  double currentError = this->getError(true);

  if (currentError <= this->bestErrorValue)
  {
    this->bestErrorValue = currentError;
    for (unsigned i = 0; i < this->hiddenNodes; i++)
    {
      for (unsigned j = 0; j <= this->inputNodes; j++)
      {
        this->bestHiddenParameters[i][j] = this->hiddenParameters[i][j];
      }
    }
    for (unsigned i = 0; i < this->outputNodes; i++)
    {
      for (unsigned j = 0; j <= this->hiddenNodes; j++)
      {
        this->bestOutputParameters[i][j] = this->outputParameters[i][j];
      }
    }
    this->validationIt = 0;
  }
  else
    this->validationIt++;
}

/**
 * @brief Returns some statistics for the current ANN.
 *
 * @return A string containing some statistics for the current ANN like the number of parameters or the size of the used memory
 */
std::string ANN::getStats() const
{
  std::ostringstream output;
  output << "Number of input nodes: " << this->inputNodes << std::endl;
  output << "Number of hidden nodes: " << this->hiddenNodes << std::endl;
  output << "Number of output nodes: " << this->outputNodes << std::endl << std::endl;

  unsigned hiddenParamNum = (this->inputNodes + 1) * this->hiddenNodes;
  output << "Number of hidden parameters: " << hiddenParamNum << std::endl;
  unsigned outputParamNum = (this->hiddenNodes + 1) * this->outputNodes;
  output << "Number of output parameters: " << outputParamNum << std::endl;
  output << "Total: " << (hiddenParamNum + outputParamNum) << std::endl << std::endl;

  if (this->learnData != NULL)
    output << "Learn dataset size: " << this->samples << std::endl;

  if (this->validationData != NULL)
    output << "Validation dataset size: " << this->validationSamples << std::endl;

  unsigned paramMem = (hiddenParamNum + outputParamNum) * sizeof(double);
  unsigned learnDataMem = 0;
  unsigned validationDataMem = 0;
  unsigned validationParamMem = 0;
  if (this->learnData != NULL)
    learnDataMem = this->samples * (this->inputNodes + this->outputNodes) * sizeof(double);
  if (this->validationData != NULL)
  {
    validationDataMem = this->validationSamples * (this->inputNodes + this->outputNodes) * sizeof(double);
    validationParamMem = paramMem;
  }
  unsigned totalMem = paramMem + learnDataMem + validationDataMem + validationParamMem;
  output << "Memory used by parameters: " << this->memSizeStr(paramMem) << std::endl;
  if (totalMem != paramMem)
  {
    output << "Memory used by learning set: " << this->memSizeStr(learnDataMem) << std::endl;
    output << "Memory used by validation set: " << this->memSizeStr(validationDataMem) << std::endl;
    output << "Memory used by validation parameters: " << this->memSizeStr(validationParamMem) << std::endl;
    output << "Total: " << this->memSizeStr(totalMem) << std::endl;
  }
  output << std::endl;

  return output.str();
}

/**
 * @brief Returns a human readable string for the given memory size.
 *
 * @param memSize An amount of memory in bytes
 * @return A human readable string for the given memory size, using the highest unit possible
 */
std::string ANN::memSizeStr(unsigned memSize) const
{
  const std::string units [] = {"B", "kB", "MB", "GB"};
  unsigned unitNum = 4;
  double divisor = 1024;

  double size = (double) memSize;
  unsigned unit = 0;
  while (size >= divisor && unit < unitNum)
  {
    size /= divisor;
    unit++;
  }

  std::ostringstream memStr;
  memStr << size << " " << units[unit] << " [" << memSize << " " << units[0] << "]";
  return memStr.str();
}

/**
 * @brief Sets the neurons' activation function due to literature recommendations.
 */
void ANN::useLiteratureValues()
{
  this->funcType = ANN::LITERATURE_FUNC_TYPE;
  this->a = ANN::LITERATURE_A_VALUE;
  this->b = ANN::LITERATURE_B_VALUE;
}

/**
 * The adjustment ist only done if <code>iterations</code> is a multiple of the adjustment frequency and if the adjusted value is not below the minimum
 * value.
 * @brief Adjusts the learning constant by multiplying it with a constant factor.
 *
 * @param iterations The number of iterations the backpropagation algorithm has passed so far
 */
void ANN::adjustLearnConstant(unsigned iterations)
{
  if (iterations % this->adjustFrequency == 0)
  {
    double newVal = this->learningConstant * this->adjustFactor;
    if (newVal >= this->adjustMin) this->learningConstant = newVal;
  }
}

/**
 * All values in the learning/validation data is bound to the intervall [0,1] after normalization.
 * @brief Function normalizing the embedded learning and validation data.
 */
void ANN::normalizeData()
{
  for (unsigned n = 0; n < this->samples; n++)
  {
    this->normalize(this->learnData[n]);
  }

  if (this->validationData != NULL)
  {
    for (unsigned n = 0; n < this->validationSamples; n++)
    {
      this->normalize(this->validationData[n]);
    }
  }
}

/**
 * The given input vector is changed during normalization.
 * @brief Function to normalize a given input vector so all values of the input vector are bound to the intervall [0,1].
 *
 * @param sample The vector to normalize
 */
void ANN::normalize(double* sample)
{
  for (unsigned i = 0; i < this->inputNodes; i++)
    sample[i] = 1.0/(1.0+exp(-sample[i]));
    //sample[i] = tanh(this->b * sample[i]);
}

/**
 * @brief Returns the value of the regularization term given the current set of paramters.
 *
 * @return The value of the regularization term
 */
double ANN::getPenality() const
{
  double value = 0;

  switch (this->regularization)
  {
    case W_ELIMINATION:
    {
      // go through hidden nodes
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        for (unsigned j = 0; j <= this->inputNodes; j++)
        {
          double fraction = pow(this->hiddenParameters[i][j]/this->w_zero, 2);
          value += fraction/(1.0 + fraction);
        }
      }

      // got through output nodes
      for (unsigned i = 0; i < this->outputNodes; i++)
      {
        for (unsigned j = 0; j <= this->hiddenNodes; j++)
        {
          double fraction = pow(this->outputParameters[i][j]/this->w_zero, 2);
          value += fraction/(1.0 + fraction);
        }
      }
      break;
    }
    case W_DECAY:
    {
      // go through hidden nodes
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        for (unsigned j = 0; j <= this->inputNodes; j++)
        {
          value += pow(this->hiddenParameters[i][j], 2);
        }
      }

      // go through output nodes
      for (unsigned i = 0; i < this->outputNodes; i++)
      {
        for (unsigned j = 0; j <= this->hiddenNodes; j++)
        {
          value += pow(this->outputParameters[i][j], 2);
        }
      }
    }
    case W_APPROXIMATE_SMOOTHER:
    {
      // calculate norms for hidden layer parameter vectors
      double* norms = new double[this->hiddenNodes];
      double p = this->power;
      bool even = (this->power % 2 == 0);
      if (even) p = p / 2;
      for (unsigned i = 0; i < this->hiddenNodes; i++)
      {
        norms[i] = 0;
        for (unsigned j = 0; j <= this->inputNodes; j++)
        {
          norms[i] += pow(this->hiddenParameters[i][j], 2);
        }
        if (even)
        {
          norms[i] = pow(norms[i], p);
        }
        else
        {
          norms[i] = sqrt(norms[i]);
          norms[i] = pow(norms[i], p);
        }
      }

      // calcultate sums
      for (unsigned o = 0; o < this->outputNodes; o++)
      {
        for (unsigned j = 0; j < this->hiddenNodes; j++)
        {
          value += (pow(this->outputParameters[o][j], 2) * norms[j]);
        }
        value += pow(this->outputParameters[o][this->hiddenNodes], 2);
      }

      delete [] norms;
      break;
    }
  }

  return (this->lambda * value);
}

/**
 * @brief Returns the derivative of the regularization term with respect to one of the <code>ANN</code>'s parameters.
 *
 * @param output Determines whether the parameter belongs to the ouput layer
 * @param i The index of the node in the given layer
 * @param j The index of the parameter among the parameters of the given node
 * @return The value of the derivative of the regularization term given the indices of a parameter in the net
 */
double ANN::getPenalityDerivative(bool output, unsigned i, unsigned j) const
{
  double retValue = 0;
  double weight = 0;
  if (output)
    weight = this->outputParameters[i][j];
  else
    weight = this->hiddenParameters[i][j];

  switch (this->regularization)
  {
    case W_ELIMINATION:
    {
      double fraction = weight/this->w_zero;
      double denominator = 1.0 + fraction * fraction;
      denominator *= denominator;
      retValue = this->lambda * 2.0 * fraction / denominator;
      break;
    }
    case W_DECAY:
    {
      retValue = this->lambda * 2.0 * weight;
      break;
    }
    case W_APPROXIMATE_SMOOTHER:
    {
      double norm = 0;
      double p = this->power;
      if (!output) p -= 2;
      bool even = (this->power % 2 == 0);
      if (even) p = p / 2;

      for (unsigned k = 0; k <= this->inputNodes; k++)
      {
        if (output)
          norm += pow(this->hiddenParameters[j][k], 2);
        else
          norm += pow(this->hiddenParameters[i][k], 2);
      }
      if (even)
      {
        norm = pow(norm, p);
      }
      else
      {
        norm = sqrt(norm);
        norm = pow(norm, p);
      }

      if (output)
      {
        retValue = this->lambda * 2 * weight * norm;
      }
      else
      {
        for (unsigned o = 0; o < this->outputNodes; o++)
          retValue += pow(this->outputParameters[o][i], 2);
        retValue *= this->lambda * weight * this->power * norm;
      }
      break;
    }
  }

  return retValue;
}
