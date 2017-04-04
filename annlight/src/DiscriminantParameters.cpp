#include "DiscriminantParameters.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "assert.h"

const double DiscriminantParameters::RANDOM_BOUND = 0.1;

/**
 * @brief Standard constructor.
 *
 * @param vecDim The dimension of the feature vectors to use with the discriminant
 * @param random Determines whether the parameters should be initialized randomly (<code>true</code>) or set to 0 by default (<code>false</code>)
 */
DiscriminantParameters::DiscriminantParameters(unsigned vecDim, bool random)
{
  this->vecDim = vecDim;
  this->matDim = vecDim;
  this->quadParNum = (vecDim * (vecDim + 1))/2;
  this->numberOfParameters = vecDim + this->quadParNum + 1;
  this->values = new double[this->numberOfParameters];

  if (random) srand48(time(NULL));
  for (unsigned i = 0; i < this->numberOfParameters; i++)
  {
    if (random)
      this->values[i] = 2.0 * DiscriminantParameters::RANDOM_BOUND * drand48() - DiscriminantParameters::RANDOM_BOUND;
    else
      this->values[i] = 0;
  }
}

/**
 * @brief Constructor craeting a new <code>DiscriminantParameters</code> object containing random initial values.
 *
 * @param vecDim The dimension of the feature vectors to use with the discriminant
 * @param seed The seed value for the random number generator
 */
DiscriminantParameters::DiscriminantParameters(unsigned vecDim, long seed)
{
  this->vecDim = vecDim;
  this->matDim = vecDim;
  this->quadParNum = (vecDim * (vecDim + 1))/2;
  this->numberOfParameters = vecDim + this->quadParNum + 1;
  this->values = new double[this->numberOfParameters];

  srand48(seed);
  for (unsigned i = 0; i < this->numberOfParameters; i++)
    this->values[i] = 2.0 * DiscriminantParameters::RANDOM_BOUND * drand48() - DiscriminantParameters::RANDOM_BOUND;
}

/**
 * @brief Standard destructor.
 */
DiscriminantParameters::~DiscriminantParameters()
{
  delete [] this->values;
}

#ifndef DISCRIMINANT_PARAMETERS_INLINES
/**
 * @brief Sets the value of the constant term parameter.
 *
 * @param value The new value to set
 */
void DiscriminantParameters::set(double value)
{
  this->values[this->numberOfParameters - 1] = value;
}

/**
 * @brief Sets the value of an entry of the parameter vector for the linear terms given the index of this entry within this parameter vector.
 *
 * @param i The index of a row in the parameter vector
 * @param value The new value to set
 */
void DiscriminantParameters::set(unsigned i, double value)
{
  this->values[this->getIndex(i)] = value;
}

/**
 * @brief Sets the value of an entry in the parameter matrix for the quadratic terms provided the row and the column of this entry.
 *
 * @param i The index of a row in the parameter matrix
 * @param j The index of a column in the parameter matrix
 * @param value The new value to set
 */
void DiscriminantParameters::set(unsigned i, unsigned j, double value)
{
  this->values[this->getIndex(i, j)] = value;
}

/**
 * @brief Gets the value of the constant term parameter.
 *
 * @return The value of the constant term parameter
 */
double DiscriminantParameters::get() const
{
  return this->values[this->numberOfParameters - 1];
}

/**
 * @brief Gets the value of an entry of the parameter vector for the linear terms given the index of this entry within this parameter vector.
 *
 * @param i The index of a row in the parameter vector
 * @return The value of the corresponding entry
 */
double DiscriminantParameters::get(unsigned i) const
{
  return this->values[this->getIndex(i)];
}

/**
 * @brief Gets the value of an entry in the parameter matrix for the quadratic terms provided the row and the column of this entry.
 *
 * @param i The index of a row in the parameter matrix
 * @param j The index of a column in the parameter matrix
 * @return The value of the corresponding entry
 */
double DiscriminantParameters::get(unsigned i, unsigned j) const
{
  return this->values[this->getIndex(i, j)];
}

/**
 * @brief Returns the dimension of the parameter vector for the linear terms.
 *
 * @return The dimension of the parameter vector
 */
unsigned DiscriminantParameters::getVecDim() const
{
  return this->vecDim;
}

/**
 * @brief Returns the dimension of the parameter matrix for the quadratic terms.
 *
 * @return The dimension of the parameter matrix
 */
unsigned DiscriminantParameters::getMatDim() const
{
  return this->matDim;
}

/**
 * @brief Computes the index of an entry of the parameter matrix for the quadratic terms given the row and the column in the parameter matrix.
 *
 * @param i The index of a row in the parameter matrix
 * @param j The index of a column in the parameter matrix
 * @return The index of the corresponding entry in the <code>parameters</code> member variable
 */
unsigned DiscriminantParameters::getIndex(unsigned i, unsigned j) const
{
  unsigned index = 0;
  if (i >= j)
    index = this->matDim * j  - ((j + 1) * j)/2 + i;
  else
    index = this->matDim * i  - ((i + 1) * i)/2 + j;
  index += this->vecDim;
  return index;
}

/**
 * @brief Calculates the index of an entry of the parameter vector for the linear terms given the index of this entry in the parameter vector.
 *
 * @param i The index of an entry of the parameter vector
 * @return The index of the corresponding entry in the <code>parameters</code> member variable
 */
unsigned DiscriminantParameters::getIndex(unsigned i) const
{
  return i;
}

/**
 * @brief Returns the number of parameters for the given discriminant.
 * @return The the number of parameters used for the current discriminant
 *
 * @param linear Determines whether to take linear terms into account
 * @param quadratic Determines whether to take quadratic terms into account
 */
unsigned DiscriminantParameters::getNumberOfParameters(bool linear, bool quadratic) const
{
  unsigned parNum = 0;
  if (linear)
    parNum += this->vecDim;
  if (quadratic)
    parNum += this->quadParNum;
  parNum++;
  return parNum;
}

/**
 * The lefthand sied operand and the righthand side operand must both be initialized and posess the same dimensions (<code>vecDim</code>,
 * <code>matDim</code>).
 * @brief Assignment operator assigning a righthand side <code>DiscriminantParameters</code> to a lefthand side <code>DiscriminantParameters</code>.
 *
 * @param other The righthand side parameter of the operator
 * @return A reference to the modified lefthand side parameter of the operator
 */
DiscriminantParameters& DiscriminantParameters::operator=(const DiscriminantParameters& other)
{
  assert(this->numberOfParameters == other.numberOfParameters);
  assert(this->vecDim = other.vecDim);
  assert(this->matDim = other.matDim);
  for (unsigned i = 0; i < this->numberOfParameters; i++)
    this->values[i] = other.values[i];
  return *this;
}

/**
 * Both objects have to be equivalent in dimensions (<code>vecDim</code>,
 * <code>matDim</code>). None of the operands is modified when performing the operation.
 * @brief Addition operator for adding two <code>DiscriminantParameters</code> objects.
 *
 * @param other The righthand side operand
 * @return A new <code>DiscriminantParameters</code> object containing the sum of both operands
 */
DiscriminantParameters DiscriminantParameters::operator+(const DiscriminantParameters& other) const
{
  DiscriminantParameters r(this->vecDim);
  for (unsigned i = 0; i < this->numberOfParameters; i++)
    r.values[i] = this->values[i] + other.values[i];
  return r;
}

/**
 * Both <code>DiscriminantParameters</code> objects have to be equivalent in dimensions (<code>vecDim</code>, <code>matDim</code>).
 * The lefthand operand is modified during the operation.
 * @brief Adds another <code>DiscriminantParameters</code> object to the lefthand operand.
 *
 * @param other The righthand side operand
 * @return A reference to the modified lefthand side operand
 */
DiscriminantParameters& DiscriminantParameters::operator+=(const DiscriminantParameters& other)
{
  assert(this->numberOfParameters == other.numberOfParameters);
  assert(this->vecDim = other.vecDim);
  assert(this->matDim = other.matDim);
  for (unsigned i = 0; i < this->numberOfParameters; i++)
    this->values[i] += other.values[i];
  return *this;
}

/**
 * The lefthand operand is not modified during the operation.
 * @brief Multiplies a <code>DiscriminantParameters</code> object with a factor by multiplying each of the parameter values with this factor.
 *
 * @param f The factor with which to multiply
 * @return A new <code>DiscriminantParameters</code> object containing the result of the multiplication
 */
DiscriminantParameters DiscriminantParameters::operator*(double f) const
{
  DiscriminantParameters r(this->vecDim);
  for (unsigned i = 0; i < this->numberOfParameters; i++)
    r.values[i] = this->values[i] * f;
  return r;
}

/**
 * The lefthand side operand is modified during this operation.
 * @brief Multiplies a <code>DiscriminantParameters</code> object with a factor by multiplying each of the parameter values with this factor.
 *
 * @param f The factor with which to multiply
 * @return A reference to the modified lefthand operand
 */
DiscriminantParameters& DiscriminantParameters::operator*=(double f)
{
  for (unsigned i = 0; i < this->numberOfParameters; i++)
    this->values[i] *= f;
  return *this;
}
#endif

/**
 * @brief Standard output operator for the <code>DiscriminantParameters</code> class.
 *
 * @param os The output stream to which the output should be performed
 * @param param The <code>DiscriminantParameters</code> object to output
 * @return <code>os</code> after performing the output operation
 */
std::ostream& operator<<(std::ostream& os, const DiscriminantParameters& param)
{
  os << "Number of parameters: " << param.numberOfParameters << std::endl;
  os << "Matrix(" << param.matDim << "):" << std::endl;
  for (unsigned i = 0; i < param.matDim; i++)
  {
    for (unsigned j = 0; j < param.matDim; j++)
    {
      os << param.get(i,j) << " ";
    }
    os << std::endl;
  }

  os << "Vector(" << param.vecDim << "):" << std::endl;
  for (unsigned i = 0; i < param.vecDim; i++)
    os << param.get(i) << " ";
  os << std::endl;

  os << "Constant: " << param.get() << std::endl;

  return os;
}

#ifndef DISCRIMINANT_PARAMETERS_INLINES
/**
 * @brief Returns the constant term parameter.
 *
 * @return The value of the constant term parameter
 */
double DiscriminantParameters::getConstant() const
{
  return this->get();
}

/**
 * This array has <code>vecDim</code> values.
 * @brief Returns the parameter vector for the linear terms in form of an array.
 *
 * @return The parameter vector for the linear terms as array
 */
double* DiscriminantParameters::getVector() const
{
  double* vec = new double[this->vecDim];

  for (unsigned i = 0; i < this->vecDim; i++)
    vec[i] = this->get(i);

  return vec;
}

/**
 * This array has <code>matDim</code> * <code>matDim</code> entries.
 * @brief Returns the parameter matrix for the quadratic terms in form of a 2-dimensional array.
 *
 * @return The parameter matrix for the quadratic terms as 2-dimensional array
 */
double** DiscriminantParameters::getMatrix() const
{
  double** mat = new double*[this->vecDim];

  for (unsigned i = 0; i < this->matDim; i++)
  {
    mat[i] = new double[this->matDim];
    for (unsigned j = 0; j < this->matDim; j++)
      mat[i][j] = this->get(i, j);
  }

  return mat;
}
#endif

/**
 * @brief Returns the mean parameter value for this <code>DiscriminantParameters</code>.
 *
 * @param linear Determines whether the linear parameters should be taken into account
 * @param quadratic Determines whether the quadratic parameters should be taken into account
 */
double DiscriminantParameters::getMeanParameterValue(bool linear, bool quadratic) const
{
  assert (linear || quadratic);
  double sum = 0;
  double paramNum = 0;
  if (linear)
  {
    for (unsigned i = 0; i < this->vecDim; i++)
    {
      sum += (this->get(i) * this->get(i));
      paramNum += 1;
    }
  }

  if (quadratic)
  {
    for (unsigned i = 0; i < this->matDim; i++)
    {
      for (unsigned j = i; j < this->matDim; j++)
      {
        sum += (this->get(i, j) * this->get(i, j));
        paramNum += 1;
      }
    }
  }

  // The constant term
  sum += this->get() * this->get();
  paramNum += 1;

  return sqrt(sum/paramNum);
}

/**
 * @brief Computes a normalized distance between two sets of parameters.
 *
 * @param first The first set of parameters
 * @param second The second set of parameters
 * @return A normalized difference value
 */
double DiscriminantParameters::compareParameters(const DiscriminantParameters& first, const DiscriminantParameters& second)
{
  double diffSum = 0;
  double firstSum = 0;
  double secondSum = 0;

  // compute sums for linear terms
  for (unsigned i = 0; i < first.getVecDim(); i++)
  {
    double firstValue = first.get(i);
    double secondValue = second.get(i);
    double difference = firstValue - secondValue;
    firstSum += (firstValue * firstValue);
    secondSum += (secondValue * secondValue);
    diffSum += (difference * difference);
  }

  // compute sums for quadratic terms
  for (unsigned i = 0; i < first.getMatDim(); i++)
  {
    for (unsigned j = i; j < first.getMatDim(); j++)
    {
      double firstValue = first.get(i, j);
      double secondValue = second.get(i, j);
      double difference = firstValue - secondValue;
      firstSum += (firstValue * firstValue);
      secondSum += (secondValue * secondValue);
      diffSum += (difference * difference);
    }
  }

  // add constant term
  double firstValue = first.get();
  double secondValue = second.get();
  double difference = (firstValue - secondValue);
  firstSum += (firstValue * firstValue);
  secondSum += (secondValue * secondValue);
  diffSum += (difference * difference);

  return ((2 * diffSum)/(firstSum - secondSum));
}

/**
 * @brief Computes a normalized distance between two sets of parameters.
 *
 * @param first A pointer to the first set of parameters
 * @param second A pointer to the second set of parameters
 * @return A normalized difference value
 */
double DiscriminantParameters::compareParameters(const DiscriminantParameters* first, const DiscriminantParameters* second)
{
  double diffSum = 0;
  double firstSum = 0;
  double secondSum = 0;

  // compute sums for linear terms
  for (unsigned i = 0; i < first->getVecDim(); i++)
  {
    double firstValue = first->get(i);
    double secondValue = second->get(i);
    double difference = firstValue - secondValue;
    firstSum += (firstValue * firstValue);
    secondSum += (secondValue * secondValue);
    diffSum += (difference * difference);
  }

  // compute sums for quadratic terms
  for (unsigned i = 0; i < first->getMatDim(); i++)
  {
    for (unsigned j = i; j < first->getMatDim(); j++)
    {
      double firstValue = first->get(i, j);
      double secondValue = second->get(i, j);
      double difference = firstValue - secondValue;
      firstSum += (firstValue * firstValue);
      secondSum += (secondValue * secondValue);
      diffSum += (difference * difference);
    }
  }

  // add constant term
  double firstValue = first->get();
  double secondValue = second->get();
  double difference = (firstValue - secondValue);
  firstSum += (firstValue * firstValue);
  secondSum += (secondValue * secondValue);
  diffSum += (difference * difference);

  return ((2 * diffSum)/(firstSum - secondSum));
}
