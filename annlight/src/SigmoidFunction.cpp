#include "SigmoidFunction.h"
#include "math.h"
#include <iostream>

/** Definition of Pi. */
const double PI = 3.14159265;

/**
 * @brief Standard constructor creating a new <code>SigmoidFunction</code>.
 *
 * @param alpha The value for the free parameter of the new <code>SigmoidFunction</code>
 * @param type The type of the <code>SigmoidFunction</code>
 */
SigmoidFunction::SigmoidFunction(double alpha, short type)
{
  this->alpha = alpha;
  this->type = type;
}

/**
 * @brief Standard destructor.
 */
SigmoidFunction::~SigmoidFunction()
{
}

/**
 * @brief Creates an array containing a number of values of the current <code>SigmoidFunction</code> for the given interval.
 *
 * @param start The lower bound of the interval
 * @param end The upper bound for the interval
 * @param steps The number of values between <code>start</code> and <code>end</code> to be computed
 * @return An array containing <code>steps</code> equidistant values between <code>start</code> and <code>end</code>
 */
double* SigmoidFunction::getValuesForRange(double start, double end, unsigned steps)
{
  double delta = (end - start)/((double)steps);
  double* values = new double[steps];
  double x = start;
  unsigned step = 0;

  while (x <= end && step < steps)
  {
    values[step] = this->getValue(x);
    x += delta;
    step++;
  }

  return values;
}

/**
 * @brief Computes a value of the current <code>SigmoidFunction</code>.
 *
 * @param x The value for which the <code>SigmoidFunction</code> should be computed
 * @return The requested value of the <code>SigmoidFunction</code>
 */
double SigmoidFunction::getValue(double x)
{
  switch (this->type)
  {
    case ARCTAN:
      return (1.0/PI * atan(this->alpha * x) + 0.5);
    case LOGISTIC:
      return (1.0/(1.0 + exp(-this->alpha * x)));
    case LINEAR:
    {
      double y = this->alpha * x + 0.5;
      if (y < 0)
        return 0;
      else if (y > 1)
        return 1;
      else
        return y;
    }
    case STEP:
    {
      if (x >= this->alpha)
        return 1;
      else
        return 0;
    }
    default:
      std::cout << "Unknown Sigmoid function type: " << this->type << std::endl;
      exit(1);
  }
}
