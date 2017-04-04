#ifndef _SIGNOIDFUNCTION_H_
#define _SIGNOIDFUNCTION_H_

/**
 * The represented Sigmoid functions intersect the y-axis at a value of 0.5. The values of the Sigmoid functions are always between 0 and 1.
 * @brief Implementation of a set of Sigmoid functions.
 *
 */
class SigmoidFunction
{
  public:
    SigmoidFunction(double alpha, short type);
    ~SigmoidFunction();

    double* getValuesForRange(double start, double end, unsigned steps);
    double getValue(double x);

    /** @brief Sigmoid function type for actustangent. */
    static const short ARCTAN = 0;
    /** @brief Sigmoid function type for a logistic function. */
    static const short LOGISTIC = 1;
    /** @brief Sigmoid function type for a capped linear function. */
    static const short LINEAR = 2;
    /** @brief Sigmoid function type for a step function. */
    static const short STEP = 3;
  protected:
    /** @brief The free parameter of the Sigmoid function. */
    double alpha;
    /** @brief The type of the Sigmoid function. */
    short type;
};

#endif
