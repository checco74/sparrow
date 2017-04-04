#include <iostream>

#ifndef _DISCRIMINANTPARAMETERS_H_
#define _DISCRIMINANTPARAMETERS_H_

//#define DISCRIMINANT_PARAMETERS_INLINES

#ifdef DISCRIMINANT_PARAMETERS_INLINES
#include "assert.h"
#endif

/**
 * The discriminant function has a form similar to the function y=xAx+wx+b.
 * Here x is an input (feature) vector, A is a parameter matrix for the quadratic term xAx, w is a parameter vector for the linear term wx and
 * b is a scalar constant. Only the entries of A and w are saved within this class. Since A is a quadratic, symmetric matrix, only the entries above
 * the diagonal as well as the entries on the diagonal need to be saved to reconstruct A.
 *
 * @brief Class representing the parameters of a quadratic discriminant function.
 */
class DiscriminantParameters
{
  public:
    DiscriminantParameters(unsigned vecDim, bool random = false);
    DiscriminantParameters(unsigned vecDim, long seed);
    ~DiscriminantParameters();

#ifndef DISCRIMINANT_PARAMETERS_INLINES
    void set(double value);
    void set(unsigned i, double value);
    void set(unsigned i, unsigned j, double value);

    double get() const;
    double get(unsigned i) const;
    double get(unsigned i, unsigned j) const;

    unsigned getVecDim() const;
    unsigned getMatDim() const;

    double getConstant() const;
    double* getVector() const;
    double** getMatrix() const;

    unsigned getNumberOfParameters(bool linear, bool quadratic) const;

    DiscriminantParameters& operator=(const DiscriminantParameters& other);
    DiscriminantParameters operator+(const DiscriminantParameters& other) const;
    DiscriminantParameters& operator+=(const DiscriminantParameters& other);
    DiscriminantParameters operator*(double f) const;
    DiscriminantParameters& operator*=(double f);
#else
    inline void set(double value) {this->values[this->numberOfParameters - 1] = value;}
    inline void set(unsigned i, double value) {this->values[this->getIndex(i)] = value;}
    inline void set(unsigned i, unsigned j, double value) {this->values[this->getIndex(i, j)] = value;}

    inline double get() const {return this->values[this->numberOfParameters - 1];}
    inline double get(unsigned i) const {return this->values[this->getIndex(i)];}
    inline double get(unsigned i, unsigned j) const {return this->values[this->getIndex(i, j)];}

    inline unsigned getVecDim() const {return this->vecDim;}
    inline unsigned getMatDim() const {return this->matDim;}

    inline double getConstant() const {return this->get();}

    inline double* getVector() const
    {
      double* vec = new double[this->vecDim];

      for (unsigned i = 0; i < this->vecDim; i++)
        vec[i] = this->get(i);

      return vec;
    }

    inline double** getMatrix() const
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

    inline unsigned getNumberOfParameters(bool linear, bool quadratic) const
    {
      unsigned parNum = 0;
      if (linear)
        parNum += this->vecDim;
      if (quadratic)
        parNum += this->quadParNum;
      parNum++;
      return parNum;
    }

    inline DiscriminantParameters& operator=(const DiscriminantParameters& other)
    {
      assert(this->numberOfParameters == other.numberOfParameters);
      assert(this->vecDim = other.vecDim);
      assert(this->matDim = other.matDim);
      for (unsigned i = 0; i < this->numberOfParameters; i++)
        this->values[i] = other.values[i];
      return *this;
    }

    inline DiscriminantParameters operator+(const DiscriminantParameters& other) const
    {
      DiscriminantParameters r(this->vecDim);
      for (unsigned i = 0; i < this->numberOfParameters; i++)
        r.values[i] = this->values[i] + other.values[i];
      return r;
    }

    inline DiscriminantParameters& operator+=(const DiscriminantParameters& other)
    {
      assert(this->numberOfParameters == other.numberOfParameters);
      assert(this->vecDim = other.vecDim);
      assert(this->matDim = other.matDim);
      for (unsigned i = 0; i < this->numberOfParameters; i++)
        this->values[i] += other.values[i];
      return *this;
    }

    inline DiscriminantParameters operator*(double f) const
    {
      DiscriminantParameters r(this->vecDim);
      for (unsigned i = 0; i < this->numberOfParameters; i++)
        r.values[i] = this->values[i] * f;
      return r;
    }

    inline DiscriminantParameters& operator*=(double f)
    {
      for (unsigned i = 0; i < this->numberOfParameters; i++)
        this->values[i] *= f;
      return *this;
    }
#endif

    double getMeanParameterValue(bool linear, bool quadratic) const;

    static double compareParameters(const DiscriminantParameters& first, const DiscriminantParameters& second);
    static double compareParameters(const DiscriminantParameters* first, const DiscriminantParameters* second);

    friend std::ostream& operator<<(std::ostream& os, const DiscriminantParameters& param);
  protected:
#ifndef DISCRIMINANT_PARAMETERS_INLINES
    unsigned getIndex(unsigned i, unsigned j) const;
    unsigned getIndex(unsigned i) const;
#else
    inline unsigned getIndex(unsigned i, unsigned j) const 
    {
      unsigned index = 0;
      if (i >= j)
        index = this->matDim * j  - ((j + 1) * j)/2 + i;
      else
        index = this->matDim * i  - ((i + 1) * i)/2 + j;
      index += this->vecDim;
      return index;
    }
    inline unsigned getIndex(unsigned i) const {return i;}
#endif

    /** This array contains both the parameters for the linear terms (parameter vector, the first entries of the array) and the parameters for the quadratic terms (parameter matrix, only the values above the diagonal and the diagonal itself). @brief The values of the parameters of the discriminant function. */
    double* values;
    /** @brief The total number of parameters saved in <code>values</code>. */
    unsigned numberOfParameters;
    /** @brief The number of quadratic parameters saved in <code>values</code>. */
    unsigned quadParNum;
    /** @brief The dimension (number of entries) of the parameter vector for the linear terms. */
    unsigned vecDim;
    /** @brief The dimension (number of rows) of the parameter matrix for the quadratic terms. */
    unsigned matDim;
    /** If parameters are initialized randomly the parameter values lie between <code>RANDOM_BOUND</code> and -<code>RANDOM_BOUND</code>. */
    static const double RANDOM_BOUND;
};

#endif
