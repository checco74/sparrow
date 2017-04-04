#include "SVDSolver.h"
#include "assert.h"
#include <stdlib.h>
#include "math.h"
#include <sstream>
#include "SVDException.h"
#include <iostream>

#define sign(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/**
 * Creates a new <code>SVDSolver</code> object with the given linear equation system data.
 * @brief Constructor for the <code>SVDSolver</code> class.
 *
 * @param covMatrix The lefthand coefficient matrix of the linear equation system (in this case a covariance matrix) as one-dimensional array
 * @param righthand The righthand side vector of the linear equation system
 * @param dim The dimension of the covariance matrix as well as the solution vector
 * @param copy Determines whether the given arrays should be copied or used directly
 */
SVDSolver::SVDSolver(double* covMatrix, double* righthand, unsigned dim, bool copy)
{
  this->dim = dim;
  this->startIndex = dim;
  
  this->eigenvalues = NULL;
  this->eigenvectors = NULL;
  this->offDiagonalElements = NULL;
  this->solution = NULL;
  
  if (copy)
  {
#ifndef SVD_SYMMETRIC_STORAGE
    unsigned matDim = this->dim * this->dim;
#else
    unsigned matDim = ((this->dim + 1) * this->dim)/2;
#endif
    this->covMatrix = new double[matDim];
    for (unsigned i = 0; i < matDim; i++)
      this->covMatrix[i] = covMatrix[i];
      
    this->righthand = new double[this->dim];
    for (unsigned i = 0; i < this->dim; i++)
      this->righthand[i] = righthand[i];
  }
  else
  {
    this->covMatrix = covMatrix;
    this->righthand = righthand;
  }
}

/**
 * Creates a new <code>SVDSolver</code> object with the given linear equation system data.
 * @brief Constructor for the <code>SVDSolver</code> class.
 *
 * @param covMatrix The lefthand coefficient matrix of the linear equation system (in this case a covariance matrix) as two-dimensional array
 * @param righthand The righthand side vector of the linear equation system
 * @param dim The dimension of the covariance matrix as well as the solution vector
 * @param copy Determines whether <code>righthand</code> should be copied or used directly
 */
SVDSolver::SVDSolver(double** covMatrix, double* righthand, unsigned dim, bool copy)
{
  this->dim = dim;
  this->startIndex = dim;
  
  this->eigenvalues = NULL;
  this->eigenvectors = NULL;
  this->offDiagonalElements = NULL;
  this->solution = NULL;
  
#ifndef SVD_SYMMETRIC_STORAGE
  unsigned matDim = this->dim * this->dim;
#else
  unsigned matDim = ((this->dim + 1) * this->dim)/2;
#endif
  this->covMatrix = new double[matDim];
  for (unsigned i = 0; i < matDim; i++)
    for (unsigned j = 0; j < matDim; j++)
      this->setCovariance(i, j, covMatrix[i][j]);
  
  if (copy)
  {
    this->righthand = new double[this->dim + 1];
    for (unsigned i = 0; i <= this->dim; i++)
      this->righthand[i] = righthand[i];
  }
  else
    this->righthand = righthand;
}

/**
 * @brief Destructor of the <code>SVDSolver</code> class.
 */
SVDSolver::~SVDSolver()
{
  delete [] this->covMatrix;
  delete [] this->righthand;
  if (this->eigenvalues != NULL) delete [] this->eigenvalues;
  if (this->eigenvectors != NULL) delete [] this->eigenvectors;
  if (this->offDiagonalElements != NULL) delete [] this->offDiagonalElements;
  if (this->solution != NULL) delete [] this->solution;
}

/**
 * The solving method describes which eigenvalues in the diagonal coefficient matrix should be set to zero. The methods are:<br>
 * <ul>
 *   <li>TAKE_BEST: Takes only the <code>value</code> highest eigenvalues into account setting all other eigenvalues to zero</li>
 *   <li>TAKE_ABOVE_THRESHOLD: Takes only those eigenvalues into account that are larger than <code>value</code>, setting all other eigenvalues to zero</li>
 *   <li>LEAVE_OUT_WORST: Leaves out the <code>value</code> lowet eigenvalues setting them to zero</li>
 *   <li>LEAVE_OUT_BELOW_THRESHOLD: Leaves out all eigenvalues that are lower than <code>value</code> by setting them to zero</li>
 *   <li>TAKE_PERCENTAGE: Takes only the highest eigenvalues into account which summed up and devided by the sum of all eigenvalues are above <code>value</code></li>
 * </ul>
 * @brief Solves the linear equation system using the given solvin method.
 *
 * @param method The method to determine which eigenvalues should be left out of the computation
 * @param value The value used to determine which eigenvalues should be left out of the computation
 * @throw SVDException If <code>value</code> is not within the specified bounds
 */
void SVDSolver::solve(t_solving_type method, double value)
{
  if (method == TAKE_BEST && (value < 0 || (unsigned)value > this->dim))
  {
    std::ostringstream errMsg;
    errMsg << "Threshold value out of bounds: min = 0, max = " << this->dim << ", current = " << value;
    throw SVDException(errMsg.str().c_str());
  }
  if (method == LEAVE_OUT_WORST && (value < 0 || (unsigned)value > this->dim))
  {
    std::ostringstream errMsg;
    errMsg << "Threshold value out of bounds: min = 0, max = " << this->dim << ", current = " << value;
    throw SVDException(errMsg.str().c_str());
  }
  if (method == TAKE_PERCENTAGE && (value < 0 || value > 1))
  {
    std::ostringstream errMsg;
    errMsg << "Threshold value must be between 0 and 1! Current value is " << value;
    throw SVDException(errMsg.str().c_str());
  }
  
  this->calculateEigenvalues();  
  
  if ((method == TAKE_BEST && (unsigned)value == this->dim || method == LEAVE_OUT_WORST && value == 0) )
    this->startIndex = 0;
  else
  {
    unsigned start = 0;
    switch (method)
    {
      case TAKE_BEST:
        start = this->dim - (unsigned)value;
        break;
      case TAKE_ABOVE_THRESHOLD:
        assert(this->eigenvalues != NULL);
        for (unsigned i = 0; i < this->dim && this->eigenvalues[i] <= value; i++) start++;
        if (start >= this->dim) start = this->dim - 1;  // At least one eigenvalue should be taken
        break;
      case LEAVE_OUT_WORST:
        start = (unsigned)value;
        break;
      case LEAVE_OUT_BELOW_THRESHOLD:
        assert(this->eigenvalues != NULL);
        for (unsigned i = 0; i < this->dim && this->eigenvalues[i] < value; i++) start++;
        if (start >= this->dim) start = this->dim - 1;
        break;
      case TAKE_PERCENTAGE:
      {
        double sum = 0;
        for (unsigned i = 0; i < this->dim; i++)
          sum += this->eigenvalues[i];
        assert(sum != 0);
        double sum2 = 0;
        start = this->dim;
        for (unsigned i = 0; i < this->dim && sum2 < value; i++)
        {
          sum2 += this->eigenvalues[this->dim - i - 1]/sum;
          start--;
        }
      }
    }
    this->startIndex = start;
  }
  
  this->solve();
}

/**
 * This function must be called before <code>tql2()</code> can be called.
 * @brief Householder reduction of a real symmetric matrix as seen in numerical recipies.
 */
void SVDSolver::tred2()
{
  unsigned i,j,k,l,ii,jp1;
  double f,g,h,hh,scale;
  //double fabs(), sign(), sqrt();

  for (i = 1; i <= this->dim; i++)
  {
    for (j = i; j <= this->dim; j++)
      this->setEigenvectorValue(j,i, this->getCovariance(j, i));
    this->setEigenvalue(i, this->getCovariance(this->dim, i));
  }
  if (this->dim!=1) 
  {
    for (ii = 2; ii <= this->dim; ii++)
    {
      i = this->dim + 2 - ii;
      l = i - 1;
      h = 0.0;
      scale = 0.0;
      if (l >= 2) 
      {
        for (k = 1; k <= l; k++)
          scale = double(scale+fabs(this->getEigenvalue(k)));
        if (scale != 0.0) 
        {
          for (k = 1; k <= l; k++)
          {
            this->setEigenvalue(k, this->getEigenvalue(k)/scale);
            h = double(h+this->getEigenvalue(k)*this->getEigenvalue(k));
          }
          f = this->getEigenvalue(l);
          g = -sign(double(sqrt(h)),f);
          this->setOffDiagonalElement(i, scale*g);
          h = double(h-f*g);
          this->setEigenvalue(l, double(f-g));
          for (j = 1; j <= l; j++)
            this->setOffDiagonalElement(j, 0.0);
          for (j = 1; j <= l; j++)
          {
            f = this->getEigenvalue(j);
            this->setEigenvectorValue(j, i, f);
            g = double(this->getOffDiagonalElement(j)+this->getEigenvectorValue(j,j)*f);
            jp1 = j+1;
            if (l>=jp1)
              for (k = jp1; k <= l; k++)
              {
                g = double(g+this->getEigenvectorValue(k,j)*this->getEigenvalue(k));
                this->setOffDiagonalElement(k, double(this->getOffDiagonalElement(k)+this->getEigenvectorValue(k,j)*f));
              }
            this->setOffDiagonalElement(j, g);
          }
          f = 0.0;
          for (j = 1; j <= l; j++)
          {
            this->setOffDiagonalElement(j, this->getOffDiagonalElement(j)/h);
            f = double(f+this->getOffDiagonalElement(j)*this->getEigenvalue(j));
          }
          hh = f/(double(h+h));
          for (j = 1; j <= l; j++)
            this->setOffDiagonalElement(j, double(this->getOffDiagonalElement(j)-hh*this->getEigenvalue(j)));
          for (j = 1; j <= l; j++)
          {
            f = this->getEigenvalue(j);
            g = this->getOffDiagonalElement(j);
            for (k = j; k <= l; k++)
              this->setEigenvectorValue(k, j, double(this->getEigenvectorValue(k,j)-f*this->getOffDiagonalElement(k)-g*this->getEigenvalue(k)));
            this->setEigenvalue(j, this->getEigenvectorValue(l,j));
            this->setEigenvectorValue(i, j, 0.0);
          }
          goto lab10;
        }
      }
      this->setOffDiagonalElement(i, this->getEigenvalue(l)); 
      for (j = 1; j <= l; j++)
      {
        this->setEigenvalue(j, this->getEigenvectorValue(l,j));
        this->setEigenvectorValue(i, j, 0.0);
        this->setEigenvectorValue(j, i, 0.0);
      }lab10:  this->setEigenvalue(i, h);
    }
    for (i = 2; i <= this->dim; i++)
    {
      l = i - 1;
      this->setEigenvectorValue(this->dim, l, this->getEigenvectorValue(l,l));
      this->setEigenvectorValue(l, l, 1.0);
      h = this->getEigenvalue(i);
      if (h!=0.0) 
      {
        for (k = 1; k <= l; k++)
        this->setEigenvalue(k, this->getEigenvectorValue(k,i)/h);
        for (j = 1; j <= l; j++)
        {
          g = 0.0;
          for (k = 1; k <= l; k++)
            g = double(g+this->getEigenvectorValue(k,i)*this->getEigenvectorValue(k,j));
            for (k = 1; k <= l; k++)
              this->setEigenvectorValue(k, j, double(this->getEigenvectorValue(k,j)-g*this->getEigenvalue(k)));
        }
      }
      for (k = 1; k <= l; k++)
        this->setEigenvectorValue(k, i, 0.0);
    }
  }
  for (i = 1; i <= this->dim; i++)
  {
    this->setEigenvalue(i, this->getEigenvectorValue(this->dim, i));
    this->setEigenvectorValue(this->dim, i, 0.0);
  }
  this->setEigenvectorValue(this->dim, this->dim, 1.0);
  this->setOffDiagonalElement(1, 0.0);
}

/**
 * The algorithm calculates the eigenvalues and eigenvectors of the coefficient matrix assuming that the function <code>tred2()</code> has been called beforehand.
 * @brief An implementation of the QL algorithm as seen in numerical recipies.
 *
 * @param ierr Variable containing information about errors that occured while executing the function
 */
void SVDSolver::tql2(int *ierr)
{
  unsigned i,j,k,l,m,ii,l1,l2,mml;
  double c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2;
  s2 = 0.0;
  c3 = 0.0;

  *ierr = 0;
  if( this->dim != 1 ) 
  {
    for (i = 2; i <= this->dim; i++)
      this->setOffDiagonalElement(i-1, this->getOffDiagonalElement(i));
    f = 0.0;
    tst1 = 0.0;
    this->setOffDiagonalElement(this->dim, 0.0);
    for (l = 1; l <= this->dim; l++)
    {
      j = 0;
      h = double(fabs(this->getEigenvalue(l))+fabs(this->getOffDiagonalElement(l)));
      if (tst1<h)
        tst1 = h;
      for (m = l; m <= this->dim; m++)
      {
        tst2 = double(tst1+fabs(this->getOffDiagonalElement(m)));
        if (tst2==tst1)
          break;
      }
      if (m!=l)
        do 
        {
          if (j==30)
            goto lab10;
          j  = j+1;
          l1 = l+1;
          l2 = l1+1;
          g  = this->getEigenvalue(l);
          p  = double((this->getEigenvalue(l1)-g))/(2.0*this->getOffDiagonalElement(l));
          r  = double(hypot(p,1.0));
          this->setEigenvalue(l, this->getOffDiagonalElement(l)/(double(p+sign(r,p))));
          this->setEigenvalue(l1, this->getOffDiagonalElement(l)*(double(p+sign(r,p))));
          dl1 = this->getEigenvalue(l1);
          h = double(g-this->getEigenvalue(l));
          if (l2 <= this->dim)
            for (i = l2; i <= this->dim; i++)
              this->setEigenvalue(i, double(this->getEigenvalue(i)-h));
          f = double(f+h);
          p = this->getEigenvalue(m);
          c = 1.0;
          c2 = c;
          el1 = this->getOffDiagonalElement(l1);
          s = 0.0;
          mml = m-l;
          for (ii = 1; ii <= mml; ii++)
          {
            c3 = c2;
            c2 = c;
            s2 = s;
            i = m-ii;
            g = c*this->getOffDiagonalElement(i);
            h = c*p;
            r = double(hypot(p,this->getOffDiagonalElement(i)));
            this->setOffDiagonalElement(i+1, s*r);
            s = this->getOffDiagonalElement(i)/r;
            c = p/r;
            p = double(c*this->getEigenvalue(i)-s*g);
            this->setEigenvalue(i+1, double(h+s*(c*g+s*this->getEigenvalue(i))));
            for (k = 1; k <= this->dim; k++)
            {
              h = this->getEigenvectorValue(k,i+1);
              this->setEigenvectorValue(k, i+1, double(s*this->getEigenvectorValue(k,i)+c*h));
              this->setEigenvectorValue(k, i, double(c*this->getEigenvectorValue(k,i)-s*h));
            }
          }
          p = -s*s2*c3*el1*this->getOffDiagonalElement(l)/dl1;
          this->setOffDiagonalElement(l, s*p);
          this->setEigenvalue(l, c*p);
          tst2 = double(tst1+fabs(this->getOffDiagonalElement(l)));  // fabs(e(l))*THRESHOLD
        }
        while(tst2>tst1);
          this->setEigenvalue(l, double(this->getEigenvalue(l)+f));
      }
      for (ii = 2; ii <= this->dim; ii++)
      {
        i = ii-1;
        k = i;
        p = this->getEigenvalue(i);
        for (j = ii; j <= this->dim; j++)
          if (this->getEigenvalue(j)<p) 
          {
            k = j;
            p = this->getEigenvalue(j);
          }
        if (k!=i) 
        {
          this->setEigenvalue(k, this->getEigenvalue(i));
          this->setEigenvalue(i, p);
          for (j = 1; j <= this->dim; j++)
          {
            p = this->getEigenvectorValue(j, i);
            this->setEigenvectorValue(j, i, this->getEigenvectorValue(j,k));
            this->setEigenvectorValue(j, k, p);
          }
        }
      }
      return;
lab10:  *ierr = l;
  }
  return;
}

/**
 * @brief Returns an entry of the covariance matrix.
 *
 * @param i The index of the row in which the entry is located
 * @param j The index of the column in which the entry is located
 * @param start1 Determines whether the matrix indices should be treated as starting with 1 or with 0
 * @return The value of the desired matrix entry
 * @throw SVDException If indices are out of bounds
 */
double SVDSolver::getCovariance(unsigned i, unsigned j, bool start1) const
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Row index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1 && j < 1 || start1 && j > this->dim || !start1 && j >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Column index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1)
#ifndef SVD_SYMMETRIC_STORAGE
    return this->covMatrix[(i - 1) + (j - 1) * this->dim];
#else
    return this->covMatrix[this->getIndex(i - 1, j - 1)];
#endif
  else
#ifndef SVD_SYMMETRIC_STORAGE
    return this->covMatrix[i + j * this->dim];
#else
    return this->covMatrix[this->getIndex(i, j)];
#endif
}

/**
 * @brief Sets an entry of the covariance matrix.
 *
 * @param i The index of the row in which the entry is located
 * @param j The index of the column in which the entry is located
 * @param value The value to which the entry should be set
 * @param start1 Determines whether the matrix indices should be treated as starting with 1 or with 0
 * @throw SVDException If passed indices are out of bounds
 */
void SVDSolver::setCovariance(unsigned i, unsigned j, double value, bool start1)
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Row index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1 && j < 1 || start1 && j > this->dim || !start1 && j >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Column index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1)
#ifndef SVD_SYMMETRIC_STORAGE
    this->covMatrix[(i - 1) + (j - 1) * this->dim] = value;
#else
    this->covMatrix[this->getIndex(i - 1, j - 1)] = value;
#endif
  else
#ifndef SVD_SYMMETRIC_STORAGE
    this->covMatrix[i + j * this->dim] = value;
#else
    this->covMatrix[this->getIndex(i, j)] = value;
#endif
}

/**
 * @brief Returns the value of the offdiagonal elements.
 *
 * @param i The index of the offdiagonal element
 * @param start1 Determines whether the array indices should be treated as starting with 1 or with 0
 * @return The value of the desired offdiagonal element
 * @throw SVDException If <code>i</code> is out of bounds
 */
double SVDSolver::getOffDiagonalElement(unsigned i, bool start1) const
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1)
    return this->offDiagonalElements[i-1];
  else
    return this->offDiagonalElements[i];
    
}

/**
 * @brief Sets the value of the offdiagonal elements.
 *
 * @param i The index of the offdiagonal element
 * @param value The value to which the offdiagonal element should be set
 * @param start1 Determines whether the array indices should be treated as starting with 1 or with 0
 * @throw SVDException If <code>i</code> is out of bounds
 */
void SVDSolver::setOffDiagonalElement(unsigned i, double value, bool start1)
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1)
    this->offDiagonalElements[i-1] = value;
  else
    this->offDiagonalElements[i] = value;
}

/**
 * @brief Returns the value of an entry of a certain eigenvector.
 *
 * @param i The index of the desired component within the eigenvector
 * @param j The index of the eigenvector the desired component belongs to
 * @param start1 Determines whether the passed indices should be treated as starting with 1 or with 0
 * @return The value of the desired eigenvecor component
 * @throw SVDException If passed indices are out of bounds
 */
double SVDSolver::getEigenvectorValue(unsigned i, unsigned j, bool start1) const
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Component index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1 && j < 1 || start1 && j > this->dim || !start1 && j >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvector index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1)
    return this->eigenvectors[(i - 1) + (j - 1) * this->dim];
  else
    return this->eigenvectors[i + j * this->dim];
}

/**
 * @brief Returns a certain eigenvector.
 *
 * @param i The index of the eigenvector
 * @param start1 Determines whether the eigenvector indices should be treated as starting with 1 or with 0
 * @return The desired eigenvector as array with length <code>dim</code>
 * @throw SVDException If the eigenvector index is out of bounds
 */
double* SVDSolver::getEigenvector(unsigned i, bool start1) const
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvector index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }

  double* vec = new double[this->dim];
  
  unsigned start = 0;
  unsigned offset = 0;
  unsigned stop = this->dim;
  if (start1)
  {
    start++;
    offset++;
    stop++;
  }
  for (unsigned h = start; h < stop; h++)
    vec[h - offset] = this->getEigenvectorValue(h, i, start1);
  
  return vec;
}

/**
 * @brief Sets the value of an entry of a certain eigenvector.
 *
 * @param i The index of the desired component within the eigenvector
 * @param j The index of the eigenvector the desired component belongs to
 * @param value The value the desired component should be set to
 * @param start1 Determines whether the passed indices should be treated as starting with 1 or with 0
 * @throw SVDException If passed indices are out of bounds
 */
void SVDSolver::setEigenvectorValue(unsigned i, unsigned j, double value, bool start1)
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Component index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (start1 && j < 1 || start1 && j > this->dim || !start1 && j >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvector index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }

  if (start1)
    this->eigenvectors[(i - 1) + (j - 1) * this->dim] = value;
  else
    this->eigenvectors[i + j * this->dim] = value;
}

/**
 * @brief Returns the value of a desired eigenvalue
 *
 * @param i The index of the desired eigenvalue
 * @param start1 Determines whether the eigenvalue indices should be treated as starting with 1 or with 0
 * @return The value of the desired eigenvalue
 * @throw SVDException If eigenvalue index is out of bounds
 */
double SVDSolver::getEigenvalue(unsigned i, bool start1) const
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvalue index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  
  if (start1)
    return this->eigenvalues[i - 1];
  else
    return this->eigenvalues[i];
}

/**
 * @brief Sets the value of a certain eigenvalue
 *
 * @param i The index of the desired eigenvalue
 * @param value The value the desired eigenvalue should be set to
 * @param start1 Determines whether the eigenvalue indices should be treated as starting with 1 or with 0
 * @throw SVDException If eigenvalue index is out of bounds
 */
void SVDSolver::setEigenvalue(unsigned i, double value, bool start1)
{
  if (start1 && i < 1 || start1 && i > this->dim || !start1 && i >= this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvalue index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }

  if (start1)
    this->eigenvalues[i - 1] = value;
  else
    this->eigenvalues[i] = value;
}

/**
 * @brief Solves the linear equation system assuming that the eigenvalues an eigenvectors of the coefficient matrix heve bee calculated.
 */
void SVDSolver::solve()
{
  // TODO Das ganze richtig proggen
  if (this->solution == NULL) this->solution = new double[this->dim];

  double *bI	= new double [this->dim];
    
  for (unsigned i = 0; i < this->dim; i++)
  {
    bI[i]	= 0.;
    for (unsigned j = 0; j < this->dim; j++)
      bI[i] += this->eigenvectors[i * this->dim + j] * this->righthand[j];
  }
  double *bIA	= new double [this->dim];
    
  for (unsigned i = 0; i < this->dim; i++)
  {
    bIA[i]	= 0.;
    if (i >= this->startIndex)
    {
      if (this->eigenvalues[i] != 0)
        bIA[i] = 1./this->eigenvalues[i] * bI[i];
      else
        bIA[i] = 0.;
    }
  }
    
  delete [] bI;
    
  for (unsigned i = 0; i < this->dim; i++)
  {
    this->solution[i]	= 0.;        
    for (unsigned j = 0; j < this->dim; j++)
      this->solution[i] += this->eigenvectors[j * this->dim + i] * bIA[j];
  }
  
  delete [] bIA;
}

/**
 * @brief Calculates the eigenvalues and eigenvectors of the coefficient matrix.
 * @throws SVDException If any error occur during the calculation of the linear equation system
 */
void SVDSolver::calculateEigenvalues()
{
  int error = 0;
  this->eigenvalues = new double[this->dim];
  this->eigenvectors = new double[this->dim * this->dim];
  this->offDiagonalElements = new double [this->dim];
  this->tred2();
  this->tql2(&error);
  
  if (error != 0)
  {
    std::ostringstream errMsg;
    errMsg << "Errors occured while calculating eigenvalues: " << error;
    throw SVDException(errMsg.str().c_str());
  }
}

/**
 * @brief Returns the solution of the linear equation system.
 * @throws SVDException If the solution has not been calculated yet
 */
double* SVDSolver::getSolution() const
{
  if (this->solution == NULL)
  {
    std::ostringstream errMsg;
    errMsg << "Solution has not been calculated!";
    throw SVDException(errMsg.str().c_str());
  }
  
  return this->solution;
}

/**
 * @brief Transforms a given vector into a representation with respect to the eigenvectors as basis.
 *
 * @param vec A vector with length equal to <code>dim</code> that should be transformed
 * @return The transformed vector with a length equal to the number of eigenvalues which were taken into account when solving the linear equation system
 */
double* SVDSolver::transform(double* vec) const
{
  double* tVec = new double[this->dim - this->startIndex];
  
  for (unsigned i = this->startIndex; i < this->dim; i++)
  {
    tVec[i - startIndex] = 0;
    for (unsigned j = 0; j < this->dim; j++)
    {
      tVec[i - startIndex] += this->getEigenvectorValue(j, i, false) * vec[j];
    }
  }
  
  return tVec;
}

/**
 * @brief Returns the effective dimension of the linear equation system that is the number of eigenvalues that have been taken into account while solving the linear equation system.
 * @return The effective dimension of the linear equation system
 */
unsigned SVDSolver::getEffectiveDim() const
{
  return (this->dim - this->startIndex);
}

#ifdef SVD_SYMMETRIC_STORAGE
unsigned SVDSolver::getIndex(unsigned i, unsigned j) const
{
  if (i > this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Row index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }
  if (j > this->dim)
  {
    std::ostringstream errMsg;
    errMsg << "Column index out of bounds!";
    throw SVDException(errMsg.str().c_str());
  }

  unsigned index = 0;
  if (i >= j)
    index = this->dim * j  - ((j + 1) * j)/2 + i;
  else
    index = this->dim * i  - ((i + 1) * i)/2 + j;
  return index;
}
#endif

/**
 * This value also represents the number of eigenvalues, the number of eigenvectors and the dimension of the eigenvectors.
 * @brief Returns the number of rows and columns of the coefficient matrix.
 *
 * @return The dimension on the coefficient matrix
 */ 
unsigned SVDSolver::getDim() const
{
  return this->dim;
}

/**
 * This function should only be called if the internal eigenvalue vector is not initialised.
 * @brief Sets the eigenvalues for this <code>SVDSolver</code>.
 *
 * @param evs A vector of eigenvalues with which the internal eigenvalue vector should be initialised
 * @throw SVDException If the eigenvalue vector is already initialised
 */
void SVDSolver::setEigenvalues(double* evs)
{
  if (this->eigenvalues != NULL)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvalue vector must not be initilised when setting new eigenvalue vector!";
    throw SVDException(errMsg.str().c_str());
  }
  this->eigenvalues = evs;
}

/**
 * This function should only be called if the internal eigenvector matrix is not initialised.
 * @brief Sets the vectors for this <code>SVDSolver</code>.
 *
 * @param evs A matrix containing the eigenvectors with which the internal eigenvector matrix should be initialised
 * @throw SVDException If the eigenvector matrix is already initialised
 */
void SVDSolver::setEigenvectors(double* evs)
{
  if (this->eigenvectors != NULL)
  {
    std::ostringstream errMsg;
    errMsg << "Eigenvectors must not be initialised when setting new eigenvalue vector!";
    throw SVDException(errMsg.str().c_str());
  }
  this->eigenvectors = evs;
}

/**
 * This function should only be called if the internal offdiagonal element vector is not initialised.
 * @brief Sets the offdiagonal elements for this <code>SVDSolver</code>.
 *
 * @param od A vector of offdiagonal elements with which the internal offdiagonal element vector should be initialised
 * @throw SVDException If the offdiagonal elemnet vector is already initialised
 */
void SVDSolver::setOffDiagonal(double* od)
{
  if (this->offDiagonalElements != NULL)
  {
    std::ostringstream errMsg;
    errMsg << "Offdiagonal must not be initialised when setting new offdiagonal!";
    throw SVDException(errMsg.str().c_str());
  }
  this->offDiagonalElements = od;
}

