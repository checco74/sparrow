#ifndef SVD_SOLVER_H_
#define SVD_SOLVER_H_

#define SVD_SYMMETRIC_STORAGE

typedef enum t_solving_type
{
  TAKE_BEST,
  TAKE_ABOVE_THRESHOLD,
  LEAVE_OUT_WORST,
  LEAVE_OUT_BELOW_THRESHOLD,
  TAKE_PERCENTAGE
};

/**
 * The class contains all data of the linear equation system as well as the eigenvalues and eigenvectors of the lefthand coefficient matrix.
 * @brief A solver for linear equation systems using principal component analysis.
 */
class SVDSolver
{
  public:
    SVDSolver(double* covMatrix, double* righthand, unsigned dim, bool copy = true);
    SVDSolver(double** covMatrix, double* righthand, unsigned dim, bool copy = true);
    ~SVDSolver();
    
    void solve(t_solving_type method, double value);
    void calculateEigenvalues();
    
    double getCovariance(unsigned i, unsigned j, bool start1 = true) const;
    void setCovariance(unsigned i, unsigned j, double value, bool start1 = true);

    double getOffDiagonalElement(unsigned i, bool start1 = true) const;
    void setOffDiagonalElement(unsigned i, double value, bool start1 = true);
    
    double getEigenvectorValue(unsigned i, unsigned j, bool start1 = true) const;
    double* getEigenvector(unsigned i, bool start1 = true) const;
    void setEigenvectorValue(unsigned i, unsigned j, double value, bool start1 = true);
    
    double getEigenvalue(unsigned i, bool start1 = true) const;
    void setEigenvalue(unsigned i, double value, bool start1 = true);
    
    void setEigenvalues(double* evs);
    void setEigenvectors(double* evs);
    void setOffDiagonal(double* od);
    
    double* getSolution() const;
    
    unsigned getEffectiveDim() const;
    
    double* transform(double* vec) const;
    
    unsigned getDim() const;
  
  protected:
    /** The lefthand coefficient matrix which in this case is a covariance matrix. */
    double* covMatrix;
    /** The offdiagonal elements calculated by the QL algorithm. */
    double* offDiagonalElements;
    /** The eigenvectors of <code>covMatrix</code>. */
    double* eigenvectors;
    /** The eigenvalues of <code>covMatrix</code>. */
    double* eigenvalues;
    /** The solution of the current linear equation system. */
    double* solution;
    /** The righthand side vector of the linear equation system. */
    double* righthand;

    /** The index of the first eigenvalue that should be taken into account. */
    unsigned startIndex;
    /** The dimension (number of rows) of <code>covMatrix</code>. */
    unsigned dim;
    
    void solve();
    void tred2();
    void tql2(int *ierr);
#ifdef SVD_SYMMETRIC_STORAGE
    /**
     * This function works with symmetric storage mode. It assumes that the matrix indices start with zero.
     * @brief Returns the index in the one-dimensional <code>covMatrix</code> array given the matrix indices of the requested entries.
     *
     * @param i The matrix row index
     * @param j The matrix column array
     * @return The index of the requested entry in <code>covMatrix</code>
     */
    unsigned getIndex(unsigned i, unsigned j) const;
#endif
};

#endif
