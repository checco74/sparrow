#include "SVDSolver.h"
#include "SVDException.h"
#include <iostream>

int main(unsigned argc, char** argv)
{
#ifndef SVD_SYMMETRIC_STORAGE
  double* R = new double[9];
  R[0] = 1.5;
  R[1] = -0.5;
  R[2] = 0;
  R[3] = -0.5;
  R[4] = 1.5;
  R[5] = 0;
  R[6] = 0;
  R[7] = 0;
  R[8] = 3;
#else
  double* R = new double[6];
  R[0] = 1.5;
  R[1] = -0.5;
  R[2] = 0;
  R[3] = 1.5;
  R[4] = 0;
  R[5] = 3;
#endif
  double* b = new double[3];
  b[0] = 1;
  b[1] = 1;
  b[2] = 3;
  
  SVDSolver* solver = new SVDSolver(R, b, 3, false);
  
  solver->solve(TAKE_BEST, 3);
  
  std::cout << "Matrix:" << std::endl;
  for (unsigned i = 0; i < 3; i++)
  {
    for (unsigned j = 0; j < 3; j++)
    {
#ifndef SVD_SYMMETRIC_STORAGE
      std::cout << " " << R[i+j*3];
#else
      unsigned index = 0;
      if (i >= j)
        index = 3 * j  - ((j + 1) * j)/2 + i;
      else
        index = 3 * i  - ((i + 1) * i)/2 + j;
      std::cout << " " << R[index];
#endif
    }
      std::cout << std::endl;
  }
  std::cout << std::endl;
  
  std::cout << "Eigenvalues: " << std::endl;
  for (unsigned i = 1; i <= 3; i++)
    std::cout << solver->getEigenvalue(i) << std::endl;
  std::cout << std::endl;

  std::cout << "EigenVectors: " << std::endl;
  for (unsigned i = 1; i <= 3; i++)
  {
    for (unsigned j = 1; j <= 3; j++)
    {
      std::cout << " " << solver->getEigenvectorValue(i, j);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  std::cout << "Solution: " << std::endl;
  double* sol = solver->getSolution();
  for (unsigned i = 0; i < 3; i++)
    std::cout << sol[i] << std::endl;
  std::cout << std::endl;
  
  std::cout << "Eigenvector(2): " << std::endl;
  double* eVec = solver->getEigenvector(2);
  for (unsigned i = 0; i < 3; i++)
    std::cout << eVec[i] << std::endl;
  std::cout << std::endl;
  
  double vec [] = {1, 0, 3};
  double* tVec = solver->transform(vec);
  std::cout << "Transformed vector: " << std::endl;
  for (unsigned i = 0; i < 3; i++)
    std::cout << tVec[i] << std::endl;
  std::cout << std::endl;
  
  delete solver;
  delete [] eVec;
  delete [] tVec;
  
#ifndef SVD_SYMMETRIC_STORAGE
  R = new double[16];
  R[0] = 3.5;
  R[1] = 0;
  R[2] = 0;
  R[3] = 0.5;
  R[4] = 0;
  R[5] = 1;
  R[6] = 0;
  R[7] = 0;
  R[8] = 0;
  R[9] = 0;
  R[10] = 2;
  R[11] = 0;
  R[12] = 0.5;
  R[13] = 0;
  R[14] = 0;
  R[15] = 3.5;
#else
  R = new double[10];
  R[0] = 3.5;
  R[1] = 0;
  R[2] = 0;
  R[3] = 0.5;
  R[4] = 1;
  R[5] = 0;
  R[6] = 0;
  R[7] = 2;
  R[8] = 0;
  R[9] = 3.5;
#endif
  b = new double[4];
  b[0] = 5.5;
  b[1] = 2;
  b[2] = 0;
  b[3] = 14.5;
  
  solver = new SVDSolver(R, b, 4, false);
  
  solver->solve(TAKE_PERCENTAGE, 0.91);
  
  std::cout << "Matrix:" << std::endl;
  for (unsigned i = 0; i < 4; i++)
  {
    for (unsigned j = 0; j < 4; j++)
    {
#ifndef SVD_SYMMETRIC_STORAGE
      std::cout << " " << R[i+j*4];
#else
      unsigned index = 0;
      if (i >= j)
        index = 4 * j  - ((j + 1) * j)/2 + i;
      else
        index = 4 * i  - ((i + 1) * i)/2 + j;
      std::cout << " " << R[index];
#endif
    }
      std::cout << std::endl;
  }
  std::cout << std::endl;
  
  std::cout << "Eigenvalues: " << std::endl;
  for (unsigned i = 1; i <= solver->getEffectiveDim(); i++)
    std::cout << solver->getEigenvalue(i) << std::endl;
  std::cout << std::endl;

  std::cout << "EigenVectors: " << std::endl;
  for (unsigned i = 1; i <= 4; i++)
  {
    for (unsigned j = 1; j <= 4; j++)
    {
      std::cout << " " << solver->getEigenvectorValue(i, j);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  std::cout << "Solution: " << std::endl;
  sol = solver->getSolution();
  for (unsigned i = 0; i < 4; i++)
    std::cout << sol[i] << std::endl;
  std::cout << std::endl;
  
  std::cout << "Eigenvector(3): " << std::endl;
  eVec = solver->getEigenvector(3);
  for (unsigned i = 0; i < 4; i++)
    std::cout << eVec[i] << std::endl;
  std::cout << std::endl;
  
  double vec2 [] = {1, 0, 3, 2};
  tVec = solver->transform(vec2);
  std::cout << "Transformed vector: " << std::endl;
  for (unsigned i = 0; i < solver->getEffectiveDim(); i++)
    std::cout << tVec[i] << std::endl;
  std::cout << std::endl;
  
  delete solver;
  delete [] eVec;
  delete [] tVec;
}

