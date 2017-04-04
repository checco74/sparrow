#include "DiscriminantParameters.h"
#include <iostream>
#include <fstream>

double constantParameter = 0;

void readParam(char* filename, DiscriminantParameters* param, bool toLSM)
{
  std::ifstream inFile(filename, std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::cout << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  // reading of linear parameters is the same for both systems
  for (unsigned i = 0; i < param->getVecDim(); i++)
  {
    double value;
    inFile.read((char*)(&value), sizeof(double));
    param->set(i, value);
  }

  if (toLSM)
  {
    // read quadratic parameters for the optimizer
    for (unsigned i = 0; i < param->getMatDim(); i++)
    {
      for (unsigned j = i; j < param->getMatDim(); j++)
      {
        double value;
        inFile.read((char*)(&value), sizeof(double));
        param->set(i, j, value);
      }
    }
  }
  else
  {
    // read quadraticv parameters for LSM
    for (unsigned i = 0; i < param->getMatDim(); i++)
    {
      for (unsigned j = 0; j <= i; j++)
      {
        double value;
        inFile.read((char*)(&value), sizeof(double));
        param->set(i, j, value);
      }
    }

    inFile.read((char*)(&constantParameter), sizeof(double));
  }
  inFile.close();
}

void writeParam(char* filename, DiscriminantParameters* param, bool toLSM)
{
  std::ofstream outFile(filename, std::fstream::out | std::fstream::binary);
  if (!outFile.is_open())
  {
    std::cout << "Could not open file \"" << filename << "\" for reading!" << std::endl;
    exit(1);
  }
  // writing of linear parameters is the same for both systems
  for (unsigned i = 0; i < param->getVecDim(); i++)
  {
    double value = param->get(i);
    outFile.write((char*)(&value), sizeof(double));
  }

  if (toLSM)
  {
    // write quadratic parameters for LSM
    for (unsigned i = 0; i < param->getMatDim(); i++)
    {
      for (unsigned j = 0; j <= i; j++)
      {
        double value = param->get(i, j) * 2.0;
        outFile.write((char*)(&value), sizeof(double));
      }
    }

    outFile.write((char*)(&constantParameter), sizeof(double));
  }
  else
  {
    // write quadratic parameters for the optimizer
    for (unsigned i = 0; i < param->getMatDim(); i++)
    {
      for (unsigned j = i; j < param->getMatDim(); j++)
      {
        double value = param->get(i, j)/2.0;
        outFile.write((char*)(&value), sizeof(double));
      }
    }
  }
  outFile.close();
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cout << "Missing conversion mode!" << std::endl;
    exit(1);
  }
  if (strcmp(argv[1], "lsm") != 0 && strcmp(argv[2], "opt") != 0)
  {
    std::cout << "Unknown conversion mode: " << argv[1] << std::endl;
    std::cout << "Possible modes are \"lsm\" or \"opt\"" << std::endl;
    exit(1);
  }
  if (argc < 3)
  {
    std::cout << "Missing window size!" << std::endl;
    exit(1);
  }
  if (argc < 4)
  {
    std::cout << "Missing input file!" << std::endl;
    exit(1);
  }
  if (argc < 5)
  {
    std::cout << "Missing output file!" << std::endl;
    exit(1);
  }

  unsigned vecDim = 20 * atoi(argv[2]);
  bool toLSM = (strcmp(argv[1], "lsm") == 0);

  DiscriminantParameters* read = new DiscriminantParameters(vecDim);

  readParam(argv[3], read, toLSM);
  writeParam(argv[4], read, toLSM);

  delete read;
}
