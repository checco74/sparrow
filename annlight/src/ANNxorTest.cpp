#include "ANN.h"
#include "time.h"
#include <iostream>

int main(int argc, char** argv)
{
  ANN* myANN = new ANN(2, 2, 1, true);

  unsigned samples = 4;
  double** learnData = new double*[4];
  double** expected = new double*[4];
  // (0,0)
  learnData[0] = new double[2];
  learnData[0][0] = 0;
  learnData[0][1] = 0;
  expected[0] = new double[1];
  expected[0][0] = -1;
  // (0,1)
  learnData[1] = new double[2];
  learnData[1][0] = 0;
  learnData[1][1] = 1;
  expected[1] = new double[1];
  expected[1][0] = 1;
  // (1,0)
  learnData[2] = new double[2];
  learnData[2][0] = 1;
  learnData[2][1] = 0;
  expected[2] = new double[1];
  expected[2][0] = 1;
  // (1,1)
  learnData[3] = new double[2];
  learnData[3][0] = 1;
  learnData[3][1] = 1;
  expected[3] = new double[1];
  expected[3][0] = -1;

  time_t t1 = time(NULL);
  //myANN->setFuncType(TANH);
  myANN->useLiteratureValues();
  myANN->setShuffleSet(false);
  myANN->initLearnData(learnData, expected, samples, true);
  myANN->learn(100000, 0);
  time_t t2 = time(NULL);

  std::cout << "Time: " << (t2-t1) << " seconds" << std::endl;

  for (unsigned i = 0; i < 4; i++)
  {
    delete [] expected[i];
    delete [] learnData[i];
  }
  delete [] expected;
  delete [] learnData;

  delete myANN;
}
