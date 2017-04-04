#include "DiscriminantParameters.h"

int main(int argc, char** argv)
{
  DiscriminantParameters* param1 = new DiscriminantParameters(3);
  DiscriminantParameters* param2 = new DiscriminantParameters(3, true);

  std::cout << "PARAM1:" << std::endl;
  std::cout << *param1 << std::endl;
  std::cout << "PARAM2:" << std::endl;
  std::cout << *param2 << std::endl;
  DiscriminantParameters param3 = ((*param2) * 10.0);
  std::cout << "Multiplikation:" << std::endl;
  std::cout << param3 << std::endl;
  std::cout << "PARAM2:" << std::endl;
  std::cout << *param2 << std::endl;
  *param2 *= 10.0;
  std::cout << "Multiplikation 2:" << std::endl;
  std::cout << *param2 << std::endl;
  std::cout << "Addition:" << std::endl;
  std::cout << (*param2 + param3) << std::endl;
  *param2 += param3;
  std::cout << "Addition 2:" << std::endl;
  std::cout << *param2 << std::endl;

  delete param1;
  delete param2;
}
