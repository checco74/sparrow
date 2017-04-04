#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "ERROR: Table name is missing!" << std::endl;
    exit(1);
  }
  if (argc < 3)
  {
    std::cerr << "ERROR: Starting index is missing!" << std::endl;
    exit(1);
  }
  if (argc < 4)
  {
    std::cerr << "ERROR: Stopping index is missing!" << std::endl;
    exit(1);
  }

  std::string tabName = argv[1];
  int start = atoi(argv[2]);
  int stop = atoi(argv[3]);
  int setSize = stop - start + 1;

  std::ofstream outFile("dumbset.dat", std::fstream::out | std::fstream::trunc);
  if (!outFile.is_open())
  {
    std::cerr << "ERROR: Could not open output file for writing!" << std::endl;
    exit(1);
  }
  outFile << tabName << std::endl;
  outFile << setSize << " 0 0" << std::endl;
  int setPos = 0;
  for (int i = start; i <= stop; i++)
  {
    outFile << setPos << " " << i << " 0" << std::endl;
    setPos++;
  }
  outFile.close();
}
