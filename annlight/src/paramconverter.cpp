#include <iostream>
#include <fstream>

bool toBinary(std::string& str)
{
  if (str == "binary")
    return true;
  else if (str == "plain")
    return false;
  else
  {
    std::cerr << "ERROR: Unknown coversion mode string: " << str << std::endl;
    exit(1);
  }
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "ERROR: Missing conversion mode!" << std::endl;
    exit(1);
  }
  std::string mode = argv[1];
  bool toBin = toBinary(mode);
  if (argc < 3)
  {
    std::cerr << "ERROR: Missing input file name!" << std::endl;
    exit(1);
  }
  std::string filename = argv[2];

  std::ifstream inFile;
  std::ofstream outFile;
  if (toBin)
  {
    inFile.open(filename.c_str(), std::fstream::in | std::fstream::binary);
    outFile.open("temp.prm", std::fstream::out | std::fstream::trunc);
  }
  else
  {
    inFile.open(filename.c_str(), std::fstream::in);
    outFile.open(filename.c_str(), std::fstream::out | std::fstream::trunc | std::fstream::binary);
  }
  if (!inFile.is_open())
  {
    std::cerr << "Parameter file \"" << filename << "\" could not be opened for reading!" << std::endl;
    exit(1);
  }
  if (!outFile.is_open())
  {
    std::cerr << "Output file could not be opened for writing!" << std::endl;
    exit(1);
  }

  if (toBin)
  {
    std::string line;
    while (getline(inFile, line, '\n'))
    {
      double value = atof(line.c_str());
      outFile.write((char*)(&value), sizeof(value));
    }
  }
  else
  {    
    while (!inFile.eof())
    {
      double value;
      inFile.read((char*)(&value), sizeof(value));
      outFile << value << std::endl;
    }
  }
  
  inFile.close();
  outFile.close();
}

