#include "KMeansClusterer.h"
#include <cmath>
#include <iostream>

KMeansException::KMeansException() throw()
{
  this->msg = new char[15];
  strcpy(this->msg,"unknown error");
}

KMeansException::KMeansException(const char *msg) throw()
{
  this->msg = new char[strlen(msg)+1];
  strcpy(this->msg,msg);
}

const char* KMeansException::what() const throw()
{
  return this->msg;
}

KMeansException::~KMeansException() throw()
{
  delete[] this->msg;
}


KMeansClusterer::KMeansClusterer()
{
  this->numberOfClusters = 0;
  this->dimension = 0;

  this->prototypes = NULL;
}

KMeansClusterer::KMeansClusterer(unsigned numberOfClusters, unsigned dimension)
{
  this->numberOfClusters = numberOfClusters;
  this->dimension = dimension;

  this->prototypes = new double*[this->numberOfClusters];
  for (unsigned i = 0; i < this->numberOfClusters; i++)
    this->prototypes[i] = new double[this->dimension];  
}

KMeansClusterer::~KMeansClusterer()
{
  if (this->prototypes != NULL)
  {
    for (unsigned i = 0; i < numberOfClusters; i++)
      delete [] this->prototypes[i];
    delete [] this->prototypes;
  }
}
    
void KMeansClusterer::writePrototypes(const std::string& filename) const
{
  if (this->prototypes == NULL)
  {
    std::ostringstream errorMessage;
    errorMessage << "Cannot write uninitialized prototypes to file!" << std::endl;
    throw KMeansException(errorMessage.str().c_str());
  }
  
  std::ofstream outFile(filename.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);
  if (!outFile.is_open())
  {
    std::ostringstream errorMessage;
    errorMessage << "Parameter file \"" << filename << "\" could not be opened for writing!";
    throw KMeansException(errorMessage.str().c_str());
  }
  
  outFile.write((char*)(&(this->numberOfClusters)), sizeof(this->numberOfClusters));
  outFile.write((char*)(&(this->dimension)), sizeof(this->dimension));
  
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++)
    for (unsigned entry = 0; entry < this->dimension; entry++)
      outFile.write((char*)(&(this->prototypes[clusterID][entry])), sizeof(this->prototypes[clusterID][entry]));
  
  outFile.close();
}

void KMeansClusterer::readPrototypes(const std::string& filename)
{
  if (this->prototypes != NULL)
  {
    for (unsigned i = 0; i < numberOfClusters; i++)
      delete [] this->prototypes[i];
  }
  
  std::ifstream inFile(filename.c_str(), std::fstream::in | std::fstream::binary);
  if (!inFile.is_open())
  {
    std::ostringstream errorMessage;
    errorMessage << "Parameter file \"" << filename << "\" could not be opened for reading!";
    throw KMeansException(errorMessage.str().c_str());
  }
  
  inFile.read((char*)(&(this->numberOfClusters)), sizeof(this->numberOfClusters));
  inFile.read((char*)(&(this->dimension)), sizeof(this->dimension));
  
  this->prototypes = new double*[this->numberOfClusters];
  for (unsigned i = 0; i < this->numberOfClusters; i++)
    this->prototypes[i] = new double[this->dimension];
  
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++)
    for (unsigned entry = 0; entry < this->dimension; entry++)
      inFile.read((char*)(&(this->prototypes[clusterID][entry])), sizeof(this->prototypes[clusterID][entry]));
  
  inFile.close();
}
    
void KMeansClusterer::learnPrototypes(double** data, unsigned* initialAssignment, unsigned setSize)
{
  if (this->prototypes == NULL)
  {
    std::ostringstream errorMessage;
    errorMessage << "Prototypes have not been initialized!" << std::endl;
    throw KMeansException(errorMessage.str().c_str());
  }

  unsigned* assignment = new unsigned[setSize];
  for (unsigned i = 0; i < setSize; i++)
  {
    if (initialAssignment[i] >= this->numberOfClusters)
    {
      std::ostringstream errorMessage;
      errorMessage << "Wrong initial assignment (" << i << "): " << initialAssignment[i] << std::endl;
      throw KMeansException(errorMessage.str().c_str());
    }
    assignment[i] = initialAssignment[i];
  }
  
  this->calculatePrototypes(data, assignment, setSize);
  
  bool redistributed = true;
  unsigned runs = 0;
  do
  {
    redistributed = false;
    for (unsigned sample = 0; sample < setSize; sample++)
    {
      double minDistance = this->getDistance(data[sample], 0);
      unsigned minCluster = 0;
      for (unsigned clusterID = 1; clusterID < this->numberOfClusters; clusterID++)
      {
        double distance = this->getDistance(data[sample], clusterID);
        if (distance < minDistance)
        {
          minDistance = distance;
          minCluster = clusterID;
        }
      }
      if (minCluster != assignment[sample])
      {
        assignment[sample] = minCluster;
        redistributed = true;
      }
    }
    
    this->calculatePrototypes(data, assignment, setSize);
    runs++;
    
  } while (redistributed);
  
  std::cout << "Iterations: " << runs << std::endl;
  
  delete [] assignment;
}

void KMeansClusterer::calculatePrototypes(double** data, unsigned* assignment, unsigned sampleSize)
{
  if (this->prototypes == NULL)
  {
    std::ostringstream errorMessage;
    errorMessage << "Prototypes have not been initialized!" << std::endl;
    throw KMeansException(errorMessage.str().c_str());
  }
  
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++)
    for (unsigned entry = 0; entry < this->dimension; entry++)
      this->prototypes[clusterID][entry] = 0;
  
  unsigned* clusterMembers = new unsigned[this->numberOfClusters];
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++) clusterMembers[clusterID] = 0;

  for (unsigned sample = 0; sample < sampleSize; sample++)
  {
    clusterMembers[assignment[sample]]++;
    for (unsigned entry = 0; entry < this->dimension; entry++)
      this->prototypes[assignment[sample]][entry] += data[sample][entry];
  }
  
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++)
    for (unsigned entry = 0; entry < this->dimension; entry++)
      if (clusterMembers[clusterID] != 0) this->prototypes[clusterID][entry] /= ((double)clusterMembers[clusterID]);
  
  delete [] clusterMembers;
}
    
double KMeansClusterer::getDistance(double* vector, unsigned clusterID) const
{
  if (this->prototypes == NULL)
  {
    std::ostringstream errorMessage;
    errorMessage << "Prototypes have not been initialized!" << std::endl;
    throw KMeansException(errorMessage.str().c_str());
  }
  
  double sumOfSquares = 0;
  
  for (unsigned i = 0; i < this->dimension; i++)
  {
    double difference = vector[i] - this->prototypes[clusterID][i];
    sumOfSquares += difference * difference; 
  }
  
  return std::sqrt(sumOfSquares);
}

double* KMeansClusterer::getDistances(double* vector)
{
  if (this->prototypes == NULL)
  {
    std::ostringstream errorMessage;
    errorMessage << "Prototypes have not been initialized!" << std::endl;
    throw KMeansException(errorMessage.str().c_str());
  }
  
  double* distances = new double[this->numberOfClusters];
  
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++)
    distances[clusterID] = this->getDistance(vector, clusterID);
  
  return distances;
}

unsigned KMeansClusterer::getDimension() const
{
  return this->dimension;
}

unsigned KMeansClusterer::getNumberOfClusters() const
{
  return this->numberOfClusters;
}

void KMeansClusterer::writeProtocol(std::ostream& s, double** dataset, unsigned setSize)
{
  s << "Number of clusters: " << this->numberOfClusters << std::endl;
  s << "Input dimension: " << this->dimension << std::endl;
  s << "Prototypes:" << std::endl;
  for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++)
  {
    s << clusterID << ":";
    for (unsigned entry = 0; entry < this->dimension; entry++)
      s << " " << this->prototypes[clusterID][entry];
    s << std::endl;
  }
  if (dataset != NULL)
  {
    unsigned* clusterMembers = new unsigned[this->numberOfClusters];
    for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++) clusterMembers[clusterID] = 0;
    s << std::endl << "Data assignment:" << std::endl;
    for (unsigned sample = 0; sample < setSize; sample++)
    {
      double minDistance = this->getDistance(dataset[sample], 0);
      unsigned minCluster = 0;
      for (unsigned clusterID = 1; clusterID < this->numberOfClusters; clusterID++)
      {
        double distance = this->getDistance(dataset[sample], clusterID);
        if (distance < minDistance)
        {
          minDistance = distance;
          minCluster = clusterID;
        }
      }
      s << "Sample " << sample << ":";
      for (unsigned entry = 0; entry < this->dimension; entry++) s << " " << dataset[sample][entry];
      s << " -> " << minCluster << std::endl;
      clusterMembers[minCluster]++;
    }
    
    s << std::endl << "Cluster member counts:" << std::endl;
    for (unsigned clusterID = 0; clusterID < this->numberOfClusters; clusterID++) s << clusterID << ": " << clusterMembers[clusterID] << std::endl;
    delete [] clusterMembers;
  }
}

unsigned KMeansClusterer::getAssignment(double* vector)
{
  double* distances = this->getDistances(vector);
  
  unsigned min = 0;
  double minDist = distances[0]; 
  for (unsigned clusterID = 1; clusterID < this->numberOfClusters; clusterID++)
  {
    if (minDist > distances[clusterID])
    {
      min = clusterID;
      minDist = distances[clusterID];
    }
  }
  
  delete [] distances;
  return min;
}

