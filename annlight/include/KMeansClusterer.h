#include <exception>
#include <stdexcept>
#include <string>
#include <list>
#include <sstream>
#include <fstream>

#ifndef _K_MEANS_CLUSTERER_H_
#define _K_MEANS_CLUSTERER_H_

class KMeansException : public std::exception
{
  public:
    KMeansException() throw();
    KMeansException(const char *msg) throw();

    virtual ~KMeansException() throw();

    const char* what() const throw();
  protected:
    char   *msg;
};

class KMeansClusterer
{
  public:
    KMeansClusterer();
    KMeansClusterer(unsigned numberOfClusters, unsigned dimension);
    ~KMeansClusterer();
    
    void writePrototypes(const std::string& filename) const;
    void readPrototypes(const std::string& filename);
    
    void learnPrototypes(double** data, unsigned* initialAssignment, unsigned setSize);
    
    double getDistance(double* vector, unsigned clusterID) const;
    double* getDistances(double* vector);
    unsigned getAssignment(double* vector);
    
    unsigned getDimension() const;
    unsigned getNumberOfClusters() const;
    
    void writeProtocol(std::ostream& s, double** dataset = NULL, unsigned setSize = 0);
    
  protected:
    void calculatePrototypes(double** data, unsigned* assignment, unsigned setSize);
    
    unsigned numberOfClusters;
    unsigned dimension;
    double** prototypes;
};

#endif
