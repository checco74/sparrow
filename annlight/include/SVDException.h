#include <exception>
#include <stdexcept>

#ifndef _SVD_EXCEPTION_H_
#define _SVD_EXCEPTION_H_

/**
 * @brief Class for <code>SVDSolver</code> exceptions.
 */
class SVDException : public std::exception
{
  public:
    SVDException() throw();
    SVDException(const char *msg) throw();

    ~SVDException() throw();

    const char* what() const throw();
  protected:
    /** @brief The exception's error message. */
    char   *msg;
};

#endif
