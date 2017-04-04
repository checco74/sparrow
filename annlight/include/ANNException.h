#include <exception>
#include <stdexcept>

#ifndef _ANNEXCEPTION_H_
#define _ANNEXCEPTION_H_

/**
 * @brief Class for ANN exceptions.
 */
class ANNException : public std::exception
{
  public:
    ANNException() throw();
    ANNException(const char *msg) throw();

    ~ANNException() throw();

    const char* what() const throw();
  protected:
    /** @brief The exception's error message. */
    char   *msg;
};

#endif
