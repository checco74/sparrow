#include <exception>
#include <stdexcept>

#ifndef _IOEXCEPTION_H_
#define _IOEXCEPTION_H_

/**
 * @brief Class for I/O exceptions.
 */
class IOException : public std::exception
{
  public:
    IOException() throw();
    IOException(const char *msg) throw();

    ~IOException() throw();

    const char* what() const throw();
  protected:
    /** @brief The exception's error message. */
    char   *msg;
};

#endif
