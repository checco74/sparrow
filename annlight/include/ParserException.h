#include <exception>
#include <stdexcept>

#ifndef _PARSER_EXCEPTION_H_
#define _PARSER_EXCEPTION_H_

/**
 * @brief Class for <code>ParameterParser</code> exceptions.
 */
class ParserException : public std::exception
{
  public:
    ParserException() throw();
    ParserException(const char *msg) throw();

    ~ParserException() throw();

    const char* what() const throw();
  protected:
    /** @brief The exception's error message. */
    char   *msg;
};

#endif
