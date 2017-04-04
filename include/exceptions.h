#include <exception>
#include <stdexcept>
#include <string>

#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

/**
 * All exceptions that might occur are derived from this class.
 * @brief A generic exception.
 */
class GenericException : public std::exception
{
  public:
    GenericException() throw();
    GenericException(const char *msg) throw();

    virtual ~GenericException() throw();

    const char* what() const throw();
    virtual std::string getType() const throw();
  protected:
    /** @brief The exception's error message. */
    char   *msg;
};

/**
 * @brief This exception is thrown if any problems occur when using the <code>Transformations</code> class' methods.
 */
class TransformationException : public GenericException
{
  public:
    TransformationException() throw();
    TransformationException(const char *msg) throw();

    ~TransformationException() throw();

    const char* what() const throw();
    std::string getType() const throw();
};

/**
 * @brief This exception is thrown if any problems occur while accessing the database.
 */
class DatabaseException : public GenericException
{
  public:
    DatabaseException() throw();
    DatabaseException(const char *msg) throw();

    ~DatabaseException() throw();

    const char* what() const throw();
    std::string getType() const throw();
};

/**
 * @brief This exception is thrown if any problems occur while preprocessing a SSE string.
 */
class PreprocessorException : public GenericException
{
  public:
    PreprocessorException() throw();
    PreprocessorException(const char *msg) throw();

    ~PreprocessorException() throw();

    const char* what() const throw();
    std::string getType() const throw();
};

/**
 * @brief This exception is thrown if any problems occur while reading profile data from a MTX file.
 */
class ProfileReaderException : public GenericException
{
  public:
    ProfileReaderException() throw();
    ProfileReaderException(const char *msg) throw();

    ~ProfileReaderException() throw();

    const char* what() const throw();
    std::string getType() const throw();
};
#endif
