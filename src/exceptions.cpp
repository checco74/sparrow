#include "exceptions.h"

/**
 * @brief Standard constructor for a <code>GenericException</code>.
 */
GenericException::GenericException() throw()
{
  this->msg = new char[15];
  strcpy(this->msg,"unknown error");
}

/**
 * @brief Constructor for the <code>GenericException</code> class.
 *
 * @param msg The error message for this <code>GenericException</code>
 */
GenericException::GenericException(const char *msg) throw()
{
  this->msg = new char[strlen(msg)+1];
  strcpy(this->msg,msg);
}

/**
 * @brief Returns the message for this <code>GenericException</code>.
 *
 * @return The message for this <code>GenericException</code>
 */
const char* GenericException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Destructor for the <code>GenericException</code> class.
 */
GenericException::~GenericException() throw()
{
  delete[] this->msg;
}

/**
 * @brief Returns the type of this exception, in this case "GenericException".
 *
 * @return The current exception's type
 */
std::string GenericException::getType() const throw()
{
  return "GenericException";
}

// ==========================================================================================================================

/**
 * @brief Standard constructor for a <code>TransformationException</code>.
 */
TransformationException::TransformationException() throw() : GenericException()
{
}

/**
 * @brief Constructor for the <code>TransformationException</code> class.
 *
 * @param msg The error message for this <code>TransformationException</code>
 */
TransformationException::TransformationException(const char *msg) throw() : GenericException(msg)
{
}

/**
 * @brief Destructor for the <code>TransformationException</code> class.
 */
TransformationException::~TransformationException() throw()
{
}

/**
 * @brief Returns the message for this <code>TransformationException</code>.
 *
 * @return The message for this <code>TransformationException</code>
 */
const char* TransformationException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Returns the type of this exception, in this case "TransformationException".
 *
 * @return The current exception's type
 */
std::string TransformationException::getType() const throw()
{
  return "TransformationException";
}

// ==========================================================================================================================

/**
 * @brief Standard constructor for a <code>DatabaseException</code>.
 */
DatabaseException::DatabaseException() throw() : GenericException()
{
}

/**
 * @brief Constructor for the <code>DatabaseException</code> class.
 *
 * @param msg The error message for this <code>DatabaseException</code>
 */
DatabaseException::DatabaseException(const char *msg) throw() : GenericException(msg)
{
}

/**
 * @brief Destructor for the <code>DatabaseException</code> class.
 */
DatabaseException::~DatabaseException() throw()
{
}

/**
 * @brief Returns the message for this <code>DatabaseException</code>.
 *
 * @return The message for this <code>DatabaseException</code>
 */
const char* DatabaseException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Returns the type of this exception, in this case "DatabaseException".
 *
 * @return The current exception's type
 */
std::string DatabaseException::getType() const throw()
{
  return "TransformationException";
}

// ==========================================================================================================================

/**
 * @brief Standard constructor for a <code>PreprocessorException</code>.
 */
PreprocessorException::PreprocessorException() throw() : GenericException()
{
}

/**
 * @brief Constructor for the <code>PreprocessorException</code> class.
 *
 * @param msg The error message for this <code>PreprocessorException</code>
 */
PreprocessorException::PreprocessorException(const char *msg) throw() : GenericException(msg)
{
}

/**
 * @brief Destructor for the <code>PreprocessorException</code> class.
 */
PreprocessorException::~PreprocessorException() throw()
{
}

/**
 * @brief Returns the message for this <code>PreprocessorException</code>.
 *
 * @return The message for this <code>PreprocessorException</code>
 */
const char* PreprocessorException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Returns the type of this exception, in this case "PreprocessorException".
 *
 * @return The current exception's type
 */
std::string PreprocessorException::getType() const throw()
{
  return "PreprocessorException";
}

// ==========================================================================================================================

/**
 * @brief Standard constructor for a <code>ProfileReaderException</code>.
 */
ProfileReaderException::ProfileReaderException() throw() : GenericException()
{
}

/**
 * @brief Constructor for the <code>ProfileReaderException</code> class.
 *
 * @param msg The error message for this <code>ProfileReaderException</code>
 */
ProfileReaderException::ProfileReaderException(const char *msg) throw() : GenericException(msg)
{
}

/**
 * @brief Destructor for the <code>ProfileReaderException</code> class.
 */
ProfileReaderException::~ProfileReaderException() throw()
{
}

/**
 * @brief Returns the message for this <code>ProfileReaderException</code>.
 *
 * @return The message for this <code>ProfileReaderException</code>
 */
const char* ProfileReaderException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Returns the type of this exception, in this case "ProfileReaderException".
 *
 * @return The current exception's type
 */
std::string ProfileReaderException::getType() const throw()
{
  return "ProfileReaderException";
}
