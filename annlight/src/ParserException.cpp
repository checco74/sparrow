#include "ParserException.h"

/**
 * @brief Standard constructor for a <code>ParserException</code>.
 */
ParserException::ParserException() throw()
{
  this->msg = new char[15];
  strcpy(this->msg,"unknown error");
}

/**
 * @brief Constructor for the <code>ParserException</code> class.
 *
 * @param msg The error message for this <code>ParserException</code>
 */
ParserException::ParserException(const char *msg) throw()
{
  this->msg = new char[strlen(msg)+1];
  strcpy(this->msg,msg);
}

/**
 * @brief Returns the message for this <code>ParserException</code>.
 *
 * @return The message for this <code>ParserException</code>
 */
const char* ParserException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Destructor for the <code>ParserException</code> class.
 */
ParserException::~ParserException() throw()
{
  delete[] this->msg;
}
