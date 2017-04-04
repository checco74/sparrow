#include "SVDException.h"

/**
 * @brief Standard constructor for a <code>SVDException</code>.
 */
SVDException::SVDException() throw()
{
  this->msg = new char[15];
  strcpy(this->msg,"unknown error");
}

/**
 * @brief Constructor for the <code>SVDException</code> class.
 *
 * @param msg The error message for this <code>SVDException</code>
 */
SVDException::SVDException(const char *msg) throw()
{
  this->msg = new char[strlen(msg)+1];
  strcpy(this->msg,msg);
}

/**
 * @brief Returns the message for this <code>SVDException</code>.
 *
 * @return The message for this <code>SVDException</code>
 */
const char* SVDException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Destructor for the <code>ParserException</code> class.
 */
SVDException::~SVDException() throw()
{
  delete[] this->msg;
}
