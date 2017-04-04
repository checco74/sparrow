#include "ANNException.h"

/**
 * @brief Standard constructor for a <code>ANNException</code>.
 */
ANNException::ANNException() throw()
{
  this->msg = new char[15];
  strcpy(this->msg,"unknown error");
}

/**
 * @brief Constructor for the <code>ANNException</code> class.
 *
 * @param msg The error message for this <code>ANNException</code>
 */
ANNException::ANNException(const char *msg) throw()
{
  this->msg = new char[strlen(msg)+1];
  strcpy(this->msg,msg);
}

/**
 * @brief Returns the message for this <code>ANNException</code>.
 *
 * @return The message for this <code>ANNException</code>
 */
const char* ANNException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Destructor for the <code>ANNException</code> class.
 */
ANNException::~ANNException() throw()
{
  delete[] this->msg;
}
