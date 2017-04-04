#include "IOException.h"

/**
 * @brief Standard constructor for a <code>IOException</code>.
 */
IOException::IOException() throw()
{
  this->msg = new char[15];
  strcpy(this->msg,"unknown error");
}

/**
 * @brief Constructor for the <code>IOException</code> class.
 *
 * @param msg The error message for this <code>IOException</code>
 */
IOException::IOException(const char *msg) throw()
{
  this->msg = new char[strlen(msg)+1];
  strcpy(this->msg,msg);
}

/**
 * @brief Returns the message for this <code>IOException</code>.
 *
 * @return The message for this <code>IOException</code>
 */
const char* IOException::what() const throw()
{
  return this->msg;
}

/**
 * @brief Destructor for the <code>IOException</code> class.
 */
IOException::~IOException() throw()
{
  delete[] this->msg;
}
