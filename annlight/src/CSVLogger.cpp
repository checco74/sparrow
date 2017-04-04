#include "CSVLogger.h"

#define STANDARD_SEPERATOR '\t'
#define STANDARD_FILENAME "data.csv"
#define STANDARD_COLUMN_NUMBER 2
#define STANDARD_HEADER1 "x"
#define STANDARD_HEADER2 "y"
#define STANDARD_TRUE_STRING "true"
#define STANDARD_FALSE_STRING "false"
#define STANDARD_PRECISION 4

/**
 * Creates a new CSVLogger with standard settings.
 * @brief Standard constructor.
 */
CSVLogger::CSVLogger()
{
  this->seperator = STANDARD_SEPERATOR;
  this->columnNumber = STANDARD_COLUMN_NUMBER;

  this->currentColumn = 0;
  this->headerWritten = false;
  this->boolTrueString = STANDARD_TRUE_STRING;
  this->boolFalseString = STANDARD_FALSE_STRING;
  this->logFileName = STANDARD_FILENAME;
  this->logFile.open(this->logFileName.c_str(), fstream::out | fstream::trunc);
  this->logFile.precision(STANDARD_PRECISION);
  this->precision = STANDARD_PRECISION;
  this->logFile.setf(ios::fixed);
  this->numberFormat = COMMA;

  this->header = new string[this->columnNumber];
  this->header[0] = STANDARD_HEADER1;
  this->header[1] = STANDARD_HEADER2;
}

/**
 * Creates a new CSVLogger with standard seperator and standard filename.
 * @brief Constructor.
 *
 * @param header An array containing the header for each column
 * @param columnNumber Number of columns for the table
 */
CSVLogger::CSVLogger(const string* header, unsigned columnNumber)
{
  this->seperator = STANDARD_SEPERATOR;
  this->columnNumber = columnNumber;

  this->currentColumn = 0;
  this->headerWritten = false;
  this->boolTrueString = STANDARD_TRUE_STRING;
  this->boolFalseString = STANDARD_FALSE_STRING;
  this->logFileName = STANDARD_FILENAME;
  this->logFile.open(this->logFileName.c_str(), fstream::out | fstream::trunc);
  this->logFile.precision(STANDARD_PRECISION);
  this->precision = STANDARD_PRECISION;
  this->logFile.setf(ios::fixed);
  this->numberFormat = COMMA;

  this->header = new string[columnNumber];
  for (unsigned i = 0; i < columnNumber; i++)
  {
    this->header[i] = header[i];
  }
}

/**
 * Creates a new CSVLogger with standard seperator.
 * @brief Constructor.
 *
 * @param filename The name of the logfile
 * @param header An array containing the header for each column
 * @param columnNumber Number of columns for the table
 */
CSVLogger::CSVLogger(string& filename, const string* header, unsigned columnNumber)
{
  this->seperator = STANDARD_SEPERATOR;
  this->columnNumber = columnNumber;

  this->currentColumn = 0;
  this->headerWritten = false;
  this->boolTrueString = STANDARD_TRUE_STRING;
  this->boolFalseString = STANDARD_FALSE_STRING;
  this->logFileName = filename;
  this->logFile.open(this->logFileName.c_str(), fstream::out | fstream::trunc);
  this->logFile.precision(STANDARD_PRECISION);
  this->precision = STANDARD_PRECISION;
  this->logFile.setf(ios::fixed);
  this->numberFormat = COMMA;

  this->header = new string[columnNumber];
  for (unsigned i = 0; i < columnNumber; i++)
  {
    this->header[i] = header[i];
  }
}

/**
 * @brief Constructor.
 *
 * @param seperator The seperator for table columns
 * @param filename The name of the logfile
 * @param header An array containing the header for each column
 * @param columnNumber Number of columns for the table
 */
CSVLogger::CSVLogger(char seperator, string& filename, string* header, unsigned columnNumber)
{
  this->seperator = seperator;
  this->columnNumber = columnNumber;

  this->currentColumn = 0;
  this->headerWritten = false;
  this->boolTrueString = STANDARD_TRUE_STRING;
  this->boolFalseString = STANDARD_FALSE_STRING;
  this->logFileName = filename;
  this->logFile.open(filename.c_str(), fstream::out | fstream::trunc);
  this->logFile.precision(STANDARD_PRECISION);
  this->precision = STANDARD_PRECISION;
  this->logFile.setf(ios::fixed);
  this->numberFormat = COMMA;

  this->header = new string[columnNumber];
  for (unsigned i = 0; i < columnNumber; i++)
  {
    this->header[i] = header[i];
  }
}

/**
 * @brief Standard destructor.
 */
CSVLogger::~CSVLogger()
{
  if(logFile.is_open()) logFile.close();
  delete [] header;
}

/**
 * @brief Adds new data to the log. Adds seperators and newline if necessary.
 * 
 * @param data The data to add
 */
void CSVLogger::addData(const char data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    logFile << data;
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}

/**
 * @brief Adds new data to the log. Adds seperators and newline if necessary.
 *
 * @param data The data to add
 */
void CSVLogger::addData(const string& data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    logFile << data;
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}

/**
 * @brief Adds new data to the log. Adds seperators and newline if necessary.
 *
 * @param data The data to add
 */
void CSVLogger::addData(const int data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    logFile << data;
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}

/**
 * @brief Adds new data to the log. Adds seperators and newline if necessary.
 *
 * @param data The data to add
 */
void CSVLogger::addData(const unsigned data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    logFile << data;
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}


/**
 * Adds seperators and newline if necessary.
 * @brief Adds new data to the log.
 * 
 * @param data The data to add
 */
void CSVLogger::addData(const double data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    logFile << changeNumberFormat(data);
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}

/**
 * Adds seperators and newline if necessary.
 * @brief Adds new data to the log.
 * 
 * @param data The data to add
 */
void CSVLogger::addData(const bool data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    if(data)
      logFile << boolTrueString;
    else
      logFile << boolFalseString;
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}

/**
 * Adds seperators and newline if necessary.
 * @brief Adds new data to the log.
 * 
 * @param data The data to add
 */
void CSVLogger::addData(const void* data)
{
  if (logFile.is_open()) {
    if(!headerWritten) writeHeader();
    if(data)
      logFile << data;
    else
      logFile << "NULL";
    if (currentColumn == columnNumber - 1) addNewLine();
    currentColumn = (currentColumn + 1) % columnNumber;
    if (currentColumn != 0) addSeperator();
  }
}

/**
 * @brief Adds a seperator to the log.
 */
void CSVLogger::addSeperator()
{
  if (logFile.is_open()) logFile << seperator;
}

/**
 * If the current column is not the last column of the table seperators are added until
 * the last column is reched.
 * @brief Adds a newline to the log.
 */
void CSVLogger::addNewLine()
{
  if (logFile.is_open())
  {
    for (unsigned i = currentColumn; i < columnNumber - 1; i++)
    {
      addSeperator();
    }
    logFile << endl;
  }
}

/**
 * @brief Returns the number Format as string.
 * 
 * @return The current number format as string
 */
const string CSVLogger::getNumberFormatAsString() const
{
  string numberFormatString;
  switch (numberFormat) 
  {
    case COMMA: numberFormatString = "COMMA";
    case POINT: numberFormatString = "POINT";
  }
  return numberFormatString;
}

/**
 * @brief Sets the string values which should be written to the log, if the value to add is a boolean.
 *
 * @param boolTrueString The string which should be written out for true boolean values
 * @param boolFalseString The string which should be written out for false boolean values
 */
void CSVLogger::setBoolStrings(const string& boolTrueString, const string& boolFalseString)
{
  this->boolTrueString = boolTrueString;
  this->boolFalseString = boolFalseString;
}

/**
 * @brief Changes the format of the given number and returns a string containing the number in this format.
 *
 * @param number The number for which the number format should be changed
 * @return A string with the newly formatted number
 */
string CSVLogger::changeNumberFormat(double number)
{
  ostringstream dataStringstream;
  dataStringstream.precision(this->precision);
  dataStringstream.setf(ios::fixed);
  dataStringstream << number;
  string dataString;
  if (this->numberFormat == COMMA)
  {
    std::string::size_type startIndex = dataStringstream.str().find_first_of('.');
    if (startIndex != string::npos)
      dataString = dataStringstream.str().replace(startIndex, 1 , ",");
    else
      dataString = dataStringstream.str();
  }
  else
  {
    dataString = dataStringstream.str();
  }

  return dataString;
}

/**
 * @brief Writes the table header to the log if it was not written yet.
 */
void CSVLogger::writeHeader()
{
  if (this->header && logFile.is_open())
  {
    this->headerWritten = true;
    addData(header[0]);
    for (unsigned i = 1; i < columnNumber; i++)
    {
      addData(header[i]);
    }
  }
}

/**
 * This is done only if the header was not written yet, which means that lines can
 * be added before the table only.
 * @brief Adds a defined line to the log.
 *
 * @param line The line to add to the log
 * @return true if the line could be added to the log, false otherwise
 */
bool CSVLogger::addLine(const string& line) 
{
  if (!this->logFile.is_open() || this->headerWritten)
    return false;
  else
  {
    logFile << line;
    if (line.find('\n') == std::string::npos) logFile << endl;
    return true;
  }
}

