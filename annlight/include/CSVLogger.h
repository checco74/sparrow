#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#ifndef _CSVLOGGER_H_
#define _CSVLOGGER_H_

using namespace std;

enum NumberFormat
{
  COMMA,
  POINT
};

/** Only data of standard C++ types can be logged. @brief Class for logging of program data into a CSV file. */
class CSVLogger
{
  public:
    CSVLogger();
    CSVLogger(const string* header, unsigned columnNumber);
    CSVLogger(string& filename, const string* header, unsigned columnNumber);
    CSVLogger(char seperator, string& filename, string* header, unsigned columnNumber);
    ~CSVLogger();

    // Data writing methods
    void addData(const char data);
    void addData(const string& data);
    void addData(const int data);
    void addData(const unsigned data);
    void addData(const double data);
    void addData(const bool data);
    void addData(const void* data);
    bool addLine(const string& line);
    void addSeperator();
    void addNewLine();

    // getter methods
    /**
     * @brief Returns the current number format.
     * @return The value of the <code>numberFormat</code> member varaible
     */
    const NumberFormat getNumberFormat() const {return numberFormat;}
    const string getNumberFormatAsString() const;
    /**
     * @brief Returns the string representation of a true bool value.
     * @return The value of the <code>boolTrueString</code> member variable
     */
    const string& getBoolTrueString() const {return boolTrueString;}
    /**
     * @brief Returns the string representation of a false bool value.
     * @return The value of the <code>boolFalseString</code> member variable
     */
    const string& getBoolFalseString() const {return boolFalseString;}
    /**
     * @brief Returns the name of the log file.
     * @return The value of the <code>logFileName</code> member variable
     */
    const string& getLogFileName() const {return logFileName;}
    /**
     * @brief Returns the precision for double values.
     * @return The value of the <code>precision</code> member variable
     */
    const unsigned getPrecision() const {return precision;}
    /**
     * @brief Determines whether the header was already written to the log.
     * @return The value of the <code>headerWritten</code> member variable
     */
    const bool getHeaderWritten() const {return headerWritten;}
    /**
     * @brief Determines whether the log file is open for writing.
     * @return <code>true</code> if the log file is open, <code>false</code> otherwise
     */
    bool isLogOpen() {return logFile.is_open();}

    // setter methods
    /**
     * @brief Sets the number format for this <code>CSVLogger</code>.
     * @param numberFormat The new number format
     */
    void setNumberFormat(const NumberFormat numberFormat) {this->numberFormat = numberFormat;}
    /**
     * @brief Sets a new string representing the true bool value.
     * @param boolTrueString The new string to represent the true bool value
     */
    void setBoolTrueString(const string& boolTrueString) {this->boolTrueString = boolTrueString;}
    /**
     * @brief Sets a new string representing the false bool value.
     * @param boolFalseString The new string to represent the false bool value
     */
    void setBoolFalseString(const string& boolFalseString) {this->boolFalseString = boolFalseString;}
    void setBoolStrings(const string& boolTrueString, const string& boolFalseString);
    /**
     * @brief Sets a new precision for double values.
     * @param precision The new precision to use
     */
    void setPrecision(const unsigned precision) {this->precision = precision; this->logFile.precision(precision);}

  protected:
    /** @brief The number of columns in the table. */
    unsigned columnNumber;

    /** @brief The index of the corrent column from the interval [0,columnNumber - 1]. */
    unsigned currentColumn;

    /** @brief Field seperator of the table. */
    char seperator;

    /** @brief The log file. */
    ofstream logFile;

    /** @brief The name of the log file. */
    string logFileName;

    /** @brief The number format for double values. */
    NumberFormat numberFormat;

    /** @brief Precision for double values. */
    unsigned precision;

    /** @brief String to write out for true bool values. */
    string boolTrueString;

    /** @brief String to write out for true bool values. */
    string boolFalseString;

    /** @brief The table header. */
    string* header;

    /** @brief Determines whether the header was already written to the log. */
    bool headerWritten;

    string changeNumberFormat(double number);
    void writeHeader();
};

#endif
