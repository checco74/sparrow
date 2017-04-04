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

/** Class for logging of program data into a CSV file. Only data of standard C++ types can be logged. */
class CSVLogger
{
	public:
	CSVLogger();
	CSVLogger(string* header, unsigned columnNumber);
	CSVLogger(string& filename, string* header, unsigned columnNumber);
	CSVLogger(char seperator, string& filename, string* header, unsigned columnNumber);
	~CSVLogger();

	// Data writing methods
	void addData(const char data);
	void addData(const string& data);
	void addData(const int data);
	void addData(const double data);
	void addData(const bool data);
	void addData(const void* data);
	bool addLine(const string& line);
	void addSeperator();
	void addNewLine();

	// getter methods
	const NumberFormat getNumberFormat() const {return numberFormat;}
	const string getNumberFormatAsString() const;
	const string& getBoolTrueString() const {return boolTrueString;}
	const string& getBoolFalseString() const {return boolFalseString;}
	const string& getLogFileName() const {return logFileName;}
	const unsigned getPrecision() const {return precision;}
	const bool getHeaderWritten() const {return headerWritten;}
	bool isLogOpen() {return logFile.is_open();}

	// setter methods
	void setNumberFormat(const NumberFormat numberFormat) {this->numberFormat = numberFormat;}
	void setBoolTrueString(const string& boolTrueString) {this->boolTrueString = boolTrueString;}
	void setBoolFalseString(const string& boolFalseString) {this->boolFalseString = boolFalseString;}
	void setBoolStrings(const string& boolTrueString, const string& boolFalseString);
	void setPrecision(const unsigned precision) {this->precision = precision; this->logFile.precision(precision);}

	private:
	/** The number of columns in the table. */
	unsigned columnNumber;

	/** The index of the corrent column from the interval [0,columnNumber - 1]. */
	unsigned currentColumn;

	/** Field seperator of the table. */
	char seperator;

	/** The log file. */
	ofstream logFile;

	/** The name of the log file. */
	string logFileName;

	/** The number format for double values. */
	NumberFormat numberFormat;

	/** Precision for double values. */
	unsigned precision;

	/** String to write out for true bool values. */
	string boolTrueString;

	/** String to write out for true bool values. */
	string boolFalseString;

	/** The table header. */
	string* header;

	/** Determines whether the header was already written to the log. */
	bool headerWritten;

	string changeNumberFormat(double number);
	void writeHeader();
};

#endif
