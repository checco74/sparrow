#include <string>
#include <fstream>
#include <iostream>
#include "assert.h"
#include <map>
//#include <pair>

using namespace std;

#define PREPROCESSOR_DEBUG_MODE

class Preprocessor
{
	public:
	Preprocessor(const char* filename, unsigned removeMax = 1, unsigned neighbourMin = 4);
	~Preprocessor();
	string getNextPreprocessedString(bool transcoded = true);
	string processString(string& s, bool transcoding);
	string transcode(string& s);
	char transcode(char c);
	bool preprocessingDone();
	void createTranscoderMap(char defaultChar);
	void addToTranscoderMap(char from, char to);

	protected:
#ifdef PREPROCESSOR_DEBUG_MODE
	/** The debugging log file. */
	ofstream log;
#endif
	/** The input database file. */
	ifstream dataBase;
	/** The maximum number of neighbouring equal characters that should be replaced. */
	unsigned removeMax;
	/** The minimum number of neighbouring equal characters necessary to use them as replacement. */
	unsigned neighbourMin;
	/** A map that enables a user defined transcoding */
	map<char, char>* transcoderMap;
	/** The character used to encode characters for which no mapping is known from the transcoder map. */
	char defaultTranscoderChar;
	/** The standard default transcoder char. */
	static const char DEFAULT_CHAR = '.';

	string cutString(string& s);
	string getNextLine();
	string createString(char c, unsigned number);
	string processString(string& s);
};
