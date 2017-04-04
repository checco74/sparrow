#include "preprocessor.h"

#ifdef PREPROCESSOR_DEBUG_MODE
const string LOG_NAME = "preprocessor.log";
const string END_OF_BLOCK = "===========================================================================================";
#endif

/**
 * Standard constructor
 *
 * @param filename The name of the database file where to read the strings from
 * @param removeMax The maximum size of a string containing equal characters, that should be replaced/removed
 * @param neighbourMin If the neighbouring strings of the currently processed string are larger than this number, the current string's characters should be replaced by the neighbouring strings character
 */
Preprocessor::Preprocessor(const char* filename, unsigned removeMax, unsigned neighbourMin)
{
	this->removeMax = removeMax;
	this->neighbourMin = neighbourMin;
	this->dataBase.open(filename, fstream::in);
	if (!this->dataBase.is_open())
	{
	cout << "Preprocessor could not open " << filename << " for reading!" << endl;
	exit(1);
	}
#ifdef PREPROCESSOR_DEBUG_MODE
	this->log.open(LOG_NAME.c_str(), fstream::out | fstream::trunc);
	this->log << "Preprocessor log" << endl;
	this->log << "First line of a block: The unprocessed string" << endl;
	this->log << "Second line of a block: Either the processed string if not using transcoded mode or transcoded string otherwise" << endl;
	this->log << "Third line of a block: The processed string if using transcoded mode" << endl;
	this->log << END_OF_BLOCK << endl;
	this->log << "BEGIN OF PREPROCESSOR LOG" << endl << endl;
#endif
	this->transcoderMap = NULL;
	this->defaultTranscoderChar = DEFAULT_CHAR;
}

/**
 * Standard destructor
 */
Preprocessor::~Preprocessor()
{
	this->dataBase.close();
	if (this->transcoderMap != NULL)
	delete this->transcoderMap;
#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << "END OF PREPROCESSOR LOG" << endl;
	this->log.close();
#endif
}

/**
 * Gets the next string from the database and preprocesses it.
 *
 * @param transcoded Determines whether to transcode the string before preprocessing
 * @return The processed string or an empty string if no more strings can be processed
 */
string Preprocessor::getNextPreprocessedString(bool transcoded)
{
	string rawString = this->getNextLine();
	if (rawString != "" && rawString != "\n")
	{
	#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << rawString << endl;
	#endif
	if (transcoded)
	{
	rawString = this->transcode(rawString);
	#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << rawString << endl;
	#endif
	}
	#ifdef PREPROCESSOR_DEBUG_MODE
	string preprocessed = this->processString(rawString);
	this->log << preprocessed << endl << END_OF_BLOCK << endl;
	return preprocessed;
	#else
	return this->processString(rawString);
	#endif
	}
	else
	return "";
}

/**
 * Cuts a string into two substrings of which the returned one contains equal characters only.
 *
 * @param s The string to cut, contains the remaining string after the operation
 * @return A substring which is taken from the front of s and consists of equal characters only
 */
string Preprocessor::cutString(string& s)
{
	string newString = "";
	unsigned index = 0;
	while (s.size() != 0 && s[index] == s[0]) newString += s[index++];
	s = s.substr(newString.size(), s.size() - newString.size());
	return newString;
}

/**
 * Gets the next '|' delimited line from the database.
 *
 * @return The next line from the database
 */
string Preprocessor::getNextLine()
{
	assert(this->dataBase.is_open());
	string line = "";
	char current;
	bool done = false;
	while (!this->dataBase.eof() && !done)
	{
	this->dataBase >> current;
	if (current != '\n' && current != '|')
	line += current;
	else if (current == '|')
	done = true;
	}
	return line;
}

/**
 * Determines whether the preprocessing of the current database is done.
 *
 * @return true if the entire database has been processed
 */
bool Preprocessor::preprocessingDone()
{
	return this->dataBase.eof();
}

/**
 * Creates a string with the length number consisting of the character c only.
 *
 * @param c The character of which the string should consist
 * @param number The length of the created string
 * @return A string consisting of number copies of c
 */
string Preprocessor::createString(char c, unsigned number)
{
	string s = "";
	for (unsigned i = 0; i < number; i++)
	s += c;
	return s;
}

/**
 * Transcodes a string by use of the transcode()-method on each of the string's characters.
 *
 * @param s The string to transcode
 * @return The transcoded string
 */
string Preprocessor::transcode(string& s)
{
	string newString = "";
	for (unsigned i = 0; i < s.size(); i++) newString += transcode(s[i]);
	return newString;
}

/**
 * Transodes a char by using a included code table or a user specified code table.
 * The standard transcoding uses the DSSP code for SSEs and puts them into 3 SSEs: H (helix), e (strand) and . (coil)
 *
 * @param c The character to transcode
 * @return The transcoded character
 */
char Preprocessor::transcode(char c)
{
	if (this->transcoderMap == NULL)
	{
	switch (c)
	{
	case 'H': // alpha-helix
	case 'G': // 3-helix(3/10 helix)
	case 'I': // 5-helix(pi helix)
	return 'H';
	case 'E': //extended strand, participates in beta ladder
	case 'B': //residue in isolated beta-bridge
	case 'e':
	return 'e';
	case ' ': // coil
	case 'T': // hydrogen bonded turn
	case 'S': // bend
	case '.':
	case ':':
	return '.';
	default:
	return '.';
	}
	}
	else
	{
	map<char, char>::iterator found = this->transcoderMap->find(c);
	if (found != this->transcoderMap->end())
	return found->second;
	else
	return this->defaultTranscoderChar;
	}
}

/**
 * Prerprocesses a given string by replacing all characters that do not fit into the string.
 * For example with removeMax = 2 and neighbouringMin = 4 the following strings are produced:
 *	HH.HHH..SSTTT => HHHHHHHTTTTT
 *	HHHS.SHHTT => HHHSSSHHTT
 *
 * @param s The string to process
 * @return The processed string
 */
string Preprocessor::processString(string& s)
{
	string rawString = s;
	string preprocessedString = "";

	string previous = "";
	string current = this->cutString(rawString);
	string next = this->cutString(rawString);

	while (current.size() > 0)
	{
	bool merge = false;
	if (current.size() <= this->removeMax && current.substr(0, 1) != previous.substr(0, 1))
	{
	if (previous.substr(0, 1) == next.substr(0, 1))
	{
	if (previous.size() + current.size() + next.size() > this->removeMax)
	{
	current = this->createString(previous[0], current.size());
	merge = true;
	}
	}
	}
	preprocessedString += current;
	if (!merge)
	previous = current;
	else
	previous += current;
	current = next;
	next = this->cutString(rawString);
	}

	rawString = preprocessedString;
	preprocessedString = "";
	previous = "";
	current = this->cutString(rawString);
	next = this->cutString(rawString);

	while (current.size() > 0)
	{
	bool merge = false;
	if (current.size() <= this->removeMax && current.substr(0, 1) != previous.substr(0, 1))
	{
	if (previous.size() >= next.size())
	{
	if (previous.size() >= this->neighbourMin || previous.size() + current.size() > this->removeMax)
	{
	current = this->createString(previous[0], current.size());
	merge = true;
	}
	}
	else
	{
	if (next.size() >= this->neighbourMin || next.size() + current.size() > this->removeMax)
	{
	current = this->createString(next[0], current.size());
	merge = true;
	}
	}
	}
	preprocessedString += current;
	if (!merge)
	previous = current;
	else
	previous += current;
	current = next;
	next = this->cutString(rawString);
	}

	return preprocessedString;
}

string Preprocessor::processString(string& s, bool transcoding)
{
	string newString = s;
#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << newString << endl;
#endif
	if (transcoding)
#ifdef PREPROCESSOR_DEBUG_MODE
	{
	newString = this->transcode(newString);
	this->log << newString << endl;
	}
	newString = this->processString(newString);
	this->log << newString << endl << END_OF_BLOCK << endl;
	return newString;
#else
	newString = this->transcode(newString);
	return this->processString(newString);
#endif
}

/**
 * Creates a new transcoder map if necessary and sets the default transcoder character.
 *
 * @param defaultChar The new default transcoder character
 */
void Preprocessor::createTranscoderMap(char defaultChar)
{
	if (this->transcoderMap == NULL)
	{
	this->transcoderMap = new map<char, char>;
#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << "Created transcoder map..." << endl;
#endif
	}
	this->defaultTranscoderChar = defaultChar;
#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << "Added transcoder default char: " << defaultChar << endl;
#endif
}

/**
 * Adds a new transcoder mapping to the transcoder map. If no transcoder map exists a new one is created.
 *
 * @param from The original, not transcoded character
 * @param to The character that from should be transcoded to
 */
void Preprocessor::addToTranscoderMap(char from, char to)
{
	if (this->transcoderMap == NULL)
	{
	this->transcoderMap = new map<char, char>;
#ifdef PREPROCESSOR_DEBUG_MODE
	this->log << "Created transcoder map..." << endl;
#endif
	}
	pair<map<char, char>::iterator, bool> inserted = this->transcoderMap->insert(make_pair(from, to));
#ifdef PREPROCESSOR_DEBUG_MODE
	if (inserted.second)
	this->log << "Added new transcoder mapping: " << from << " -> " << to << endl;
	else
	this->log << "Could not add transcoder mapping: " << from << " -> " << to << endl;
#endif
}
