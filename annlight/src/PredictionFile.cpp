#include "PredictionFile.h"
#include "IOException.h"
#include <sstream>
#include <fstream>

/**
 * @brief Standard constructor of the <code>PredictionFileEntry</code> class.
 *
 * @param line A string representation of the current prediction file line as read from a prediction file
 */
PredictionFileEntry::PredictionFileEntry(std::string& line)
{
  if(line.size() == 0)
    throw IOException("Passed empty string to PredictionFileEntry::PredictionFileEntry()");
  this->aa = line[0];
  if (line.size() == 0)
    throw IOException("Missing secondary structure in prediction file line!");
  this->type = PredictionFileEntry::toType(line[2]);

  this->scores = new double[PredictionFileEntry::NUMBER_OF_SCORES];
  std::string::size_type previousPos = 4;
  for (unsigned i = 0; i < PredictionFileEntry::NUMBER_OF_SCORES; i++)
  {
    if (i < PredictionFileEntry::NUMBER_OF_SCORES - 1)
    {
      std::string::size_type currentPos = line.find(' ', previousPos);
      if (currentPos == std::string::npos)
        throw IOException("Missing score in prediction file line!");
      this->scores[i] = atof((line.substr(previousPos, currentPos - previousPos)).c_str());
      previousPos = currentPos + 1;
    }
    else
    {
      this->scores[i] = atof((line.substr(previousPos, line.size() - previousPos)).c_str());
    }
  }
}

/**
 * Creates a new <code>PredictionFileEntry</code> according to the given data.
 * @brief Constructor of the <code>PredictionFileEntry</code>.
 *
 * @param aa The amino acid one letter code for the prediction file line
 * @param type The secondary structure type for the prediction file line
 * @param scores The scores for the prediction file line
 * @param copy Determines whether to save a copy of the <code>scores</code> parameter in the <code>PredictionFileEntry</code> or to save the pointer directly
 */
PredictionFileEntry::PredictionFileEntry(char aa, SSEEntryType type, double* scores, bool copy)
{
  this->aa = aa;
  this->type = type;

  if (copy)
  {
    this->scores = new double[PredictionFileEntry::NUMBER_OF_SCORES];
    for (unsigned i = 0; i < PredictionFileEntry::NUMBER_OF_SCORES; i++)
      this->scores[i] = scores[i];
  }
  else
    this->scores = scores;
}

/**
 * @brief Standard destructor for the <code>PredictionFileEntry</code> class.
 */
PredictionFileEntry::~PredictionFileEntry()
{
  delete [] this->scores;
}

/**
 * @brief Returns the one letter code for the represented residue.
 *
 * @return The one letter code for the represented residue
 */
char PredictionFileEntry::getAA() const
{
  return this->aa;
}

/**
 * @brief Returns the secondary structure type of the represented residue.
 *
 * @return The secondary structure type of the represented residue
 */
SSEEntryType PredictionFileEntry::getType() const
{
  return this->type;
}

/**
 * @brief Returns the character representation of the secondary structure type of the represented residue.
 *
 * @return The character representation of the structure type of the represented residue
 */
char PredictionFileEntry::getTypeChar() const
{
  switch (this->type)
  {
    case SSE_HELIX:
      return PredictionFileEntry::HELIX_CHAR;
    case SSE_STRAND:
      return PredictionFileEntry::STRAND_CHAR;
    case SSE_COIL:
      return PredictionFileEntry::COIL_CHAR;
    default:
      return PredictionFileEntry::COIL_CHAR;
  }
}

/**
 * @brief Returns one of the scores of the represented residue
 *
 * @param i The index of the score of interest (must be less than <code>PredictionFileEntry::NUMBER_OF_SCORES</code>)
 */
double PredictionFileEntry::getScore(unsigned i) const
{
  if (i >= PredictionFileEntry::NUMBER_OF_SCORES)
  {
    std::ostringstream msg;
    msg << "Requested score number " << i << " does not exist!";
    throw IOException(msg.str().c_str());
  }
  return this->scores[i];
}

/**
 * @brief Returns all scores of the prepresented residue
 *
 * @return A copy of the nested <code>scores</code> class member
 */
double* PredictionFileEntry::getScores() const
{
  double* retScores = new double[PredictionFileEntry::NUMBER_OF_SCORES];
  for (unsigned i = 0; i < PredictionFileEntry::NUMBER_OF_SCORES; i++)
    retScores[i] = this->scores[i];
  return retScores;
}

/**
 * @brief Converts a character representation of a <code>SSEEntryType</code> into a <code>SSEEntryType</code>
 *
 * @param c The character representation of the <code>SSEEntryType</code>
 * @return The <code>SSEEntryType</code> corresponding to <code>c</code>
 */
SSEEntryType PredictionFileEntry::toType(char c)
{
  switch (c)
  {
    case PredictionFileEntry::HELIX_CHAR:
      return SSE_HELIX;
    case PredictionFileEntry::STRAND_CHAR:
      return SSE_STRAND;
    case PredictionFileEntry::COIL_CHAR:
      return SSE_COIL;
    default:
      std::ostringstream msg;
      msg << "Unknown secondary structure letter type: " << c;
      throw IOException(msg.str().c_str());
  }
}

/**
 * @brief Converts an score index into a <code>SSEEntryType</code>
 *
 * @param u The index of a score within the prediction
 * @return The <code>SSEEntryType</code> corresponding to <code>c</code>
 */
SSEEntryType PredictionFileEntry::toType(unsigned u)
{
  switch (u)
  {
    case PredictionFileEntry::HELIX_INDEX:
      return SSE_HELIX;
    case PredictionFileEntry::STRAND_INDEX:
      return SSE_STRAND;
    case PredictionFileEntry::COIL_INDEX:
      return SSE_COIL;
    default:
      std::ostringstream msg;
      msg << "Unknown secondary structure index type: " << u;
      throw IOException(msg.str().c_str());
  }
}

std::ostream& operator<<(std::ostream& s, const PredictionFileEntry& entry)
{
  s << entry.getAA() << " " << entry.getTypeChar();
  for (unsigned i = 0; i < PredictionFileEntry::NUMBER_OF_SCORES; i++)
    s << " " << entry.getScore(i);
  return s;
}

/**
 * If the passed file name does not correspond to any existing file no data is read.
 * @brief Standard constructor for the <code>PredictionFile</code> class.
 *
 * @param filename The name of the prediction file from which to read the data
 */
PredictionFile::PredictionFile(std::string& filename)
{
  std::ifstream inFile(filename.c_str(), std::fstream::in);

  if (inFile.is_open())
  {
    std::string line;
    getline(inFile, line);
    this->proteinName = line;
    getline(inFile, line);
    this->proteinLength = atoi(line.c_str());

    unsigned read = 0;
    while (getline(inFile, line) && read < this->proteinLength)
    {
      if (line != "")
      {
        PredictionFileEntry* entry = new PredictionFileEntry(line);
        this->entries.push_back(entry);
        read++;
      }
    }

    if (read != this->proteinLength)
    {
      std::ostringstream msg;
      msg << "Missing entries in file \"" << filename << "\": read " << read << " of " << this->proteinLength;
      throw IOException(msg.str().c_str());
    }

    inFile.close();
  }
  else
  {
    this->proteinName = "";
    this->proteinLength = 0;
  }
}

/**
 * @brief The standard destructor for the <code>PredictionFile</code> class.
 */
PredictionFile::~PredictionFile()
{
  for (std::list<PredictionFileEntry*>::iterator it = this->entries.begin(); it != this->entries.end(); it++) delete (*it);
}

/**
 * @brief Returns the length of the protein read from the represented prediction output file.
 *
 * @return The value of the <code>proteinLength</code> member variable
 */
unsigned PredictionFile::getProteinLength() const
{
  return this->proteinLength;
}

/**
 * @brief Returns the name of the protein read from the represented prediction output file.
 *
 * @return The value of the <code>proteinName</code> member variable
 */
std::string PredictionFile::getProteinName() const
{
  return this->proteinName;
}

/**
 * @brief The sequence of the protein read from the represented prediction output file.
 *
 * @return A string representation of the sequence of the read protein
 */
std::string PredictionFile::getSequence() const
{
  std::ostringstream seq;
  for(std::list<PredictionFileEntry*>::const_iterator it = this->entries.begin(); it != this->entries.end(); it++)
    seq << (*it)->getAA();
  return seq.str();
}

/**
 * @brief The structure of the protein read from the represented prediction output file.
 *
 * @return A string representation of the structure of the read protein
 */
std::string PredictionFile::getStructure() const
{
  std::ostringstream str;
  for(std::list<PredictionFileEntry*>::const_iterator it = this->entries.begin(); it != this->entries.end(); it++)
    str << (*it)->getTypeChar();
  return str.str();
}

/**
 * @brief Returns a prediction file entry corresponding to a specific residue.
 *
 * @param i The ID of the residue for which the prediction file entry should be returned
 * @return The requested prediction file entry
 */
PredictionFileEntry* PredictionFile::getEntry(unsigned i) const
{
  unsigned k = 0;
  std::list<PredictionFileEntry*>::const_iterator it = this->entries.begin();
  while (k < i) it++;
  return *it;
}

/**
 * @brief Returns a list of all entries in the read prediction file.
 *
 * @return A copy of the <code>entries</code> member variable
 */
std::list<PredictionFileEntry*>* PredictionFile::getEntries() const
{
  std::list<PredictionFileEntry*>* copy = new std::list<PredictionFileEntry*>;
  for(std::list<PredictionFileEntry*>::const_iterator it = this->entries.begin(); it != this->entries.end(); it++)
    copy->push_back(*it);
  return copy;
}

/**
 * @brief Sets the name of the represented protein.
 *
 * @param name The name of the represented protein
 */
void PredictionFile::setProteinName(std::string& name)
{
  this->proteinName = name;
}

/**
 * @brief Sets the length of the represented protein.
 *
 * @param length The length of the represented protein
 */
void PredictionFile::setProteinLength(unsigned length)
{
  this->proteinLength = length;
}

/**
 * @brief Adds an entry for a new residue to the represented prediction file.
 *
 * @param aa The amino acid one letter code of the new residue
 * @param type The secondary structure type of the new residue
 * @param scores The scores of the new residue
 * @param copy Determines whether to save a copy of <code>scores</code> or save it directly
 */
void PredictionFile::addEntry(char aa, SSEEntryType type, double* scores, bool copy)
{
  PredictionFileEntry* entry = new PredictionFileEntry(aa, type, scores, copy);
  this->entries.push_back(entry);
}

/**
 * @brief Writes the data of the <code>PredictionFile</code> into a file.
 *
 * @param filename The name of the file to write
 */
void PredictionFile::write(std::string& filename) const
{
  std::ofstream outFile(filename.c_str(), std::fstream::out | std::fstream::trunc);

  if (!outFile.is_open())
  {
    std::ostringstream msg;
    msg << "Could not open file \"" << filename << "\" for writing!";
    throw IOException(msg.str().c_str());
  }

  outFile << *this;

  outFile.close();
}

/**
 * @brief Determines whether this <code>PredictionFile</code> contains any data (so if any data has been read)
 *
 * @return <code>true</code> if the <code>PredictionFile</code> contains any data, <code>false</code> otherwise
 */
bool PredictionFile::read() const
{
  return (this->proteinName != "" && this->proteinLength != 0 && this->entries.size() != 0);
}

std::ostream& operator<<(std::ostream& s, const PredictionFile& file)
{
  s << file.getProteinName() << std::endl;
  s << file.getProteinLength() << std::endl;
  std::list<PredictionFileEntry*>* entries = file.getEntries();
  unsigned i = 1;
  for (std::list<PredictionFileEntry*>::const_iterator it = entries->begin(); it != entries->end(); it++)
  {
    s << i << " " << *(*it) << std::endl;
  }
  delete entries;
  return s;
}
