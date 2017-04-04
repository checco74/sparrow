#ifndef _PREDICTION_FILE_H_
#define _PREDICTION_FILE_H_

#include <string>
#include <list>
#include <iostream>

typedef enum SSEEntryType
{
  SSE_HELIX,
  SSE_STRAND,
  SSE_COIL
};

/**
 * A prediction file line represents a single residue in a predicted protein and contains an amino acid one letter code for this residue as
 * well as its structure type and all scores generated during the prediction.
 * @brief Class representation of a line representing a residue in a prediction output file.
 */
class PredictionFileEntry
{
  public:
    PredictionFileEntry(std::string& line);
    PredictionFileEntry(char aa, SSEEntryType type, double* scores, bool copy = true);
    ~PredictionFileEntry();

    char getAA() const;
    SSEEntryType getType() const;
    char getTypeChar() const;
    double getScore(unsigned i) const;
    double* getScores() const;

    static SSEEntryType toType(char c);
    static SSEEntryType toType(unsigned u);

    /** @brief The number scores belonging to a prediction file entry. */
    static const unsigned NUMBER_OF_SCORES = 3;

  protected:

    /** @brief The amino acid one letter code for the residue represented by this prediction file line. */
    char aa;
    /** @brief The secondary structure type of the residue represented by this prediction file line. */
    SSEEntryType type;
    /** @brief The scores for the residue represented by this prediction file line. */
    double* scores;
    /** @brief Character representation of a residue predicted as helix. */
    static const char HELIX_CHAR = 'H';
    /** @brief The index of the prediction entry corresponding to the helix score. */
    static const unsigned HELIX_INDEX = 0;
    /** @brief Character representation of a residue predicted as strand. */
    static const char STRAND_CHAR = 'E';
    /** @brief The index of the prediction entry corresponding to the strand score. */
    static const unsigned STRAND_INDEX = 1;
    /** @brief Character representation of a residue predicted as coil. */
    static const char COIL_CHAR = 'C';
    /** @brief The index of the prediction entry corresponding to the coil score. */
    static const unsigned COIL_INDEX = 2;
};

std::ostream& operator<<(std::ostream& s, const PredictionFileEntry& entry);

/**
 * A prediction output file consists of two lines of header containing the name an the length of the predicted protein as well as one entry per residue
 * describing its amino acid, its secondary structure type and the scores it obtained from the scoring functions. The entries representing residues are
 * saved in <code>PredictionFileEntry</code> objects.
 * @brief A class representation of a prediction output file.
 * @see PredictionFileEntry
 */
class PredictionFile
{
  public:
    PredictionFile(std::string& filename);
    ~PredictionFile();

    unsigned getProteinLength() const;
    std::string getProteinName() const;
    std::string getSequence() const;
    std::string getStructure() const;
    PredictionFileEntry* getEntry(unsigned i) const;
    std::list<PredictionFileEntry*>* getEntries() const;

    void setProteinName(std::string& name);
    void setProteinLength(unsigned length);
    void addEntry(char aa, SSEEntryType type, double* scores, bool copy = true);

    void write(std::string& filename) const;
    bool read() const;

  protected:
    /** @brief The name of the file represented by this <code>PredictionFile</code> object. */
    std::string proteinName;
    /** @brief The length of the represented protein. */
    unsigned proteinLength;
    /** @brief A list containing all entries representing the residues of the represented protein. */
    std::list<PredictionFileEntry*> entries;
};

std::ostream& operator<<(std::ostream& s, const PredictionFile& file);

#endif
