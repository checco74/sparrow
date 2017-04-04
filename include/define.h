#ifndef _STNDEFIN_
#define _STNDEFIN_

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>

#include <sys/wait.h>
#include <sys/types.h>

#define NUM_OF_MOTIVES 3
#define NUM_OF_AMINO_ACIDS 20
#define NO_TRANSITIONS 1
#define RND drand48()
#define TINY 1.0e-20
#define MAXLENGTH 21
#define DEBUG_MODE

const char SEPARATOR = '/';
const std::string DEFAULT_PROFILE = "EEEEEEEEEEEEEEEEEEEE";

typedef enum detailLevel // the higher the level the less detailed the secondary structure distinction
{
	FULL,	// Level 0 - 21 classes - full definition: all StructureTypes are used
	ANY_COIL,	// Level 1 - 20 classes - coils are not separated between short and long
	NO_BENDS,	// Level 2 - 19 classes - coils (of any length) AND bends are grouped together
	NO_TURNS,	// Level 3 - 18 classes - coils (of any length) AND bends AND turns are grouped together
	NO_REGIONS,	// Level 4 - 13 classes - motive patterns are not separated into regions
	NO_REGIONS_ANY_COIL,	// Level 5 - 12 classes - like the previous one MINUS coil length distinction
	NO_REGIONS_NO_BENDS,	// Level 6 - 11 classes - like the previous one MINUS bends distinction
	NO_REGIONS_NO_TURNS,	// Level 7 - 10 classes - like the previous one MINUS turns distinction
	DSSP_PLUS_SHORTS,	// Level 8 - 9 classes - DSSP standard definitios PLUS short coils
	DSSP_STANDARD,	// Level 9 - 8 classes - length and region aren't taken into account
	DSSP_NO_BENDS,	// Level 10 - 7 classes - DSSP without bends
	DSSP_NO_TURNS,	// Level 11 - 6 classes - DSSP without bends or turns
	DSSP_NO_EXOTIC,	// Level 12 - 5 classes - DSSP without rare stuff
	Q3_PLUS_TURNS,	// Level 13- 4 classes
	Q3_STANDARD	// Level 14 - 3 classes
};

/**
 * @enum StructureType Classification of the structural type of a residue.
 */
typedef enum StructureType
{
  SHORT_HELIX,             // short helix residue
  START_HELIX,             // helix residue at the beginning of a helix which is not short
  CENTER_HELIX,            // helix residue in the center of a helix which is not short
  END_HELIX,               // helix residue at the end of a helix which is not short
  SHORT_STRAND,            // short strand residue
  START_STRAND,            // strand residue at the beginning of a strand which is not short
  CENTER_STRAND,           // strand residue in the center of a strand which is not short
  END_STRAND,              // strand residue at the end of a strand which is not short
  ISOLATED_STRAND,         // strand residue being part of an isolated strand
  RANDOM_COIL,             // random coil residue
  SHORT_T_HELIX,           // short 3-helix residue
  START_T_HELIX,           // helix residue at the beginning of a 3-helix which is not short
  CENTER_T_HELIX,          // helix residue in the center of a 3-helix which is not short
  END_T_HELIX,             // helix residue at the end of a 3-helix which is not short
  SHORT_P_HELIX,           // short 5-helix residue
  START_P_HELIX,           // helix residue at the beginning of a 5-helix which is not short
  CENTER_P_HELIX,          // helix residue in the center of a 5-helix which is not short
  END_P_HELIX,             // helix residue at the end of a 5-helix which is not short
  TURN,                    // residue on an hydrogen-bonded turn
  BEND,                    // residue on a bend
  SHORT_COIL,              // short coil residue
  NUM_OF_STRUCTURE_TYPES   // This entry is not a structural type itself, but represents the number of different
                           // structural types existing (new types should be inserted before this entry)
};

typedef enum DSSPStructureType
{
	ALPHA_HELIX,	// long helix residue
	SHORT_ALPHA_HELIX,	// short helix residue
	BETA_STRAND,	// long strand residue
	SHORT_BETA_STRAND,	// short strand residue
	THIN_HELIX,	// long 3-helix residue
	SHORT_THIN_HELIX,	// short 3-helix residue
	FAT_HELIX,	// long 5-helix residue
	SHORT_FAT_HELIX,	// short 5-helix residue
	BETA_BRIDGE,	// residue in an isolated strand
	GENERIC_COIL,	// residue in a long coil
	TURN_LOOP,	// residue on an hydrogen-bonded turn
	BEND_LOOP,	// residue on a bend
	SHORT_GENERIC_COIL	// residue in a short coil
};

typedef enum DSSPStandardType
{
	DSSP_HELIX,	// helix residue
	DSSP_STRAND,	// strand residue
	DSSP_T_HELIX,	// 3-helix residue
	DSSP_P_HELIX,	// 5-helix residue
	DSSP_BRIDGE,	// residue in an isolated strand
	DSSP_COIL,	// residue in a long coil
	DSSP_TURN,	// residue on an hydrogen-bonded turn
	DSSP_BEND,	// residue on a bend
	DSSP_SHORT	// residue in a short coil
};

typedef enum Q3ExtraType
{
	_HELIX,	// long helix residue
	_STRAND,	// long strand residue
	_COIL,	// residue in a long coil
	_TURN,	// residue on an hydrogen-bonded turn
	_BEND	// residue on a bend
};

typedef enum Category
{
	HELIX,	// Residue is a helical residue
	STRAND,	// Residue is a beta strand residue
	COIL,	// Residue is a random coil
	T_HELIX,	// Residue is a 3-helix
	P_HELIX,	// Residue is a 5-helix
	TURN_COIL,	// Residue is a turn coil
	BEND_COIL,	// Residue is a bend coil
	ISOLATED_BRIDGE	// Residue is an isolated beta bridge
};

typedef enum PositionType
{
	C_TERM, // Residue belongs to the C-Terminus of a protein
	N_TERM, // Residue belongs to the N-Terminus of a protein
	CENTER	// Residue belongs to the center of a protein
};

typedef enum Position
{
	CHAIN_C_TERM, // Residue belongs to the C-terminal region of a protein
	CHAIN_N_TERM, // Residue belongs to the N-terminal region of a protein
	CHAIN_CENTER	// Residue belongs to the center of a protein
};

typedef enum inputMode // identifies the type of input
{
	FASTA_MODE, PSIBLAST_MODE
};

typedef enum readMode // the configuration file reading modes
{
	BASIC, STD_MODE, AUXN_MODE, AUXA_MODE
};

typedef enum CATHphase
{
	LEARN, ASSESS, EMPLOY
};

/**
 * @brief Struct to hold the data of a whole protein.
 */
typedef struct Protein
{
  /** @brief The database ID of the protein. */
  unsigned id;
  /** @brief The sequence data for the protein. */
  std::string sequence;
  /** @brief The structure of the protein as string of <code>StructureType</code>s. */
  std::string structure;
  /** @brief The position types of the protein's residues */
  std::string positions;
  /** @brief The profile for this protein. */
  std::string profile;
};

/**
 * All exceptions that might occur are derived from this class.
 * @brief A generic exception.
 */
class GenericException : public std::exception
{
  public:
    GenericException() throw();
    GenericException(const char *msg) throw();

    virtual ~GenericException() throw();

    const char* what() const throw();
    virtual std::string getType() const throw();
  protected:
    /** @brief The exception's error message. */
    char   *msg;
};

#endif

