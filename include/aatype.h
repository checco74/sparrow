#ifndef _AMINO_ACID
#define _AMINO_ACID

#include <math.h>
#include "define.h"

using namespace std;

class aminoacid
{
	public:
		int numerical;
	char oneletter;
	static const string onelettercodes;
	static const uint maxValue = 111;
	static double* logValue;
	static void init();
	static void kill();

	aminoacid();
	~aminoacid();

	static int astype(char letter);
	static int *astype(char *letters, int length);
	static char ascode(int number);
	static char *ascode(int *numbers, int length);
};

#endif
