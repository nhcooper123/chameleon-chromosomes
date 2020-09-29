#ifndef ___CHR_COUNT_FORMAT
#define ___CHR_COUNT_FORMAT

#include "sequenceContainer.h"
#include "finiteIntAlphabet.h"

class chrCountFormat {
public:
	static sequenceContainer read(istream &infile, finiteIntAlphabet* alph);
	//readUnAligned: the input sequences do not need to be aligned (not all sequences are the same length).
	static sequenceContainer readUnAligned(istream &infile, finiteIntAlphabet* alph);
	//returns the minimum and maximum counts observed in the file
	static pair<int,int> getMinMaxCountInFile(const string &infilName);
	static void write(ostream &out, const sequenceContainer& sd);
};

#endif

/* This format is similar to FASTA with 2 exceptions:
1. only one position is allowed
2. Each sequence may be represented by more than one count.
In this case the possible counts are seperated by '_': 
	3_5 : 3 or 5 BOTH may be present. 
	3=0.5_4=0.5 : Counts are sampled from a popuation. There is 50% chance that the count for that taxa is 3 and 50% that it is 5.
EXAMPLE OF chrCount FORMAT:
>Langur
18
>Baboon
18
>Human
23=0.4_12=0.6
>Rat
23_24
>Cow
25
>Horse
25
*/
