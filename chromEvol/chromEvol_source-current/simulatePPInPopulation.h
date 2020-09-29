#ifndef ___SIMULATE_PP__
#define ___SIMULATE_PP__

#include "definitions.h"

#include <vector>
using namespace std;

/******************************************************************
This class simulates the number of PP events within a population of a specified size that started from n diploid individuals. 
The inputs are: (1) speciation rate (2) extinction rate (3) polyploid speciation frequency (4) final population size
The out put: the number of extant taxa in each ploidy level
*******************************************************************/

class simulatePPInPopulation  {
public:
	simulatePPInPopulation(MDOUBLE speciationRate, MDOUBLE extinctionRate, MDOUBLE ppSpeciation, int finalPopSize);
	virtual ~simulatePPInPopulation();
	void run(int startDiploidNum, const string& outfile);
	void run(int startDiploidNum, string const& statsFile, const string& outfile);

private:
	void init(int startDiploidNum);
	int drawBin();
	void printResults(const string& outfile);
	void speciateTaxa(int bin);
    void extinctTaxa(int bin);
	void readDiversificationFile(const string& diversiifcationFile);
	//in case diversification distribution is given - decide which rates to use. Otherwise - return _specRate + _extRate 
	pair<MDOUBLE, MDOUBLE> getDiversificationRates(); 

	
private:
	MDOUBLE _specRate; //speciation rate
	MDOUBLE _extRate; //extinction rate
	Vdouble _relativeExtRate;
	Vdouble _speciesRichnessDist;
	MDOUBLE _ppSpeciation;
	int _curTaxaNum; 
	int _finalTaxaNum;
public: 
	Vint _ppDist;
	MDOUBLE _avg;
};

#endif
