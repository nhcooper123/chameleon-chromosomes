#include "simulatePPInPopulation.h"
#include "talRandom.h"
#include "errorMsg.h"
#include "someUtil.h"
#include "ConversionUtils.h"

#define MAX_PP 30



simulatePPInPopulation::simulatePPInPopulation(MDOUBLE speciationRate, MDOUBLE extinctionRate, MDOUBLE ppSpeciation, int finalPopSize)
: _specRate(speciationRate), _extRate(extinctionRate), _ppSpeciation(ppSpeciation), _finalTaxaNum(finalPopSize)	
{
}

simulatePPInPopulation::~simulatePPInPopulation()
{
}

void simulatePPInPopulation::init(int startDiploidNum)
{
	_curTaxaNum = startDiploidNum;
	_ppDist.resize(MAX_PP, 0);
	_ppDist[0] = _curTaxaNum;
}

void simulatePPInPopulation::run(int startDiploidNum, string const& statsFile, const string& outfile)
{
	readDiversificationFile(statsFile);
	run(startDiploidNum, outfile);
}
void simulatePPInPopulation::run(int startDiploidNum, const string& outfile)
{
	init(startDiploidNum);
	for (int loop = 0 ; loop < 100000 ; ++loop){
		while (_curTaxaNum < _finalTaxaNum) {
			int bin = drawBin();
			//decide if the taxa will speciate or go extinct
			pair<MDOUBLE, MDOUBLE> rates = getDiversificationRates(); //in case diversification distribution is given - decide which rates to use 
			MDOUBLE rand = talRandom::giveRandomNumberBetweenZeroAndEntry(rates.first+ rates.second);
			if (rand < rates.first)
				speciateTaxa(bin);
			else extinctTaxa(bin);
			if (_curTaxaNum == 0)
				break;
		}
		if (_curTaxaNum == 0) {
			init(startDiploidNum);
            continue;
		}
		printResults(outfile);
		return;
	}

	errorMsg::reportError("simulatePPInPopulation::run could not reach requested population size with the given speciation and extinction rates. The reason is unknown.");
	return;
}

void simulatePPInPopulation::speciateTaxa(int bin)
{
	//decide if pp speciation or homoploid
	MDOUBLE rand = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
	if (rand < _ppSpeciation) {
		//pp speciation
		if (bin == _ppDist.size() -1) {
			//the last bin represents equal or bigger
			++_ppDist[_ppDist.size()-1]; 
		}
		else
			++_ppDist[bin+1];
	}
	else {
		//homoplid speciation
		++_ppDist[bin];
	}
	++_curTaxaNum;
}


void simulatePPInPopulation::extinctTaxa(int bin)
{
	if (_ppDist[bin] < 1)
		errorMsg::reportError("bin is empty in simulatePPInPopulation::extinctTaxa");
	else
		--_ppDist[bin];
	--_curTaxaNum;
}


//draw a bin from which an individiual will be speciate or extinct.
int simulatePPInPopulation::drawBin()
{
	for (int loop = 0 ; loop < 100000 ; ++loop) 
	{
		int theRandNum = 1 + talRandom::giveIntRandomNumberBetweenZeroAndEntry(_curTaxaNum); //adding 1 since draw does not include last int
		int sum = 0;
		for (int ppLevel = 0; ppLevel < _ppDist.size(); ++ppLevel) 
		{
			sum += _ppDist[ppLevel];
			if (theRandNum <= sum) 
				return ppLevel;
		}
	}
	errorMsg::reportError("simulatePPInPopulation::drawBin could not give a random bin. The reason is unknown.");
	return -1;
}


void simulatePPInPopulation::printResults(const string& outfile)
{
	ofstream outF(outfile.c_str());
	//calc average
	outF<<"###"<<endl;
	MDOUBLE sum = 0.0;
	for (int i = 0; i < _ppDist.size(); ++i) {
		outF<<i+1<<" = "<<_ppDist[i]<<endl;
		sum += (i+1) * _ppDist[i];
	}
	_avg = sum / _curTaxaNum;
	outF<<"#AVG = "<<_avg<<endl;
	outF.close();
}

void simulatePPInPopulation::readDiversificationFile(const string& diversiifcationFile)
{
	_relativeExtRate.clear();
	_speciesRichnessDist.clear();
	if (diversiifcationFile.empty())
		return;
	ifstream inFile(diversiifcationFile.c_str());
	vector<string> data;
	putFileIntoVectorStringArray(inFile, data);
	inFile.close();
	if (data.empty()){
		errorMsg::reportError("unable to open file, or file is empty in readDiversificationFile");
	}
	MDOUBLE sum = 0;
	vector<string>::const_iterator it;
	for (it = data.begin(); it != data.end(); ++it) {
		if (it->empty()) {++it;continue; }// empty line continue
		vector<string> strVec;
		splitString(*it, strVec, " \t"); 
		if (strVec.size() != 4)
			errorMsg::reportError("error reading line: " + *it);
		
		MDOUBLE sr = string2double(trim(strVec[1]).c_str());
		MDOUBLE rate = string2double(trim(strVec[3]));
		_speciesRichnessDist.push_back(sr);
		_relativeExtRate.push_back(rate);
		sum +=sr;
	}
	MDOUBLE sumFreq = 0.0;
	for (int i = 0; i < _speciesRichnessDist.size(); ++i) {
		_speciesRichnessDist[i] /= sum;
		MDOUBLE x = _speciesRichnessDist[i];
		sumFreq += _speciesRichnessDist[i];
	}
	if (!DEQUAL(sumFreq, 1.0, 0.0001))
		errorMsg::reportError("sum frequencies should be 1.0");
}

//in case diversification distribution is given - decide which rates to use. Otherwise - return _specRate + _extRate 
pair<MDOUBLE, MDOUBLE> simulatePPInPopulation::getDiversificationRates()
{
	if (_speciesRichnessDist.empty()) {
		return pair<MDOUBLE, MDOUBLE>(_specRate, _extRate);
	}
	for (int loop = 0 ; loop < 100000 ; ++loop) 
	{
		MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0); 
		MDOUBLE sum = 0.0;
		for (int i = 0; i < _ppDist.size(); ++i) 
		{
			sum += _speciesRichnessDist[i];
			if (theRandNum <= sum) 
				return pair<MDOUBLE, MDOUBLE>(1.0, _relativeExtRate[i]);
		}
	}
	errorMsg::reportError("simulatePPInPopulation::getDiversificationRates could not diversiifcation rates. The reason is unknown.");
	return pair<MDOUBLE, MDOUBLE>(-1.0, -1.0);


}
