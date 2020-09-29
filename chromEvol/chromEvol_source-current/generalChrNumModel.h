#ifndef __GENERAL_CHROMOSOME_NUMBER_MODEL
#define __GENERAL_CHROMOSOME_NUMBER_MODEL

#include "definitions.h"
#include "chrNumModel.h"


/*
This model includes all possible parameters. 
the vector _validParams indicate the parameters that are included in the model. 
Currently all included params are also optimized
*/


class generalChrNumModel : public chrNumModel {
public:
	explicit generalChrNumModel(const finiteIntAlphabet* pAlph, const Vdouble& freq, rootFreqType freqType, vector<paramObj>& initParams, int maxBaseTransition);
	explicit generalChrNumModel(const finiteIntAlphabet* pAlph, const Vdouble& freq, rootFreqType freqType, vector<paramObj>& initParams, const vector<pair<int, MDOUBLE> >& compProbs);
	
	//creates a random model:
	//explicit generalChrNumModel(int maxChrAllowed, const Vdouble& freq, rootFreqType freqType);
	
	generalChrNumModel(const generalChrNumModel& other);
	virtual ~generalChrNumModel();
	virtual generalChrNumModel& operator=(const generalChrNumModel &other);
	virtual replacementModel* clone() const { return new generalChrNumModel(*this); }
	//reset the rate parameters to random values
	virtual void setRandomRateParams(bool bUpdateQ);


	virtual void getParams(vector<paramObj>& params) const; 
	virtual void getOptimizedParams(vector<paramObj>& params) const; 
	virtual bool setGainConstR(MDOUBLE inVal); //returns true if a change was done and false if the current val is equal to the inVal
	MDOUBLE getGainConstR() const {return _gainConstR;} 
	virtual bool setLossConstR(MDOUBLE inVal);
	MDOUBLE getLossConstR() const {return _lossConstR;}
	virtual bool setDuplConstR(MDOUBLE inVal); 
	MDOUBLE getDuplConstR() const {return _duplConstR;}
	virtual bool setBaseNumberR(MDOUBLE inVal); 
	MDOUBLE getBaseNumberR() const { return DEQUAL(_baseNumberR, IGNORE_PARAM) ? 0.0 : _baseNumberR; }
	virtual bool setBaseNumber(int inVal);
	virtual int getBaseNumber() const {return _baseNumber;};

	virtual bool setDuplLinearR(MDOUBLE inVal); 
	virtual bool setLossLinearR(MDOUBLE inVal); 
	virtual bool setGainLinearR(MDOUBLE inVal); 

	////get rate parameters for a specific chromosome. 
	//e.g: using a linear model the total gain param is gain_const + fromState*gain_linear
	virtual MDOUBLE getGainR(int fromState) const;
	virtual MDOUBLE getLossR(int fromState) const;
	virtual MDOUBLE getDuplR(int fromState) const;
	virtual MDOUBLE getDemiR(int fromState) const;
	


	MDOUBLE getParamUpperBound(paramType type) const;
	MDOUBLE getParamLowerBound(paramType type) const;



protected:
	virtual void updateQ();
	virtual void updateQBaseNumber();
	virtual void updateQBaseNumberTransProb();

	virtual void validateParams() {};
	void copy(const generalChrNumModel* pOther);
	void init(const vector<paramObj>& initParams);
	MDOUBLE getParamVal(paramType type);

protected:
	vector<paramObj> _validParams; //specify the parameters that should be included in the model
	//rate parameters:
	//constant = does not depend on the current number of chromosomes
	//linear = rate depends linearly on the current chromosome number
	MDOUBLE _gainConstR; 
	MDOUBLE _gainLinearR; 
	MDOUBLE _lossConstR; 
	MDOUBLE _lossLinearR; 
	MDOUBLE _duplConstR; 
	MDOUBLE _duplLinearR; 
	MDOUBLE _baseNumberR; 
	int _baseNumber; 
	int _maxBaseTransition; //used only when using the base number parameter - the largest base number transition allowed
	vector<pair<int, MDOUBLE> > _baseTranProbs; 
};

#endif
