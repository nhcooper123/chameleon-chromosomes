#ifndef ___CHROMOSOME_NUMBER_MODEL
#define ___CHROMOSOME_NUMBER_MODEL

#include "definitions.h"
#include "replacementModel.h"
#include "errorMsg.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "finiteIntAlphabet.h"
class stochasticProcess;


/*
This is the abstract class for a model that describe the change in chromosome number along evolution.
In the derived classes updateQ must be defined.
The model allows at one time step for one gain, one loss, a duplication of the whole genome. 
In the derived classes the the rates of gain,loss, duplication 
may depend (i.e., a linear dependence) or not on the current number of chromosomes.
THE MODEL IS NOT REVERSIBLE!
*/
#define lowerValueOfRateParam 0.00
#define upperValueOfRateParam 100.0
#define lowerValueOfDistrParam 0.05
#define upperValueOfDistrParam 100.00
#define upperValueOfLinearRateParam 5.00
#define upperValueOfEpsRateParam 0.50
#define upperValueOfScaleBranchParam 1.0
//#define lowerBoundBaseNumber 3
#define DEMI_EQUAL_DUPL -2
#define IGNORE_PARAM -999

enum paramType {LOSS_CONST, GAIN_CONST, LOSS_LINEAR, GAIN_LINEAR, DUPL, DUPL_LINEAR, HALF_DUPL, BASE_NUMBER_R, BASE_NUMBER, SCALE_BRANCH};

class paramObj {
public:
	paramObj(paramType type, MDOUBLE value) {_type = type; _value = value;}
	paramObj() {_value = -1000;}
	paramType _type;
	MDOUBLE _value;
};

class chrNumModel : public replacementModel {
public:
	enum modelType {CONST_RATE_NO_DUPL, LINEAR_RATE, LINEAR_RATE_NO_DUPL, CONST_RATE_TRUNCATE, CONST_RATE_BASE_NUMBER, BASE_NUMBER_NO_DUPL, CONST_RATE_EPS, CONST_RATE_EPS_NO_DUPL, GENERAL_CHR_MODEL}; 
	enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED}; //The root frequencies: UNIFORM: all counts in the range [1,maxCromosome] are equally likely.
													//ROOT_LL: the frequencies are according to the likelihood at the root: f(i) = L(i)/SIGMA_i{L(i)}
													//FIXED: the root frequencies are pre-determined from a user input file
	enum jumpType {GAIN_J=0, LOSS_J, DUPL_J, DEMI_J, BASE_J, MAX_CHR_J, JUMP_TYPE_MAX}; //ENUM_LIST_MAX is used just for allocation purposes
	//branchModelType: SPECIATIONAL: all branch lengths are equal; 
	//GRADUAL: change in proporion to given branch lengths; 
	//COMBINED: add a parameter to account for the relative contribution of cladogenesis to the amount of character change;
	//COMBINED_INTERNALS: same as COMBINED but change only internal branches. This makes sense if changes occur at or just after speciation
	enum branchModelType {SPECIATIONAL, GRADUAL, COMBINED, COMBINED_INTERNALS}; 
public:
	explicit chrNumModel(const finiteIntAlphabet* pAlph, const Vdouble& freq, rootFreqType freqType, MDOUBLE halfDuplR);
	explicit chrNumModel();
	
	chrNumModel(const chrNumModel& other);
	virtual ~chrNumModel();
	virtual chrNumModel& operator=(const chrNumModel &other);
	virtual replacementModel* clone() const = 0;
	//virtual chrNumModel* cloneRandomModel() const = 0;
	//set the rate parameters to random values
	virtual void setRandomRateParams(bool bUpdateQ) = 0;

	const int alphabetSize() const {return _freq.size();} 
	const MDOUBLE err_allow_for_pijt_function() const {return 1e-4;} // same as q2p definitions
	const MDOUBLE freq(const int i) const {
		if (i >= _freq.size()) 
			errorMsg::reportError("Error in chrNumModel::freq, i > size of frequency vector");
		return _freq[i];
	}
	void setFreq(int let, MDOUBLE val) {_freq[let] = val;}
	
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE Pij_t2(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;
	
	virtual void printModelParams(ostream &out) const; 
	void printQmatrix(ostream &out) const;
	virtual void getParams(vector<paramObj>& params) const = 0; 
	virtual void getOptimizedParams(vector<paramObj>& params) const = 0; 
	virtual void setParam(MDOUBLE val, paramType type); 
	virtual string getParamName(paramType type) const; 
	virtual MDOUBLE getParamLowerBound(paramType type) const; 
	virtual MDOUBLE getParamUpperBound(paramType type) const;

	virtual bool setGainConstR(MDOUBLE inVal) = 0; //returns true if a change was done and false if the current val is equal to the inVal
	virtual MDOUBLE getGainConstR() const = 0;
	virtual bool setLossConstR(MDOUBLE inVal) = 0;
	virtual MDOUBLE getLossConstR() const = 0;
	virtual bool setDuplConstR(MDOUBLE inVal) = 0; 
	virtual MDOUBLE getDuplConstR() const = 0;
	virtual bool setLossLinearR(MDOUBLE inVal); 
	virtual bool setGainLinearR(MDOUBLE inVal); 
	virtual bool setBaseNumberR(MDOUBLE inVal) = 0; 
	virtual MDOUBLE getBaseNumberR() const = 0;
	virtual bool setBaseNumber(int inVal) = 0; 
	virtual int getBaseNumber() const = 0;

	bool setHalfDuplR(MDOUBLE inVal);
	virtual MDOUBLE getScaleFactor() {return _scaleFactor;}
	virtual bool setScaleFactor(MDOUBLE inVal);
	
	////get rate parameters for a specific chromosome. 
	//e.g: using a linear model the total gain param is gain_const + fromState*gain_linear
	virtual MDOUBLE getGainR(int fromState) const {errorMsg::reportError("getGainR() not implemented"); return 0.0;};
	virtual MDOUBLE getLossR(int fromState) const {errorMsg::reportError("getLossR() not implemented");return 0.0;};
	virtual MDOUBLE getDuplR(int fromState) const {errorMsg::reportError("getDuplR() not implemented");return 0.0;};
	virtual MDOUBLE getDemiR(int fromState) const {errorMsg::reportError("getDemiR() not implemented");return 0.0;};
	
	virtual void copy(const chrNumModel* other);

	void mutiplyQMatrixByScalar(MDOUBLE sc);
	MDOUBLE getNormQ () const {return _normQ;};
protected:
	virtual void init();
	virtual void validateParams() = 0;
	virtual void updateQ() = 0;
	//update only those cells of the the Q matrix that are relavant to "half duplications": i->1.5i
	virtual void updateQHalfDupl(MDOUBLE halfDuplRate, MDOUBLE duplRate);
	virtual bool pijt_is_prob_value(MDOUBLE val) const;

	//calculate exp(inMat*t)
	void exp_Mt(const VVdouble& inMat, MDOUBLE t, VVdoubleRep& outMat, bool bScaleSquare) const;
	
protected:
	const finiteIntAlphabet* _pAlph;
	Vdouble _freq;
	VVdouble _Q;
	rootFreqType _rootFreqType;
	MDOUBLE _normQ; //the first norm of the Q matrix

	MDOUBLE _halfDuplR; //specify i->1.5i transitions. If ==-1 then ignore. If ==DEMI_EQUAL_DUPL then equal to _duplConstR
	MDOUBLE _scaleFactor; //Used only in COMBINED branch model: add a fixed length to all branches 


	mutable bool _bQchanged; //indicates whether the Q matrix was changed after the last Pij_t call
	mutable VVdoubleRep _lastPtCalculated;
	mutable MDOUBLE _lastTcalculated;
};

#endif

