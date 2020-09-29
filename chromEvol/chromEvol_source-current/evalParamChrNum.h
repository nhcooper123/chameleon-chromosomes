#ifndef __EVAL_PARAMS_CHR_NUM
#define __EVAL_PARAMS_CHR_NUM

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "chrNumModel.h"

class evalParamChrNum {
public:
	explicit evalParamChrNum(stochasticProcess* pSp, const tree& tr, const sequenceContainer &sc, paramType type, chrNumModel::rootFreqType freqType);
	MDOUBLE operator() (MDOUBLE x);

private:
	const tree& _tr;
	const sequenceContainer& _sc;	
	stochasticProcess* _pSp;
	paramType _paramType;
	chrNumModel::rootFreqType _freqType;
};

class evalParamSpeciationalProp{
public:
	explicit evalParamSpeciationalProp(stochasticProcess* pSp, const tree& originalTree, const sequenceContainer &sc, paramType type, chrNumModel::rootFreqType freqType);
	MDOUBLE operator() (MDOUBLE x);

private:
	const tree& _originalTree;
	const sequenceContainer& _sc;	
	stochasticProcess* _pSp;
	paramType _paramType;
	chrNumModel::rootFreqType _freqType;
};

#endif
