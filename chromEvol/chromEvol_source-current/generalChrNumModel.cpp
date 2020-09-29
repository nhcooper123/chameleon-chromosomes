#include "generalChrNumModel.h"
#include "matrixUtils.h"
#include "someUtil.h"
#include "talRandom.h"
#include "chrNumberOptions.h"
#include <algorithm>



generalChrNumModel::generalChrNumModel(const finiteIntAlphabet* pAlph, const Vdouble& freq, rootFreqType freqType, vector<paramObj>& initParams, int maxBaseTransition)
: chrNumModel(pAlph, freq, freqType, IGNORE_PARAM), _maxBaseTransition(maxBaseTransition) 
{
	init(initParams);
	updateQ();
	LOGDO(10,printMatrix(_Q, myLog::LogFile()));
}


generalChrNumModel::generalChrNumModel(const finiteIntAlphabet* pAlph, const Vdouble& freq, rootFreqType freqType, vector<paramObj>& initParams, const vector<pair<int, MDOUBLE> >& compProbs)
: chrNumModel(pAlph, freq, freqType, IGNORE_PARAM), _baseTranProbs(compProbs)
{
	_maxBaseTransition = _baseTranProbs[compProbs.size()-1].first; 
	init(initParams);
	updateQ();
	LOGDO(10,printMatrix(_Q, myLog::LogFile()));
}


generalChrNumModel::generalChrNumModel(const generalChrNumModel& other)
{
	copy(&other);
}

generalChrNumModel::~generalChrNumModel()
{
}

generalChrNumModel& generalChrNumModel::operator=(const generalChrNumModel&other)
{
	copy(&other);
	return *this;
}

void generalChrNumModel::copy(const generalChrNumModel* pOther)
{
	_gainConstR = pOther->_gainConstR;
	_gainLinearR = pOther->_gainLinearR;
	_lossConstR = pOther->_lossConstR;
	_lossLinearR = pOther->_lossLinearR;
	_duplConstR = pOther->_duplConstR;
	_duplLinearR = pOther->_duplLinearR;
	_baseNumberR = pOther->_baseNumberR;
	_baseNumber = pOther->_baseNumber;
	_validParams = pOther->_validParams;
	_maxBaseTransition = pOther->_maxBaseTransition;
	_baseTranProbs = pOther->_baseTranProbs;
	chrNumModel::copy(pOther);
}

void generalChrNumModel::init(const vector<paramObj>& initParams)
{
	_gainConstR = IGNORE_PARAM; 
	_gainLinearR = IGNORE_PARAM; 
	_lossConstR = IGNORE_PARAM; 
	_lossLinearR = IGNORE_PARAM; 
	_duplConstR = IGNORE_PARAM; 
	_duplLinearR = IGNORE_PARAM;
	_halfDuplR = IGNORE_PARAM;
	_baseNumberR = IGNORE_PARAM;
	_baseNumber = IGNORE_PARAM;
	_scaleFactor = IGNORE_PARAM;
	bool bValidateParam = false;
	bool bDemiEqualDup = false;
	for (int i = 0; i < initParams.size(); ++i) 
	{
		if (!DEQUAL(initParams[i]._value,IGNORE_PARAM))
		{
			_validParams.push_back(initParams[i]);
			bValidateParam = true;
			switch (initParams[i]._type)
			{
				case LOSS_CONST:
					if ((initParams[i]._value < getParamUpperBound(LOSS_CONST)) && (initParams[i]._value > getParamLowerBound(LOSS_CONST)))
						_lossConstR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _lossConstR is out of bound");
					break;
				case LOSS_LINEAR:
					if ((initParams[i]._value < getParamUpperBound(LOSS_LINEAR)) && (initParams[i]._value > getParamLowerBound(LOSS_LINEAR)))
						_lossLinearR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _lossLinearR is out of bound");
					break;
				case GAIN_CONST: 
					if ((initParams[i]._value < getParamUpperBound(GAIN_CONST)) && (initParams[i]._value > getParamLowerBound(GAIN_CONST)))
						_gainConstR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _gainConstR is out of bound");
					break;
				case GAIN_LINEAR: 
					if ((initParams[i]._value < getParamUpperBound(GAIN_LINEAR)) && (initParams[i]._value > getParamLowerBound(GAIN_LINEAR)))
						_gainLinearR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _gainLinearR is out of bound");
					break;
				case DUPL: 
					if ((initParams[i]._value < getParamUpperBound(DUPL)) && (initParams[i]._value > getParamLowerBound(DUPL)))
						_duplConstR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _duplConstR is out of bound");
					break;
				case DUPL_LINEAR: 
					if ((initParams[i]._value < getParamUpperBound(DUPL_LINEAR)) && (initParams[i]._value > getParamLowerBound(DUPL_LINEAR)))
						_duplLinearR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _duplLinearR is out of bound");
					break;
				case HALF_DUPL: 
					if (initParams[i]._value == DEMI_EQUAL_DUPL) {
						bDemiEqualDup = true; 
						continue;
					}
					if ((initParams[i]._value < getParamUpperBound(HALF_DUPL)) && (initParams[i]._value > getParamLowerBound(HALF_DUPL)))
						_halfDuplR = initParams[i]._value;
					else
						errorMsg::reportError("error in generalChrNumModel::init _halfDuplR is out of bound");
					break;
				case BASE_NUMBER_R: 
					if ((initParams[i]._value < getParamUpperBound(BASE_NUMBER_R)) && (initParams[i]._value > getParamLowerBound(BASE_NUMBER_R)))
						_baseNumberR = initParams[i]._value;
					else 
						errorMsg::reportError("error in generalChrNumModel::init _baseNumberR is out of bound");
					break;
				case BASE_NUMBER: 
					if ((initParams[i]._value <= getParamUpperBound(BASE_NUMBER)) && (initParams[i]._value >= getParamLowerBound(BASE_NUMBER)))
						_baseNumber = static_cast<int>(initParams[i]._value);
					else 
						errorMsg::reportError("error in generalChrNumModel::init _baseNumber is out of bound");
					break;
				case SCALE_BRANCH: 
					_scaleFactor = initParams[i]._value;
					break;
				default: 
					errorMsg::reportError("error in generalChrNumModel::init unknown param type");
					break;
			}
		}
	}
	if (bDemiEqualDup)
		_halfDuplR = DEMI_EQUAL_DUPL;

	if (!bValidateParam)
		errorMsg::reportError("error in generalChrNumModel::init there are no parameters to be included in the model");
	if ((_halfDuplR != IGNORE_PARAM) && (_duplLinearR != IGNORE_PARAM))
		errorMsg::reportError("not implemented yet: duplication rate linearly dependent on the current count with demiDuplication");
	if ((_baseNumber != IGNORE_PARAM) && (_baseNumberR == IGNORE_PARAM))
		errorMsg::reportError("cannot use a model with a baseNumber without using the _baseNumberR rate parameter");
	if ((_baseNumber == IGNORE_PARAM) && (_baseNumberR != IGNORE_PARAM))
		errorMsg::reportError("cannot use a model with a rate for baseNumber without using the _baseNumber parameter");
	if ((_baseNumber != IGNORE_PARAM) && (_maxBaseTransition < _baseNumber))
		errorMsg::reportError("cannot use a model with a baseNumber while the maximum base transition is lower than the base number");

}

//reset the rate parameters to random values
void generalChrNumModel::setRandomRateParams(bool bUpdateQ)
{
	//note: _gainLinearR and _lossLinearR can be negative but larger than -_gainConst/maxChr so that all entries in the Q are positive
	//in such a case _gainConstR is also restricted such that the lower bound is max{0, g_l(1-maxNumber)} 
	if (!DEQUAL(_gainConstR, IGNORE_PARAM))
		_gainConstR = talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(GAIN_CONST), getParamUpperBound(GAIN_CONST));
	if (!DEQUAL(_gainLinearR, IGNORE_PARAM))
        _gainLinearR = talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(GAIN_LINEAR), getParamUpperBound(GAIN_LINEAR));
	if (!DEQUAL(_lossConstR, IGNORE_PARAM))
        _lossConstR = talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(LOSS_CONST), getParamUpperBound(LOSS_CONST));
	if (!DEQUAL(_lossLinearR, IGNORE_PARAM))
        _lossLinearR = talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(LOSS_LINEAR), getParamUpperBound(LOSS_LINEAR));
	if (!DEQUAL(_duplConstR, IGNORE_PARAM))
        _duplConstR = talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(DUPL), getParamUpperBound(DUPL));
	if (!DEQUAL(_baseNumberR, IGNORE_PARAM))
        _baseNumberR = talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(BASE_NUMBER_R), getParamUpperBound(BASE_NUMBER_R));
	if (!DEQUAL(_baseNumber, IGNORE_PARAM) && chrNumberOptions::_bOptBaseNumber !=0)
        _baseNumber = static_cast<int>(talRandom::giveRandomNumberBetweenTwoPoints(getParamLowerBound(BASE_NUMBER), getParamUpperBound(BASE_NUMBER)));
		

	chrNumModel::setRandomRateParams(false);
	if (bUpdateQ == true)
		updateQ();
}

MDOUBLE generalChrNumModel::getParamUpperBound(paramType type) const
{
	MDOUBLE res = 0.0;
	switch (type)
	{
	case GAIN_LINEAR:
	case LOSS_LINEAR:
	case DUPL_LINEAR:
		res = upperValueOfLinearRateParam;
		break;
	case GAIN_CONST:
	case LOSS_CONST:
	case DUPL:
	case BASE_NUMBER_R:
		res = upperValueOfRateParam;
		break;
	case BASE_NUMBER:
		res = _maxBaseTransition;
		break;
	default:
		res = chrNumModel::getParamUpperBound(type);
	}
	return res;
}


MDOUBLE generalChrNumModel::getParamLowerBound(paramType type) const
{
	MDOUBLE res = 0.0;
	switch (type)
	{
	case GAIN_LINEAR:
		if (DEQUAL(_gainConstR, IGNORE_PARAM))
			res = lowerValueOfRateParam;
		else
			res = -_gainConstR / (_pAlph->max()-1);
		break;
	case LOSS_LINEAR:
		if (DEQUAL(_lossConstR, IGNORE_PARAM))
			res = lowerValueOfRateParam;
		else
			res = -_lossConstR / (_pAlph->max()-1);
		break;
	case DUPL_LINEAR:
		if (DEQUAL(_duplConstR, IGNORE_PARAM))
			res = lowerValueOfRateParam;
		else
			res = -_duplConstR / (_pAlph->max()-1);
		break;
	case LOSS_CONST:
		if (DEQUAL(_lossLinearR, IGNORE_PARAM))
            res = lowerValueOfRateParam;
		else
			res = (_lossLinearR > 0) ? lowerValueOfRateParam : (_lossLinearR * (1-_pAlph->max()));
		break;
	case GAIN_CONST:
		if (DEQUAL(_gainLinearR, IGNORE_PARAM))
            res = lowerValueOfRateParam;
		else
			res = (_gainLinearR > 0) ? lowerValueOfRateParam : (_gainLinearR * (1-_pAlph->max()));
        break;
	case DUPL:
		if (DEQUAL(_duplLinearR, IGNORE_PARAM))
            res = lowerValueOfRateParam;
		else
			res = (_duplLinearR > 0) ? lowerValueOfRateParam : (_duplLinearR * (1-_pAlph->max()));
        break;
	case BASE_NUMBER_R:
		res = lowerValueOfRateParam;
		break;
	case BASE_NUMBER:
		res = chrNumberOptions::_minBaseTransition;
		break;
	default:
		res = chrNumModel::getParamLowerBound(type);
	}
	return res;
}

bool generalChrNumModel::setGainConstR(MDOUBLE val) 
{ 
	if (val == _gainConstR)
		return false;
	if ((val <= getParamLowerBound(GAIN_CONST)) || (val >= getParamUpperBound(GAIN_CONST))) {
		errorMsg::reportError("in generalChrNumModel::setLossConstR val is out of bound");
		return false;
	}

    _gainConstR = val; 
	updateQ();
	return true;
}


MDOUBLE generalChrNumModel::getGainR(int fromState) const
{
	MDOUBLE gainLinearR = DEQUAL(_gainLinearR, IGNORE_PARAM) ? 0.0 : _gainLinearR;
	MDOUBLE gainConstR = DEQUAL(_gainConstR, IGNORE_PARAM) ? 0.0 : _gainConstR;
	return gainConstR + fromState*gainLinearR;
}

MDOUBLE generalChrNumModel::getLossR(int fromState) const
{
    MDOUBLE lossLinearR = DEQUAL(_lossLinearR, IGNORE_PARAM) ? 0.0 : _lossLinearR;
	MDOUBLE lossConstR = DEQUAL(_lossConstR, IGNORE_PARAM) ? 0.0 : _lossConstR;
	return lossConstR + fromState*lossLinearR;
}

MDOUBLE generalChrNumModel::getDuplR(int fromState) const
{
 	MDOUBLE duplLinearR = DEQUAL(_duplLinearR, IGNORE_PARAM) ? 0.0 : _duplLinearR;
	MDOUBLE duplConstR = DEQUAL(_duplConstR, IGNORE_PARAM) ? 0.0 : _duplConstR;
	return duplConstR + fromState*duplLinearR;
}

MDOUBLE generalChrNumModel::getDemiR(int fromState) const
{
	if (DEQUAL(_halfDuplR, IGNORE_PARAM))
		return 0;
	else if (_halfDuplR == DEMI_EQUAL_DUPL)
		return getDuplR(fromState);
	else
        return _halfDuplR;
}


bool generalChrNumModel::setLossConstR(MDOUBLE val) 
{ 
	if (val == _lossConstR)
		return false;
	if ((val <= getParamLowerBound(LOSS_CONST)) || (val >= getParamUpperBound(LOSS_CONST))) {
		errorMsg::reportError("in generalChrNumModel::setLossConstR val is out of bound");
		return false;
	}
    _lossConstR= val; 
	updateQ();
	return true;
}

bool generalChrNumModel::setDuplConstR(MDOUBLE val) 
{ 
	if (val == _duplConstR)
		return false;
	if ((val <= getParamLowerBound(DUPL)) || (val >= getParamUpperBound(DUPL))) {
		errorMsg::reportError("in generalChrNumModel::setDuplConstR val is out of bound");
		return false;
	}

	_duplConstR= val; 
	updateQ();
	return true;
}

bool generalChrNumModel::setDuplLinearR(MDOUBLE val)
{
	if (val == _duplLinearR)
		return false;
	if ((val <= getParamLowerBound(DUPL_LINEAR)) || (val >= getParamUpperBound(DUPL_LINEAR))) {
		errorMsg::reportError("in generalChrNumModel::setDuplLinearR val is out of bound");
		return false;
	}
    _duplLinearR= val; 
	updateQ();
	return true;
}
bool generalChrNumModel::setLossLinearR(MDOUBLE val) 
{
	if (val == _lossLinearR)
		return false;
	if ((val <= getParamLowerBound(LOSS_LINEAR)) || (val >= getParamUpperBound(LOSS_LINEAR))) {
		errorMsg::reportError("in generalChrNumModel::setLossLinearR val is out of bound");
		return false;
	}
    _lossLinearR = val; 
	updateQ();
	return true;
}

bool generalChrNumModel::setGainLinearR(MDOUBLE val) 
{
	if (val == _gainLinearR)
		return false;
	if ((val <= getParamLowerBound(GAIN_LINEAR)) || (val >= getParamUpperBound(GAIN_LINEAR))) {
		errorMsg::reportError("in generalChrNumModel::setGainLinearR val is out of bound");
		return false;
	}
    _gainLinearR = val; 
	updateQ();
	return true;
}

bool generalChrNumModel::setBaseNumberR(MDOUBLE val) 
{
	if (val == _baseNumberR)
		return false;
	if ((val <= getParamLowerBound(BASE_NUMBER_R)) || (val >= getParamUpperBound(BASE_NUMBER_R))) {
		errorMsg::reportError("in generalChrNumModel::setBaseNumberR val is out of bound");
		return false;
	}
    _baseNumberR = val; 
	updateQ();
	return true;
}

bool generalChrNumModel::setBaseNumber(int val) 
{
	if (val == _baseNumber)
		return false;
	if ((val < getParamLowerBound(BASE_NUMBER)) || (val > getParamUpperBound(BASE_NUMBER))) {
		errorMsg::reportError("in generalChrNumModel::setBaseNumber val is out of bound");
		return false;
	}
    _baseNumber = val; 
	updateQ();
	return true;
}


void generalChrNumModel::getParams(vector<paramObj>& params) const
{
	params.clear();
	if (!DEQUAL(IGNORE_PARAM, _lossConstR)) {
		params.push_back(paramObj(LOSS_CONST, _lossConstR));
	}
	if (!DEQUAL(IGNORE_PARAM, _gainConstR)) {
		params.push_back(paramObj(GAIN_CONST, _gainConstR));
	}
	if (!DEQUAL(IGNORE_PARAM, _duplConstR)) {
		params.push_back(paramObj(DUPL, _duplConstR));
	}
	if (!DEQUAL(IGNORE_PARAM, _lossLinearR)) {
		params.push_back(paramObj(LOSS_LINEAR, _lossLinearR));
	}
	if (!DEQUAL(IGNORE_PARAM, _gainLinearR)) {
		params.push_back(paramObj(GAIN_LINEAR, _gainLinearR));
	}
	if (!DEQUAL(IGNORE_PARAM, _duplLinearR)) {
		params.push_back(paramObj(DUPL_LINEAR, _duplLinearR));
	}
	if (!DEQUAL(IGNORE_PARAM, _baseNumberR)) {
		params.push_back(paramObj(BASE_NUMBER_R, _baseNumberR));
	}
	if (!DEQUAL(IGNORE_PARAM, _baseNumber)) {
		params.push_back(paramObj(BASE_NUMBER, _baseNumber));
	}

	chrNumModel::getParams(params);
}



void generalChrNumModel::getOptimizedParams(vector<paramObj>& params) const
{
	params.clear();
	if (!DEQUAL(IGNORE_PARAM, _baseNumber) && chrNumberOptions::_bOptBaseNumber) {
		params.push_back(paramObj(BASE_NUMBER, _baseNumber));
	}
	if (!DEQUAL(IGNORE_PARAM, _baseNumberR)) {
		params.push_back(paramObj(BASE_NUMBER_R, _baseNumberR));
	}
	if (!DEQUAL(IGNORE_PARAM, _duplConstR)) {
		params.push_back(paramObj(DUPL, _duplConstR));
	}
	if (!DEQUAL(IGNORE_PARAM, _lossConstR)) {
		params.push_back(paramObj(LOSS_CONST, _lossConstR));
	}
	if (!DEQUAL(IGNORE_PARAM, _gainConstR)) {
		params.push_back(paramObj(GAIN_CONST, _gainConstR));
	}
	if (!DEQUAL(IGNORE_PARAM, _lossLinearR)) {
		params.push_back(paramObj(LOSS_LINEAR, _lossLinearR));
	}
	if (!DEQUAL(IGNORE_PARAM, _gainLinearR)) {
		params.push_back(paramObj(GAIN_LINEAR, _gainLinearR));
	}
	if (!DEQUAL(IGNORE_PARAM, _duplLinearR)) {
		params.push_back(paramObj(DUPL_LINEAR, _duplLinearR));
	}

	chrNumModel::getOptimizedParams(params);
}


//void generalChrNumModel::updateQ()
//{
//	zeroMatrix(_Q);
//	//use local rates instead of asking if should IGNORE every time;
//    MDOUBLE lossLinearR = DEQUAL(_lossLinearR, IGNORE_PARAM) ? 0.0 : _lossLinearR;
//	MDOUBLE gainLinearR = DEQUAL(_gainLinearR, IGNORE_PARAM) ? 0.0 : _gainLinearR;
//	MDOUBLE duplLinearR = DEQUAL(_duplLinearR, IGNORE_PARAM) ? 0.0 : _duplLinearR;
//	MDOUBLE lossConstR = DEQUAL(_lossConstR, IGNORE_PARAM) ? 0.0 : _lossConstR;
//	MDOUBLE gainConstR = DEQUAL(_gainConstR, IGNORE_PARAM) ? 0.0 : _gainConstR;
//	MDOUBLE duplConstR = DEQUAL(_duplConstR, IGNORE_PARAM) ? 0.0 : _duplConstR;
//
//
//	//update first row (chrNum = 1): no loss possible. Q[1][2] = (g_c)+d
//	_Q[0][0] = - (gainConstR + duplConstR);
//	_Q[0][1] = gainConstR + duplConstR;
//
//    //update upper half of the matrix
//	int chr = 2;
//	MDOUBLE lossR, gainR, duplR;
//	for (; chr <= alphabetSize()/2; ++chr)
//	{
//		lossR = lossConstR + ((chr-1) * lossLinearR);
//		gainR = gainConstR + ((chr-1) * gainLinearR);
//		duplR = duplConstR + ((chr-1) * duplLinearR); 
//		_Q[chr-1][chr-2] = lossR;
//		_Q[chr-1][chr-1] = -lossR - gainR - duplR;
//		_Q[chr-1][chr] = gainR;
//		_Q[chr-1][2*chr -1] = duplR;
//	}
//	//update lower half of the matrix: duplication always go to last column
//	for (; chr < alphabetSize()-1; ++chr)
//	{
//		lossR = lossConstR + ((chr-1) * lossLinearR);
//		gainR = gainConstR + ((chr-1) * gainLinearR);
//		duplR = duplConstR + ((chr-1) * duplLinearR); 
//		_Q[chr-1][chr-2] = lossR;
//		_Q[chr-1][chr-1] = -lossR - gainR - duplR;
//		_Q[chr-1][chr] = gainR;
//		_Q[chr-1][alphabetSize()-1] = duplR;
//	}
//	//update one row before last: moving up is by gains+duplication 
//	lossR = lossConstR + ((chr-1) * lossLinearR);
//	gainR = gainConstR + ((chr-1) * gainLinearR);
//	duplR = duplConstR + ((chr-1) * duplLinearR);
//	_Q[chr-1][chr-2] = lossR;
//	_Q[chr-1][chr-1] = -lossR - gainR - duplR;
//	_Q[chr-1][chr] = gainR + duplR;
//
//	++chr;
//	lossR = lossConstR + ((chr-1) * lossLinearR);
//	gainR = gainConstR + ((chr-1) * gainLinearR);
//	duplR = duplConstR + ((chr-1) * duplLinearR); 
//
//	//update last row: no gains+duplication 
//	_Q[chr-1][chr-2] = lossR;
//	_Q[chr-1][chr-1] = -lossR;
//	if (_halfDuplR != IGNORE_PARAM)
//	{
//		if (_halfDuplR == DEMI_EQUAL_DUPL)
//            updateQHalfDupl(_duplConstR, _duplConstR);
//		else
//			updateQHalfDupl(_halfDuplR, _duplConstR);
//	}
//	_bQchanged = true;
//	_normQ = getMatrixNorm(_Q);
//}


/*
In case we want to allow for addition of a base number, so that
additions of _baseNumber*c chromosome are allowed (c = integer--> any multiplication is also allowed)
we have to add 2 parameters: the base number and its transition rate. In this case 
Q_ij = 
		g if j=i+1
		l if j=i-1
		d if j=i*2
		b if (j-i)%baseNumber == 0

If both a duplication (or demi duplication) and addition of the base number are possible then the rate is max{d,b, demi}
Have to think if we want an additional parameter for the baseRate

*/

void generalChrNumModel::updateQ()
{
	zeroMatrix(_Q);
	//use local rates instead of asking if should IGNORE every time;
    MDOUBLE lossLinearR = DEQUAL(_lossLinearR, IGNORE_PARAM) ? 0.0 : _lossLinearR;
	MDOUBLE gainLinearR = DEQUAL(_gainLinearR, IGNORE_PARAM) ? 0.0 : _gainLinearR;
	MDOUBLE duplLinearR = DEQUAL(_duplLinearR, IGNORE_PARAM) ? 0.0 : _duplLinearR;
	MDOUBLE lossConstR = DEQUAL(_lossConstR, IGNORE_PARAM) ? 0.0 : _lossConstR;
	MDOUBLE gainConstR = DEQUAL(_gainConstR, IGNORE_PARAM) ? 0.0 : _gainConstR;
	MDOUBLE duplConstR = DEQUAL(_duplConstR, IGNORE_PARAM) ? 0.0 : _duplConstR;
		
	int maxCount = _pAlph->max();
    //update upper half of the matrix
	MDOUBLE lossR, gainR, duplR;
	int id = 0;
	int count = _pAlph->id2Count(id);
	//update upper half of the matrix
	for (; count <= maxCount/2; ++id,count = _pAlph->id2Count(id))
	{
        lossR = lossConstR + ((count-1) * lossLinearR);
		gainR = gainConstR + ((count-1) * gainLinearR);
		duplR = duplConstR + ((count-1) * duplLinearR); 
		if (id == 0)
		{//update first row (minimum count): no loss possible. Q[id][2id] = (g_c)+d
			_Q[id][id] = - (gainR + duplR);
			_Q[0][id+1] = gainR;
			_Q[0][_pAlph->count2Id(count*2)] += duplR;
		}
		else 
		{
			_Q[id][id-1] = lossR;
			_Q[id][id] = -lossR - gainR - duplR;
			_Q[id][id+1] = gainR;
			_Q[id][_pAlph->count2Id(count*2)] = duplR;
		}
	}
	//update lower half of the matrix: duplication always go to last column
	for (; count <= maxCount-2; ++id,count = _pAlph->id2Count(id))
	{
		lossR = lossConstR + ((count-1) * lossLinearR);
		gainR = gainConstR + ((count-1) * gainLinearR);
		duplR = duplConstR + ((count-1) * duplLinearR); 
		_Q[id][id] = -gainR - duplR;
		_Q[id][id+1] = gainR;
		_Q[id][alphabetSize()-1] = duplR;
		if (id > 0) {
			_Q[id][id-1] = lossR;
			_Q[id][id] -= lossR;
		}
	}
	//update one row before last: moving up is by gains+duplication 
	lossR = lossConstR + ((count-1) * lossLinearR);
	gainR = gainConstR + ((count-1) * gainLinearR);
	duplR = duplConstR + ((count-1) * duplLinearR);
	_Q[id][id-1] = lossR;
	_Q[id][id] = -lossR - gainR - duplR;
	_Q[id][id+1] = gainR + duplR;

	++id;
	count = _pAlph->id2Count(id);
	lossR = lossConstR + ((count-1) * lossLinearR);
	
	//update last row: no gains+duplication 
	_Q[id][id-1] = lossR;
	_Q[id][id] = -lossR;
	if (_halfDuplR != IGNORE_PARAM)
	{
		if (_halfDuplR == DEMI_EQUAL_DUPL)
            updateQHalfDupl(_duplConstR, _duplConstR);
		else
			updateQHalfDupl(_halfDuplR, _duplConstR);
	}
	if ((_baseNumberR != IGNORE_PARAM) && (_baseNumber != IGNORE_PARAM))
	{
		updateQBaseNumber();
	}
	_bQchanged = true;
	_normQ = getMatrixNorm(_Q);
}


//Note that we allow both duplication and baseNumber transition, so in case base=5 then transition 5-->10 will be baseR+duplR
//If we have a vector of possible base transitions (if b=5 this can be P(5 transition)=0.4, P(15 transition)=0.2, P(20 transition)=0.4), 
//then we should first construct a vector of possible base number transitions where each entry is the probability of the jump from (toCol - row). 
//we should then multiply the _baseNumberr by this probability when updating the Q matrix
void generalChrNumModel::updateQBaseNumber()
{
	if (!_baseTranProbs.empty())
		return updateQBaseNumberTransProb();
	int maxCount = _pAlph->max();
	//1. update the last column. 
	int row = 0;
	int count = _pAlph->id2Count(row);
	//update upper half of the matrix
	int lastCol = alphabetSize()-1;
	for (; count < maxCount; ++row,count = _pAlph->id2Count(row))
	{
		if ((maxCount - count) <= _maxBaseTransition) {
			_Q[row][lastCol] += _baseNumberR;
			_Q[row][row] -= _baseNumberR;
		}
	}

	
	row = 0;
	count = _pAlph->id2Count(row);
	for (; count < maxCount; ++row,count = _pAlph->id2Count(row))
	{
		int toCol = row + _baseNumber;
		for (; toCol < lastCol; toCol += _baseNumber)
		{
			if ((toCol - row) > _maxBaseTransition)
				break;
			_Q[row][toCol] += _baseNumberR;
			_Q[row][row] -= _baseNumberR;
		}
	}
}

//Note that we allow both duplication and baseNumber transition, so in case base=5 then transition 5-->10 will be baseR+duplR
//If we have a vector of possible base transitions (if b=5 this can be P(5 transition)=0.4, P(15 transition)=0.2, P(20 transition)=0.4), 
//then we should first construct a vector of possible base number transitions where each entry is the probability of the jump from (toCol - row). 
//we should then multiply the _baseNumberR by this probability when updating the Q matrix. 
//This is needed because the rate extracted from the inference step was applied for all possible baseNumber transitions per row. 
//However, in the simulations we want to restrict the range of possible transitions. 
void generalChrNumModel::updateQBaseNumberTransProb()
{
	int lastCol = alphabetSize()-1;
	int maxCount = _pAlph->max();
	int row = 0;
	int count = _pAlph->id2Count(row);
	for (; count < maxCount; ++row,count = _pAlph->id2Count(row))
	{
		//compute a row-specific baseNumberR - that takes into account the number of possible baseNumber transitions in that row. 
		int maxTransitions = min(_maxBaseTransition, lastCol-row);
		int possibleTransitions = maxTransitions/_baseNumber;
		if ((maxTransitions % _baseNumber) != 0) //in case we have a reminder - there is another 
			++possibleTransitions;
		MDOUBLE baseNumberR = possibleTransitions * _baseNumberR;

		MDOUBLE sum = 0.0;
		for (int i = 0; i < _baseTranProbs.size(); ++i)
		{
			int toCol = row + _baseTranProbs[i].first;
			MDOUBLE prob = _baseTranProbs[i].second;
			if (((toCol - row) > _maxBaseTransition) || (toCol >= maxCount))
				break;
			MDOUBLE rate = baseNumberR * prob;
			_Q[row][toCol] += rate;
			_Q[row][row] -= rate;
			sum += prob;
		}
		//update the last column
			//1. update the last column. 
		MDOUBLE lastColRate = (1.0 - sum) * baseNumberR;
		_Q[row][lastCol] += lastColRate;
		_Q[row][row] -= lastColRate;
	}
}
