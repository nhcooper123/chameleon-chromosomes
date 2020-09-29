#include "chrNumModel.h"
#include "matrixUtils.h"
#include "someUtil.h"
#include "evaluateCharacterFreq.h"
#include "talRandom.h"
#include "chrNumberOptions.h"
#include <cmath>
#include <algorithm>


chrNumModel::chrNumModel(const finiteIntAlphabet* pAlph, const Vdouble& freq, rootFreqType freqType, MDOUBLE halfDuplR)
: _pAlph(pAlph),_freq(freq), _rootFreqType(freqType), _halfDuplR(halfDuplR), _scaleFactor(IGNORE_PARAM)
{
	init();
}


chrNumModel::chrNumModel()
: _rootFreqType(UNIFORM), _halfDuplR(IGNORE_PARAM), _scaleFactor(IGNORE_PARAM)
{
	_freq.resize(0, -1);
}


chrNumModel::chrNumModel(const chrNumModel &other)
{
	copy(&other);
}

chrNumModel::~chrNumModel()
{
}


void chrNumModel::init()
{
    resizeMatrix(_Q, alphabetSize(), alphabetSize());
	resizeMatrix(_lastPtCalculated, alphabetSize(), alphabetSize());
	makeSureNoZeroFreqs(_freq);
}

chrNumModel& chrNumModel::operator=(const chrNumModel &other)
{
	copy(&other);
	return *this;
}

void chrNumModel::copy(const chrNumModel* pOther)
{
	_freq = pOther->_freq;
	_Q = pOther->_Q;
	_bQchanged = pOther->_bQchanged;
	_lastPtCalculated = pOther->_lastPtCalculated;
	_lastTcalculated = pOther->_lastTcalculated;
	_rootFreqType = pOther->_rootFreqType;
	_normQ = pOther->_normQ;
	_halfDuplR = pOther->_halfDuplR;
	_scaleFactor = pOther->_scaleFactor;
	_pAlph = pOther->_pAlph;
}


//set the rate parameters to random values. 
void chrNumModel::setRandomRateParams(bool bUpdateQ)
{
	if ((_halfDuplR != DEMI_EQUAL_DUPL) && (_halfDuplR != IGNORE_PARAM))
		_halfDuplR = talRandom::giveRandomNumberBetweenTwoPoints(lowerValueOfRateParam, upperValueOfRateParam);
	if (_scaleFactor != IGNORE_PARAM)
		_scaleFactor = talRandom::giveRandomNumberBetweenTwoPoints(lowerValueOfRateParam, upperValueOfScaleBranchParam);
	if (bUpdateQ == true)
		errorMsg::reportError("unexpected updateQ in chrNumModel::setRandomRateParams");
}



const MDOUBLE chrNumModel::dPij_dt(const int i,const int j, const MDOUBLE d) const {
	if (d==0.0)
        return _Q[i][j];
    errorMsg::reportError("Error in chrNumModel, dPij_dt called");
    return 0.0; // not supposed to be here
}

const MDOUBLE chrNumModel::d2Pij_dt2(const int i,const int j, const MDOUBLE d) const {
    errorMsg::reportError("Error in chrNumModel, d2Pij_dt2 called");
    return 0.0; // not supposed to be here
}


bool chrNumModel::pijt_is_prob_value(MDOUBLE val) const 
{
    if ((val+err_allow_for_pijt_function()<0) || (val>1+err_allow_for_pijt_function()))
        return false;
    else
        return true;
}

const MDOUBLE chrNumModel::Pij_t(const int i,const int j, const MDOUBLE d) const
{
	if (!_bQchanged && DEQUAL(d, _lastTcalculated))
	{
		//LOG(10, <<"PijT = "<<convert(_lastPtCalculated[i][j])<<" i = "<<i<<" j = "<<j<<" d = "<<d<<endl;);
		MDOUBLE val = convert(_lastPtCalculated[i][j]);
		if (!pijt_is_prob_value(val))
		{
			LOG(3, <<"Error in chrNumModel::Pij_t"<<endl<<" pijt= "<<val<<" i="<<i<<" j="<<j<< "t = "<<d<<endl;);
			LOGDO(3, printModelParams(myLog::LogFile())); 
			LOGDO(6, printMatrix(_Q, myLog::LogFile()));
			errorMsg::reportError("Error in chrNumModel::Pij_t, pijt <0 or >1 for t = " + double2string(d));
		}
		if (val<0.0)
			val = EPSILON; // absolute zero creates a problem later on in computations
		if (val>1.0)
			val = 1.0;
		return val;

        //return convert(_lastPtCalculated[i][j]);
	}
	zeroMatrix(_lastPtCalculated);

	exp_Mt(_Q, d, _lastPtCalculated, true);
	//VVdoubleRep Qt = multiplyMatrixByScalar(QdblRep, d);
	//VVdoubleRep unit;
	//unitMatrix(unit,_Q.size());
	//_lastPtCalculated = add(unit,Qt) ; // I + Qt
	//VVdoubleRep Qt_power = Qt;
	//VVdoubleRep prevIter_matrix = _lastPtCalculated;
	//VVdoubleRep diffM = _lastPtCalculated; //init to whatever
	//int n=2;
	//bool bConverged = false;
	//while (bConverged == false) 
	//{
	//	prevIter_matrix = _lastPtCalculated;
	//	VVdoubleRep tempQ = multiplyMatrixByScalar(Qt,1.0/n);
	//	Qt_power = multiplyMatrixes(Qt_power,tempQ);
	//	_lastPtCalculated = add(_lastPtCalculated,Qt_power); // I + Qt + Qt^2/2! + ....  + Qt^n/n!
	//	//check if the difference between the cur and prev iteration is smaller than the allowed error of all matrix entries
	//	bConverged = true;
	//	for (int row = 0; row < _lastPtCalculated.size(); ++row) {
	//		for (int col = 0; col < _lastPtCalculated.size(); ++col)
	//		{
	//			MDOUBLE pij = convert(_lastPtCalculated[row][col]);
	//			MDOUBLE diff = fabs(convert(_lastPtCalculated[row][col] - prevIter_matrix[row][col]));
	//			if ((diff > err_allow_for_pijt_function()) || (!pijt_is_prob_value(convert(_lastPtCalculated[i][j]))))
	//				bConverged = false;
	//		}
	//	}
	//	n++;
	//	if (n>250) { 
	//		string err = "Error in chrNumModel::Pij_t, too many iterations for t = " + double2string(d);
	//		errorMsg::reportError(err);
	//	}
	//}
	MDOUBLE val = convert(_lastPtCalculated[i][j]);
	if (!pijt_is_prob_value(val))
	{
		LOG(3, <<"Error in chrNumModel::Pij_t"<<endl<<" pijt= "<<val<<" i="<<i<<" j="<<j<< "t = "<<d<<endl;);
		LOGDO(3, printModelParams(myLog::LogFile())); 
		LOGDO(6, printMatrix(_Q, myLog::LogFile()));
		errorMsg::reportError("Error in chrNumModel::Pij_t, pijt <0 or >1 for t = " + double2string(d));
	}

	if (val<0.0)
		val = EPSILON; // absolute zero creates a problem later on in computations
	if (val>1.0)
		val = 1.0;
	_bQchanged = false;
	_lastTcalculated = d;
	//LOG(10, <<"PijT = "<<val<<" i = "<<i<<" j = "<<j<<" d = "<<d<<endl;);
	return val; 
}


const MDOUBLE chrNumModel::Pij_t2(const int i,const int j, const MDOUBLE d)  const {
	// converting Q into doubleRep format
	VVdoubleRep QdblRep; 
	resizeMatrix(QdblRep,_Q.size(),_Q.size());
	for (int row=0;row<_Q.size();row++){
			for (int col=0;col<_Q[row].size();col++)
				QdblRep[row][col]=convert(_Q[row][col]);
	}
	
	VVdoubleRep Qt = multiplyMatrixByScalar(QdblRep,d);
	VVdoubleRep unit;
	unitMatrix(unit,_Q.size());
	VVdoubleRep Pt = add(unit,Qt) ; // I + Qt
	VVdoubleRep Qt_power = Qt;
	MDOUBLE old_val = convert(Pt[i][j]);
	MDOUBLE diff(1.0);
	int n=2;
	while ((diff>err_allow_for_pijt_function()) || (!pijt_is_prob_value(convert(Pt[i][j])))){//(abs(old_val-new_val) > err_allow_for_pijt_function()){
		old_val = convert(Pt[i][j]);
		VVdoubleRep tempQ = multiplyMatrixByScalar(Qt,1.0/n);
		
		Qt_power = multiplyMatrixes(Qt_power,tempQ);
		Pt= add(Pt,Qt_power); // I + Qt + Qt^2/2! + ....  + Qt^n/n!
		diff = convert(Pt[i][j])-old_val; // difference is measured by diff between P[0][0] vals (a little primitive...)
		if (diff<0) diff=-diff;
		n++;
		if (n>150) { 
			string err = "Error in chrNumModel::Pij_t, too many iterations for t = " + double2string(d) + " diff = " + double2string(convert(diff));
			cerr.precision(10);
			cerr<<diff<<endl;
			LOG(3, <<err<<endl;);
			LOGDO(3, printModelParams(myLog::LogFile())); 
			//err = "";
			//for(int p = 0; p < params.size(); ++p) 
			//	err += getParamName(types[p]) + "\t" + double2string(params[p]) + "\n";
			errorMsg::reportError(err);
		}
	}
	MDOUBLE val = convert(Pt[i][j]);
	if (!pijt_is_prob_value(val))
	{
		LOG(3, <<"Error in chrNumModel::Pij_t"<<endl<<" pijt= "<<val<<" i="<<i<<" j="<<j<< "t = "<<d<<endl;);
		errorMsg::reportError("Error in chrNumModel::Pij_t, pijt <0 or >1 for t = " + double2string(d));
	}
	if (val<0.0)
		val = EPSILON; // absolute zero creates a problem later on in computations
	if (val>1.0)
		val = 1.0;

	//LOG(10, <<"PijT = "<<val<<" i = "<<i<<" j = "<<j<<" d = "<<d<<endl;);
	return val; 
}

void chrNumModel::printModelParams(ostream &out) const
{
	vector<paramObj> params;
	getParams(params);
	out<<"model params:"<<endl;
	for(int p = 0; p < params.size(); ++p) 
		out <<getParamName(params[p]._type) <<"="<<params[p]._value<<"\t";
	out<<endl;
}

void chrNumModel::printQmatrix(ostream &out) const
{
	printMatrix(_Q, out);
}


void chrNumModel::getParams(vector<paramObj>& params) const
{
	if ((_halfDuplR != DEMI_EQUAL_DUPL) && (_halfDuplR != IGNORE_PARAM))
		params.push_back(paramObj(HALF_DUPL, _halfDuplR));
	if (_scaleFactor != IGNORE_PARAM)
		params.push_back(paramObj(SCALE_BRANCH, _scaleFactor));
}

void chrNumModel::getOptimizedParams(vector<paramObj>& params) const
{
	if ((_halfDuplR != DEMI_EQUAL_DUPL) && (_halfDuplR != IGNORE_PARAM))
		params.push_back(paramObj(HALF_DUPL, _halfDuplR));
	if (_scaleFactor != IGNORE_PARAM)
		params.push_back(paramObj(SCALE_BRANCH, _scaleFactor));
}


string chrNumModel::getParamName(paramType type) const
{
	string res = "";
	switch (type)
	{
		case LOSS_CONST:
			res = "LOSS_CONST";
			break;
		case LOSS_LINEAR:
			res = "LOSS_LINEAR";
			break;
		case GAIN_CONST: 
			res = "GAIN_CONST";
			break;
		case GAIN_LINEAR: 
			res = "GAIN_LINEAR";
			break;
		case DUPL: 
			res = "DUPL";
			break;
		case HALF_DUPL: 
			res = "HALF_DUPL";
			break;
		case BASE_NUMBER_R: 
			res = "BASE_NUMBER_R";
			break;
		case BASE_NUMBER: 
			res = "BASE_NUMBER";
			break;
		case SCALE_BRANCH: 
			res = "SCALE_BRANCH";
			break;
		default: 
			errorMsg::reportError("error in chrNumModel::getParamName unknown param type");
			break;
	}
	return res;
}


MDOUBLE chrNumModel::getParamLowerBound(paramType type) const
{
	MDOUBLE res = lowerValueOfRateParam;
	return res;
}

MDOUBLE chrNumModel::getParamUpperBound(paramType type) const
{
	MDOUBLE res = 0.0;
	switch (type)
	{
	case SCALE_BRANCH:
		res = upperValueOfScaleBranchParam;
		break;
	default:
		res = upperValueOfRateParam;
	}
	return res;
}


void chrNumModel::setParam(MDOUBLE val, paramType type)
{
	switch (type)
	{
		case LOSS_CONST:
			setLossConstR(val);
			break;
		case LOSS_LINEAR:
			setLossLinearR(val);
			break;
		case GAIN_CONST: 
			setGainConstR(val);
			break;
		case GAIN_LINEAR: 
			setGainLinearR(val);
			break;
		case DUPL: 
			setDuplConstR(val);
			break;
		case HALF_DUPL: 
			setHalfDuplR(val);
			break;
		case BASE_NUMBER_R: 
			setBaseNumberR(val);
			break;
		case BASE_NUMBER: 
			setBaseNumber(static_cast<int>(floor(val)));
			break;
		case SCALE_BRANCH: 
			setScaleFactor(val);
			break;
		default: 
			errorMsg::reportError("error in chrNumModel::setParam unknown param type");
			break;
	}
}

bool chrNumModel::setLossLinearR(MDOUBLE inVal)
{
	errorMsg::reportError("cannot setLossLinearR");
	return false;
}

bool chrNumModel::setGainLinearR(MDOUBLE inVal)
{
	errorMsg::reportError("cannot setGainLinearR");
	return false;
}


bool chrNumModel::setHalfDuplR(MDOUBLE inVal)
{
	if ((_halfDuplR == DEMI_EQUAL_DUPL) || (_halfDuplR == IGNORE_PARAM))
		errorMsg::reportError("cannot setHalfDuplR - rate is fixed");
	if (inVal == _halfDuplR)
		return false;
	if ((inVal < lowerValueOfRateParam) || (inVal > upperValueOfRateParam))
        errorMsg::reportError("the halfDupl rate is not within its bounds. inVal=" + double2string(inVal));
    _halfDuplR = inVal; 
	updateQ();
	return true;
}

bool chrNumModel::setScaleFactor(MDOUBLE inVal)
{
	if (_scaleFactor== IGNORE_PARAM)
		errorMsg::reportError("cannot setScaleFactor - rate is fixed");
	if (inVal == _scaleFactor)
		return false;
	if ((inVal < lowerValueOfRateParam) || (inVal > upperValueOfScaleBranchParam))
        errorMsg::reportError("the scale factor rate is not within its bounds. inVal=" + double2string(inVal));
	_scaleFactor = inVal; 
	return true;
}


void chrNumModel::mutiplyQMatrixByScalar(MDOUBLE sc)
{
	_Q = multiplyMatrixByScalar(_Q, sc);
	LOGDO(3, printMatrix(_Q, myLog::LogFile()));
}

//calculate exp(inMat*t). The calculation is done in doubleRep to increase accuracy
//if bScaleSquare==true then perform scaling and squaring
void chrNumModel::exp_Mt(const VVdouble& inM, MDOUBLE t, VVdoubleRep& outMat, bool bScaleSquare) const
{
	int scaleF = 0;
	MDOUBLE d = t;
	if (bScaleSquare)
	{
		//scale if norm > 2^x)
		MDOUBLE norm = _normQ * t;
		//choose s such that s>(log_2{Norm} -x (and so add +1 to get the lowest int that satisfy this criterion
		scaleF = std::max<int>(0, static_cast<int>((log(norm) / log(2.0) - chrNumberOptions::_pow2Scale + 1)));
		d /= pow(2.0, scaleF);
	}
	// converting inM into doubleRep format
	VVdoubleRep MdblRep; 
	resizeMatrix(MdblRep, inM.size(), inM.size());
	for (int row=0; row<inM.size();row++){
			for (int col=0; col<inM[row].size();col++)
				MdblRep[row][col]=convert(inM[row][col]);
	}
	VVdoubleRep Mt = multiplyMatrixByScalar(MdblRep, d);
	VVdoubleRep unit;
	unitMatrix(unit, MdblRep.size());
	outMat = add(unit, Mt) ; // I + Mt
	VVdoubleRep mat_power = Mt;
	VVdoubleRep prevIter_matrix = outMat;
	VVdoubleRep diffM = outMat; //init to whatever
	int n=2;
	bool bConverged = false;
	while (bConverged == false) 
	{
		prevIter_matrix = outMat;
		VVdoubleRep tempQ = multiplyMatrixByScalar(Mt, 1.0/n);
		mat_power = multiplyMatrixes(mat_power, tempQ);
		outMat = add(outMat, mat_power); // I + Qt + Qt^2/2! + ....  + Qt^n/n!
		//check if the difference between the cur and prev iteration is smaller than the allowed error of all matrix entries
		bConverged = true;
		for (int row = 0; row < outMat.size(); ++row) {
			for (int col = 0; col < outMat.size(); ++col)
			{
				MDOUBLE diff = fabs(convert(outMat[row][col]) - convert(prevIter_matrix[row][col]));
				if ((diff > err_allow_for_pijt_function()) || (!pijt_is_prob_value(convert(outMat[row][col]))))
					bConverged = false;
			}
		}
		n++;
		if (n>250) { 
			string err = "Error in chrNumModel::Pij_t, too many iterations for t = " + double2string(t);
			LOG(3, <<err<<endl;);
			LOGDO(3, printModelParams(myLog::LogFile())); 
			LOGDO(6, printMatrix(_Q, myLog::LogFile()));
			errorMsg::reportError(err);
		}
	}

	//square matrix s times if scaling was done
	for(int i = 0; i < scaleF; ++i)
	{
		outMat = multiplyMatrixes(outMat, outMat);
	}
}


////update only those cells of the the Q matrix that are relavant to "half duplications": i->1.5i
//void chrNumModel::updateQHalfDupl(MDOUBLE halfDuplRate, MDOUBLE duplRate)
//{
//	//update chromosome=1: Q1->2 = halfDuplRate+gain+dupl
//	_Q[0][1] += halfDuplRate;
//	_Q[0][0] -= halfDuplRate;
//	//update chromosome=2: Q2->3 = halfDuplRate+gain
//	_Q[1][2] += halfDuplRate;
//	_Q[1][1] -= halfDuplRate;
//	int alphSize = alphabetSize();
//	if (alphSize < 6)
//    {//special cases
//		if (alphSize == 4)
//		{
//			_Q[2][3] += halfDuplRate;
//			_Q[2][2] -= halfDuplRate;
//		}
//		else if (alphSize == 5)
//		{
//			_Q[2][3] += halfDuplRate/2;
//			_Q[2][4] += halfDuplRate/2;
//			_Q[2][2] -= halfDuplRate;
//			_Q[3][4] += halfDuplRate;
//			_Q[3][3] -= halfDuplRate;
//		}
//		return;
//	}
//	//update chromosome=3: Q3->4 = halfDuplRate/2 + gain. Q3->5=halfDuplRate
//	_Q[2][2] -= (halfDuplRate-_Q[2][4]);
//	_Q[2][3] += halfDuplRate/2;
//	_Q[2][4] = halfDuplRate/2;
//	//update rows for even chromosomes
//	int maxChr = alphSize * 2 / 3;
//	int chr;
//	for (chr = 4; chr < maxChr; chr+=2)
//	{
//		int toChr = static_cast<int>(1.5 * chr); 
//		_Q[chr-1][chr -1] += _Q[chr-1][toChr -1];
//        _Q[chr-1][toChr -1] = halfDuplRate; 
//		_Q[chr-1][chr -1] -= halfDuplRate;
//	}
//	//update rows for odd chromosomes
//	for (chr = 5; chr < maxChr; chr+=2)
//	{
//		int toChr1 = static_cast<int>(1.5 * chr); 
//		int toChr2 = toChr1+1; 
//		_Q[chr-1][chr -1] += (_Q[chr-1][toChr1 -1]+_Q[chr-1][toChr2-1])  ;
//		_Q[chr-1][toChr1 -1] = (halfDuplRate/2);
//		_Q[chr-1][toChr2 -1] = (halfDuplRate/2);
//		_Q[chr-1][chr -1] -= halfDuplRate;
//	}
//	//update the row of maxChr
//	chr = maxChr;
//	int lastCol = alphSize-1;
//	if (alphSize%3 ==0) 
//	{//-->maxChr is always even and go to last column which already has dupl rate
//		_Q[chr-1][lastCol] += halfDuplRate;
//		_Q[chr-1][chr-1] -= halfDuplRate;
//	}
//	else if (alphSize%3 ==1) 
//	{//-->maxChr is always even and go to one column before last.
//	//This column may have been assigned in case of epsilon model
//		_Q[chr-1][chr-1] -= (halfDuplRate-_Q[chr-1][lastCol-1]);
//		_Q[chr-1][lastCol-1] = halfDuplRate;
//	}
//	else 
//	{//-->maxChr is always odd and halfDupl splits between the last 2 columns 
//		_Q[chr-1][chr-1] -= (halfDuplRate-_Q[chr-1][lastCol-1]);
//		_Q[chr-1][lastCol-1] = halfDuplRate/2; //should be empty column
//		_Q[chr-1][lastCol] += halfDuplRate/2;//should be dupl+halfDupl
//	}
//	for (chr = maxChr+1; chr < alphSize; ++chr)
//	{
//		MDOUBLE xx = _Q[chr-1][lastCol];
//		_Q[chr-1][lastCol] += halfDuplRate;
//		_Q[chr-1][chr-1] -= halfDuplRate;
//	}
//}


void chrNumModel::updateQHalfDupl(MDOUBLE halfDuplRate, MDOUBLE duplRate)
{
	int maxCount = _pAlph->max();
	if (maxCount < 6)
    {//special cases
		if (maxCount == 4)
		{
			_Q[_pAlph->count2Id(3)][_pAlph->count2Id(4)] += halfDuplRate;
			_Q[_pAlph->count2Id(3)][_pAlph->count2Id(3)] -= halfDuplRate;
		}
		else if (maxCount == 5)
		{
			_Q[_pAlph->count2Id(3)][_pAlph->count2Id(4)] += halfDuplRate/2;
			_Q[_pAlph->count2Id(3)][_pAlph->count2Id(5)] += halfDuplRate/2;
			_Q[_pAlph->count2Id(3)][_pAlph->count2Id(3)] -= halfDuplRate;
			_Q[_pAlph->count2Id(4)][_pAlph->count2Id(5)] += halfDuplRate;
			_Q[_pAlph->count2Id(4)][_pAlph->count2Id(4)] -= halfDuplRate;
		}
		return;
	}
	//update rows for counts smaller than maxCount*2/3
	int maxChrToUpdate = maxCount * 2 / 3;
	int id = 0;
	int count = _pAlph->id2Count(id);
	for (; count <= maxChrToUpdate; ++id,count = _pAlph->id2Count(id))
	{
		if ((count%2) == 0)
		{//even counts
			int toCount = static_cast<int>(1.5 * count); 
			_Q[id][_pAlph->count2Id(toCount)] += halfDuplRate; 
			_Q[id][id] -= halfDuplRate;
		}
		else
		{//odd counts
			if (count == 1)
			{
				MDOUBLE xx1 = _Q[id][id]; 
				_Q[id][_pAlph->count2Id(2)] += halfDuplRate; 
				MDOUBLE xx = _Q[id][id]; 
				_Q[id][id] -= halfDuplRate;
				MDOUBLE xx2 = _Q[id][id]; 
			}
			else 
			{
				int toChr1 = static_cast<int>(1.5 * count); 
				int toChr2 = toChr1+1; 
				_Q[id][_pAlph->count2Id(toChr1)] += (halfDuplRate/2);
				_Q[id][_pAlph->count2Id(toChr2)] += (halfDuplRate/2);
				_Q[id][id] -= halfDuplRate;
			}

		}

	}

	int lastCol = alphabetSize()-1;
	for (; count < maxCount; ++id,count = _pAlph->id2Count(id))
	{
		_Q[id][lastCol] += halfDuplRate;
		_Q[id][id] -= halfDuplRate;
	}
}

