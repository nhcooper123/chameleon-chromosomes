#include "evalParamChrNum.h"
#include "errorMsg.h"
#include "chrNumberMng.h"
using namespace std;

evalParamChrNum::evalParamChrNum(stochasticProcess* pSp, const tree& tr, const sequenceContainer &sc, paramType type, chrNumModel::rootFreqType freqType)
:_tr(tr), _sc(sc), _pSp(pSp), _paramType(type), _freqType(freqType)
{
}


MDOUBLE evalParamChrNum::operator()(MDOUBLE x)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel());
	pModel->setParam(x, _paramType);
	MDOUBLE LL = chrNumberMng::compLogLikelihood(_pSp, _tr, _sc, _freqType);
	return -LL;
}

evalParamSpeciationalProp::evalParamSpeciationalProp(stochasticProcess* pSp, const tree& originalTree, const sequenceContainer &sc, paramType type, chrNumModel::rootFreqType freqType)
:_originalTree(originalTree), _sc(sc), _pSp(pSp), _paramType(type), _freqType(freqType)
{
}

MDOUBLE evalParamSpeciationalProp::operator()(MDOUBLE x)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel());
	pModel->setScaleFactor(x);
	tree tr  = _originalTree;
	chrNumberMng::scaleTree(tr, x);
	MDOUBLE LL = chrNumberMng::compLogLikelihood(_pSp, tr, _sc, _freqType);
	return -LL;
}



