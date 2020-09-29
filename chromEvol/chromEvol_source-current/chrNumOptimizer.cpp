#include "chrNumOptimizer.h"
#include "chrNumModel.h"
#include "logFile.h"
#include "chrNumberOptions.h"
#include "numRec.h"
#include "evalParamChrNum.h"
#include "chrNumberMng.h"
#include <algorithm>
using namespace std;

chrNumOptimizer::chrNumOptimizer(const tree &tr, const sequenceContainer &sc)
:_tree(tr), _sc(sc)
{
	_diffs.clear();
	_nonDiffs.clear();
}

chrNumOptimizer::~chrNumOptimizer()
{}

MDOUBLE chrNumOptimizer::optimizeModel(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel());
	if (pModel->getScaleFactor() != IGNORE_PARAM)
		return optimizeCombinedSpeciationalModel(pSp, epsilonLLImprovement, maxIter, tolParamOptimization);


	cerr<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
	LOG(3, <<"starting optimization:"<<endl;); 
	LOGDO(5, pModel->printModelParams(myLog::LogFile())); 

	vector<paramObj> paramVec;
	_bestL = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
	LOG(5, <<"ll before optimization = "<<_bestL<<endl;);

	MDOUBLE prevIterL = VERYSMALL;
	int iter;
	for (iter = 0; iter < maxIter; ++iter) 
	{
        if (_bestL < prevIterL + epsilonLLImprovement)
			return _bestL; //likelihood converged
		//optimize exon free params
		prevIterL = _bestL;
		pModel->getOptimizedParams(paramVec); 
		LOG(3,<<"iteration: "<<iter<<" begin"<<endl;);
		for(int p = 0; p < paramVec.size(); ++p) {
			MDOUBLE paramFound;
			string paramName = pModel->getParamName(paramVec[p]._type);
			LOG(3,<<"optmizing "<<paramName<<endl;);
			MDOUBLE newL;
            MDOUBLE lowerBound = EPSILON + pModel->getParamLowerBound(paramVec[p]._type);
			MDOUBLE upperBound = pModel->getParamUpperBound(paramVec[p]._type);
			if (paramName == "BASE_NUMBER")
				newL = optimizeBaseNum(pSp, epsilonLLImprovement, maxIter, tolParamOptimization, paramFound);
			else 
				newL = -brent(lowerBound, paramVec[p]._value, upperBound, evalParamChrNum(pSp, _tree, _sc, paramVec[p]._type, chrNumberOptions::_rootFreqType),tolParamOptimization, &paramFound); 	

            if (DBIG_EQUAL(newL, _bestL) )
			{
                _bestL = newL;
                pModel->setParam(paramFound, paramVec[p]._type);
			}
            else
			{//likelihood went down!
				LOG(3,<<"inside chrNumOptimizer::optimizeModel"<<endl;);
				LOG(3,<<"like went down when optimizing: "<<paramName<<endl;);
				LOG(3,<<"param Found = "<<paramFound<<"   LL ="<<newL<<"..."<<endl;);
				LOG(3,<<"param  old = "<<paramVec[p]._value<<"   LL="<<_bestL<<"..."<<endl;);
				LOGDO(3, pModel->printModelParams(myLog::LogFile())); 
				LOGDO(3, pModel->printQmatrix(myLog::LogFile())); 

				pModel->setParam(paramVec[p]._value, paramVec[p]._type);
				errorMsg::reportError("lL went down!!!");
			}
            LOG(3,<<" LL= "<<_bestL<<" new = "<<paramFound<<" old="<<paramVec[p]._value<<endl;);
		}
	}
	if (chrNumberOptions::_rootFreqType == chrNumModel::ROOT_LL)
	{
		//the freqs should be adjusted to the current model parameters since it is not gurantee that the last LL computation is for the best params
		MDOUBLE xx = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
		if (! DEQUAL(xx, _bestL, 0.000001))
			errorMsg::reportError("error in LL computations");

	}
	return _bestL;
}

//optimize baseNumber by going over 80% of the numbers in _diffs and 20%$ of the numbers in _nonDiffs
MDOUBLE chrNumOptimizer::optimizeBaseNum(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization, MDOUBLE& paramFound)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel());
	if (_diffs.empty())
		initBaseNumberOptimization(pModel);

	random_shuffle(_diffs.begin(), _diffs.end());
	random_shuffle(_nonDiffs.begin(), _nonDiffs.end());
	int diffsLimit = ceil((_diffs.size() < 6)? _diffs.size(): 0.8 * static_cast<MDOUBLE>(_diffs.size()));
	int nonDiffsLimit = ceil(0.2 * static_cast<MDOUBLE>(_nonDiffs.size()));
	MDOUBLE maxL = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
	paramFound = pModel->getBaseNumber();
	for (int i = 0; i < diffsLimit; ++i) {
		pModel->setParam(_diffs[i], BASE_NUMBER);
		MDOUBLE ll = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
		LOG(10,<<"optimizeBaseNum, f("<<_diffs[i]<<")="<<ll<<endl);
		if (ll > maxL) {
			maxL = ll;
			paramFound = _diffs[i];
		}
	}
	for (int i = 0; i < nonDiffsLimit; ++i) {
		pModel->setParam(_nonDiffs[i], BASE_NUMBER);
		MDOUBLE ll = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
		LOG(10,<<"optimizeBaseNum, f("<<_nonDiffs[i]<<")="<<ll<<endl);
		if (ll > maxL) {
			maxL = ll;
			paramFound = _nonDiffs[i];
		}
	}

	return maxL;

}
	
void chrNumOptimizer::initBaseNumberOptimization(const chrNumModel* pModel) {
	const finiteIntAlphabet* pAlph = static_cast<const finiteIntAlphabet*>(_sc.getAlphabet());
	_diffs.clear();
	_nonDiffs.clear();
	//1. get the number of time each number appears in the seq container
	Vint counts(_sc.alphabetSize());
	int alphSize = _sc.alphabetSize();
	for (int seq = 0; seq < _sc.numberOfSeqs();++seq) 
	{
		int id = _sc.placeToId(seq);
		int charType = _sc[id][0];
		if (charType == pAlph->unknown())
			continue;
		if (pAlph->isSpecific(charType))
			++counts[charType];
		else { //composite
			Vint compIds;
			Vdouble compProbs;
			pAlph->getCompositeParts(charType, compIds, compProbs);
			for (int c = 0; c < compIds.size(); ++c)
			{
				int charType = compIds[c];
				++counts[charType];
			}
		}
	}

	//2. create a map of the of possible base transition
	map<int, int> diffs_map;
	int lowBase = static_cast<int>(pModel->getParamLowerBound(BASE_NUMBER));
	int highBase = static_cast<int>(pModel->getParamUpperBound(BASE_NUMBER));
	_diffs.clear();
	for (int i = 0; i < counts.size(); ++i) {
		if (counts[i] > 0) {
			int count_i = pAlph->id2Count(i);
			for (int j = i+1; j < counts.size(); ++j) {
				if (counts[j] > 0) {
					int count_j = pAlph->id2Count(j);
					int diff = abs(count_j - count_i);
					if ((diff < lowBase) || (diff> highBase))
						continue;
					if (diffs_map.count(diff) > 0)
						continue;
					diffs_map[diff]++;
					_diffs.push_back(diff);
				}
			}
		}
	}
	//3. make a vector of possible numbers that do not appear in the diffs vector
	for (int d = lowBase; d <= highBase; ++d) {
		if (find(_diffs.begin(), _diffs.end(), d) == _diffs.end())
			_nonDiffs.push_back(d);

	}
}

MDOUBLE chrNumOptimizer::optimizeModelManyStart(stochasticProcess* pSp, const Vint& pointsNum, const Vint& iterNum, const Vdouble& tols)
{
	//make sure that the number of points in each cycle is not bigger than the previous cycle.
	int i;
	for (i = 0; i < pointsNum.size()-1; ++i) {
		if (pointsNum[i] < pointsNum[i+1])
			errorMsg::reportError("input error in chrNumOptimizer::optimizeModelManyStart()");
	}
	//create starting models
	vector<stochasticProcess*> spVec;
	for (int i = 0; i < pointsNum[0]; ++i)
	{
		//the first model is identical to the current one
		if (i == 0)
			spVec.push_back(pSp->clone());
		else
			spVec.push_back(getRandomProcess(pSp)); 
		//LOGDO(5, static_cast<chrNumModel*>(spVec[i]->getPijAccelerator()->getReplacementModel())->printModelParams(myLog::LogFile())); 
	}
	
	int numOfOptCycles = pointsNum.size();
	Vdouble likelihoodVec;
	for (int n = 0; n < numOfOptCycles; ++n)
	{
		LOG(5, <<"=====Cycle======= " <<n<<endl;);
		if (n != 0)
		{
			//sort results and continue optimization only with the best (pointsNum[i]) points
			vector<stochasticProcess*> tmpSpVec(0); //store temporarily the best pointsNum[i] sp* 
			Vdouble sortedL = likelihoodVec;
			sort(sortedL.begin(),sortedL.end());
			MDOUBLE threshold = sortedL[sortedL.size()- pointsNum[n]];
			for (int j = 0; j < likelihoodVec.size(); ++j)
			{
				if (likelihoodVec[j] >= threshold) 
					tmpSpVec.push_back(spVec[j]);
				else
					delete spVec[j];
			}
			spVec.clear();
			spVec = tmpSpVec;
		} 

		likelihoodVec.clear();
		likelihoodVec.resize(pointsNum[n]); 
		for (int c = 0; c < pointsNum[n]; ++c)
		{
			LOG(3, <<"=====optimizing point======= " <<c<<endl;);
			stochasticProcess* ptempSp = spVec[c];
			MDOUBLE ll = optimizeModel(spVec[c], tols[n], iterNum[n], tols[n]);
			LOG(3, <<"point: "<<c<<"  likelihood = "<<ll<<endl<<endl;);
			likelihoodVec[c] = ll;
		}
	}
	//finish optimization - get best model
	Vdouble sortedL = likelihoodVec;
	sort(sortedL.begin(),sortedL.end());
	MDOUBLE _bestL = sortedL[likelihoodVec.size() - 1];
	LOG(3, <<endl<<"FINAL LIKELIHOODS++++++++++++++"<<endl;);
	for (int i = 0; i < likelihoodVec.size(); ++i)
	{
		if (_bestL == likelihoodVec[i]) 
		{
			//this operator= delete the current pointers of pSp and clone new ones
			*pSp = *spVec[i];
		}
		delete spVec[i];
		spVec[i] = NULL;
		LOG(chrNumberOptions::_logValue, <<"point "<<i<<" likelihood = "<<likelihoodVec[i]<<endl;);
		
	}	
	spVec.clear();
	return _bestL;
}



stochasticProcess* chrNumOptimizer::getRandomProcess(const stochasticProcess* pBaseSp) const
{
	stochasticProcess* pRes = pBaseSp->clone();
	static_cast<chrNumModel*>(pRes->getPijAccelerator()->getReplacementModel())->setRandomRateParams(true);
	return pRes;
}

MDOUBLE chrNumOptimizer::optimizeCombinedSpeciationalModel(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel());
	//first scale tree
	tree myTree(_tree);
	chrNumberMng::scaleTree(myTree, pModel->getScaleFactor());

	cerr<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
	LOG(3, <<"starting optimization:"<<endl;); 
	LOGDO(5, pModel->printModelParams(myLog::LogFile())); 

	vector<paramObj> paramVec;
	_bestL = chrNumberMng::compLogLikelihood(pSp, myTree, _sc, chrNumberOptions::_rootFreqType);
	LOG(5, <<"ll before optimization = "<<_bestL<<endl;);

	MDOUBLE prevIterL = VERYSMALL;
	int iter;
	for (iter = 0; iter < maxIter; ++iter) 
	{
        if (_bestL < prevIterL + epsilonLLImprovement)
			return _bestL; //likelihood converged
		//optimize exon free params
		prevIterL = _bestL;
		pModel->getParams(paramVec); 
		LOG(3,<<"iteration: "<<iter<<" begin"<<endl;);
		for(int p = 0; p < paramVec.size(); ++p) {
			MDOUBLE paramFound;
			string paramName = pModel->getParamName(paramVec[p]._type);
			LOG(3,<<"optmizing "<<paramName<<endl;);
			MDOUBLE newL;
            MDOUBLE lowerBound = pModel->getParamLowerBound(paramVec[p]._type);
			MDOUBLE upperBound = pModel->getParamUpperBound(paramVec[p]._type);
			if (paramVec[p]._type == SCALE_BRANCH)
			{
				newL = -brent(lowerBound, paramVec[p]._value, upperBound, evalParamSpeciationalProp(pSp, _tree, _sc, paramVec[p]._type, chrNumberOptions::_rootFreqType),tolParamOptimization, &paramFound); 	
				myTree = _tree;
				chrNumberMng::scaleTree(myTree, paramFound);
			}
			else 
			{
				newL = -brent(lowerBound, paramVec[p]._value, upperBound, evalParamChrNum(pSp, myTree, _sc, paramVec[p]._type, chrNumberOptions::_rootFreqType),tolParamOptimization, &paramFound); 	
			}
			if (newL >= _bestL) 
			{
				_bestL = newL;
				pModel->setParam(paramFound, paramVec[p]._type);
			}
            else
			{//likelihood went down!
				LOG(3,<<"inside chrNumOptimizer::optimizeCombinedSpeciationalModel"<<endl;);
				LOG(3,<<"like went down when optimizing: "<<paramName<<endl;);
				LOG(3,<<"param Found = "<<paramFound<<"   LL ="<<newL<<"..."<<endl;);
				LOG(3,<<"param  old = "<<paramVec[p]._value<<"   LL="<<_bestL<<"..."<<endl;);
				LOGDO(3, pModel->printModelParams(myLog::LogFile())); 
				LOGDO(3, pModel->printQmatrix(myLog::LogFile())); 
				if (paramVec[p]._type == SCALE_BRANCH) {
					myTree = _tree;
					chrNumberMng::scaleTree(myTree, paramVec[p]._value);
				}
				pModel->setParam(paramVec[p]._value, paramVec[p]._type);
				errorMsg::reportError("lL went down!!!");
			}
            LOG(3,<<" LL= "<<_bestL<<" new = "<<paramFound<<" old="<<paramVec[p]._value<<endl;);
		}
	}
	if (chrNumberOptions::_rootFreqType == chrNumModel::ROOT_LL)
	{
		//the freqs should be adjusted to the current model parameters since it is not gurantee that the last LL computation is for the best params
		MDOUBLE xx = chrNumberMng::compLogLikelihood(pSp, myTree, _sc, chrNumberOptions::_rootFreqType);
		if (xx != _bestL)
			errorMsg::reportError("error in LL computations");
	}
	return _bestL;
}

/*Optimization using gradient ascent. 
This method works only if the function is  differentiable in a neighborhood of a point, and thus cannot be applied to integer parameters such as the base number.
At each point we want to find the gradients of the likelihood function (at the current point) relative to each optmizaed parameter and then take steps at the direction of the positive of the gradient.
The step size is chosen based on a line search: first we try the farthest point (the one that will bring us to the boundary of the parameter space), and then we iteratively half the step size.
*/
//MDOUBLE chrNumOptimizer::optimizeModel_GA(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization)
//{
//	chrNumModel* pModel = static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel());
//
//	cerr<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
//	LOG(3, <<"starting optimization gradient ascent:"<<endl;); 
//	LOGDO(5, pModel->printModelParams(myLog::LogFile())); 
//
//	_bestL = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
//	LOG(5, <<"ll before optimization = "<<_bestL<<endl;);
//
//	stochasticProcess leftSp = *pSp;
//	stochasticProcess rightSp = *pSp;
//
//	MDOUBLE prevIterL = VERYSMALL;
//	vector<MDOUBLE> gradient;
//	int iter;
//	bool hit_the_edge = false;
//	for (iter = 0; iter < maxIter; ++iter) 
//	{
//        if ((_bestL < prevIterL + epsilonLLImprovement) && !hit_the_edge) //if we hit the boundary then we still want to do one more iteration since we can still find a better maxima along the ridge
//			break; //likelihood converged
//		hit_the_edge = true; // so the next iteration will have to change it back 
//		prevIterL = _bestL;
//
//		LOG(3,<<"iteration: "<<iter<<" begin"<<endl;);
//		have to go over the maximization
//		vector<paramObj> paramLeft, paramRight;
//		static_cast<chrNumModel*>(leftSp.getPijAccelerator()->getReplacementModel())->getOptimizedParams(paramLeft);  //we should not get baseNUmber!!
//		find_gradient(paramLeft, _bestL, &leftSp, tolParamOptimization, gradient);
//		chrNumModel* pRightModel = static_cast<chrNumModel*>(rightSp.getPijAccelerator()->getReplacementModel());
//		pRightModel->getOptimizedParams(paramRight);
//		for(int p = 0; p < paramRight.size(); ++p) {
//			MDOUBLE newVal = paramLeft[p]._value + gradient[p];
//			pRightModel->setParam(newVal, paramRight[p]._type); //should make sure newVal is non-negative
//		}
//		MDOUBLE leftPointLikelihood = _bestL, rightPointLikelihood = VERYSMALL;
//		string valid = "left";
//
//		while (getDistance(paramRight, paramLeft) > tolParamOptimization) {
//			_bestL = gradient_walk_step(paramLeft, paramRight, &leftSp, leftPointLikelihood, rightPointLikelihood, valid);  
//			if (valid=="left") // we don't go with the gradient all the way.
//				hit_the_edge = false;
//		}
//
//        //print new params to log
//		//LOG(3,<<" LL= "<<_bestL<<" new = "<<paramFound<<" old="<<paramVec[p]._value<<endl;);
//	}
//	
//	if (chrNumberOptions::_rootFreqType == chrNumModel::ROOT_LL)
//	{
//		//the freqs should be adjusted to the current model parameters since it is not gurantee that the last LL computation is for the best params
//		MDOUBLE xx = chrNumberMng::compLogLikelihood(pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
//		if (! DEQUAL(xx, _bestL, 0.000001))
//			errorMsg::reportError("error in LL computations");
//
//	}
//	return _bestL;
//}

//practically - makes one point of line maximization in order to find the best step size.
//calcualtes the middle point between left and right. Then set the middle point to the left/right according to the point with lower likelihood 
//Valid: tells whether the left or right was caluclated before (so if middle was set to the right then the "new right" should be calcualted again while left is valid)
//returns the new likelihood 
MDOUBLE chrNumOptimizer::gradient_walk_step(vector<paramObj>& paramLeft, vector<paramObj>& paramRight, const stochasticProcess* pSp, MDOUBLE &left_likelihood, MDOUBLE &right_likelihood, string &valid) {
	//have to make sure that the pSp that is used here does not change the pSp that was called
	MDOUBLE res = VERYSMALL;
	stochasticProcess sp = *pSp;
	chrNumModel* pModel = static_cast<chrNumModel*>(sp.getPijAccelerator()->getReplacementModel());
	vector<paramObj> middlePoint = paramLeft;
	if(valid != "left")
	{
		for (int p = 0; p < paramLeft.size(); ++p) {
			pModel->setParam(paramLeft[p]._value, paramLeft[p]._type);
		}
		left_likelihood = chrNumberMng::compLogLikelihood(&sp, _tree, _sc, chrNumberOptions::_rootFreqType);
	}
	if(valid != "right")
	{
		for (int p = 0; p < paramRight.size(); ++p) {
			pModel->setParam(paramRight[p]._value, paramRight[p]._type);
		}
		left_likelihood = chrNumberMng::compLogLikelihood(&sp, _tree, _sc, chrNumberOptions::_rootFreqType);
	}
	for(int p = 0; p < paramRight.size(); ++p)
	{
		middlePoint[p]._value = ((paramRight[p]._value + paramLeft[p]._value)/2);	
	}

	if(left_likelihood > right_likelihood)
	{
		res = left_likelihood;
		paramRight = middlePoint;
		valid = "left";
	}
	else
	{
		res = right_likelihood;
		paramLeft = middlePoint;
		valid = "right";
	}
	return res;
}

//find_gradient: the gradient is not only the direction of the vector but also its length (size). 
//Here, we want to find the gradient vector that will hit the nearest boundary of the parameter space. 
void chrNumOptimizer::find_gradient(vector<paramObj>& paramVec, MDOUBLE curL, const stochasticProcess* pSp, MDOUBLE epsilon, vector<MDOUBLE>& gradient) {
	gradient.clear();
	stochasticProcess sp = *pSp;
	chrNumModel* pModel = static_cast<chrNumModel*>(sp.getPijAccelerator()->getReplacementModel());
	MDOUBLE mul = 1000/epsilon; //the multiplier of the initial gradient found that will bring the point to the nearest edge
	MDOUBLE temp_mul = 1000/epsilon;

	for (int p = 0; p < paramVec.size(); ++p) {
		MDOUBLE lowerBound = 0.0001 + pModel->getParamLowerBound(paramVec[p]._type);
		MDOUBLE upperBound = -0.0001 + pModel->getParamUpperBound(paramVec[p]._type);

		MDOUBLE origVal = paramVec[p]._value;
		MDOUBLE newVal = origVal + epsilon;
		pModel->setParam(newVal, paramVec[p]._type);
		MDOUBLE tempL = chrNumberMng::compLogLikelihood(&sp, _tree, _sc, chrNumberOptions::_rootFreqType);
		//if tempL>curL but we are at the boundary of the parameter space - then we don't want to proceed more at that direction.
		//Similarly, if tempL<curL then the next steps will be at the opposite direction of current+epsilon, but we cannot make any additional steps at that direction if we are already very close to the boundary. 
		if ((tempL>curL && newVal>upperBound) || (tempL<curL && (origVal-epsilon)<lowerBound)) {//so we won't go off the edge
			gradient.push_back(0);
		}
		else {
			//here we want to find the shortest addition of the gradient from the current point that will bring us to the edge of the parameter space. 
			//size is the multiplier we want to find (we will multiply the gradient by this number).
			//for each parameter we find the multiplier that will hit the edge in its direction. The final multiplier will be the smallest
			MDOUBLE grad = (tempL-curL)/epsilon;
			gradient.push_back(grad); 
			MDOUBLE step_wanted;
			if (grad > 0)
				step_wanted = upperBound-origVal;
			else if (grad < 0)
				step_wanted = lowerBound-origVal;
			temp_mul = step_wanted/gradient[p]; //divide the step wanted by the current value of the gradient 
			if(temp_mul<mul && temp_mul>0)
				mul = temp_mul;
		}
		pModel->setParam(origVal, paramVec[p]._type);
	}
	if (paramVec.size() != gradient.size())
		errorMsg::reportError("error in find_gradient: the paramVec must be the same size as gradientVec");
	//resize gradient
	for (int p = 0; p < paramVec.size(); ++p) {
		gradient[p] = gradient[p]*mul;
	}
}

MDOUBLE chrNumOptimizer::getDistance(const vector<paramObj>& paramLeft, const vector<paramObj>& paramRight) {
	int dimension = paramLeft.size();
	MDOUBLE dist=0;
	for(int i=0; i < dimension; ++i){
		dist += pow(paramLeft[i]._value - paramRight[i]._value,2);
	}
	dist=sqrt(dist);
	return dist;
}
