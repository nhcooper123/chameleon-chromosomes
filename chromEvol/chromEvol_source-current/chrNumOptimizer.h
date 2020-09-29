#ifndef __CHROMOSOME_NUMBER_OPTIMIZER
#define __CHROMOSOME_NUMBER_OPTIMIZER

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "chrNumModel.h"



/***********************************************************************
**********************************************************************/
class chrNumOptimizer{
public:
	chrNumOptimizer(const tree &tr, const sequenceContainer &sc);
	~chrNumOptimizer();
	//optimize the model param and return best log-likelihood
	MDOUBLE optimizeModel(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization);
	MDOUBLE optimizeBaseNum(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization, MDOUBLE& paramFound);

	MDOUBLE optimizeModel_GA(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization);
    MDOUBLE optimizeModelManyStart(stochasticProcess* pSp, const Vint& pointsNum, const Vint& iterNum, const Vdouble& tols);

private:
	stochasticProcess* getRandomProcess(const stochasticProcess* pBaseSp) const;
	MDOUBLE optimizeCombinedSpeciationalModel(stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, int maxIter, MDOUBLE tolParamOptimization);

	void find_gradient(vector<paramObj>& paramVec, MDOUBLE curL, const stochasticProcess* pSp, MDOUBLE epsilonLLImprovement, vector<MDOUBLE>& gradient);
	MDOUBLE gradient_walk_step(vector<paramObj>& paramLeft, vector<paramObj>& paramRight, const stochasticProcess* pSp, MDOUBLE &left_likelihood, MDOUBLE &right_likelihood, string &valid);
	MDOUBLE getDistance(const vector<paramObj>& paramLeft, const vector<paramObj>& paramRight); 
	void initBaseNumberOptimization(const chrNumModel* pModel);

private:
	tree _tree;
	sequenceContainer _sc;

	//for base number optimizations
	Vint _nonDiffs;
	Vint _diffs;

	MDOUBLE _bestL;
};
#endif

