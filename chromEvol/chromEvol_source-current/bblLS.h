#ifndef ___BBL_LS
#define ___BBL_LS

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "exonIntronModel.h"
#include "codonPositionModel.h"
#include "codonPositionModelLight.h"
using namespace std;

#define MAX_BRANCH_LENGTH 10.0

/*
This class optimize the branches using "naive" line search methodology.
go over each branch and optimize it using brent.
In one iteration it optimze seperatly all branches.
This procedure continues until convergence is reached or until the maximum number of iteration is reached.
*/
class bblLS {
public:
	
	explicit bblLS();
	~bblLS() {};
	MDOUBLE optimizeBranches(tree& et, const sequenceContainer &intronSc, const sequenceContainer &exonSc,
						exonIntronModel& sp, int maxIter=50, MDOUBLE epsilon=0.05, MDOUBLE curL = NULL);
	MDOUBLE optimizeBranches(tree& et, const sequenceContainer &sc,
						codonPositionModel& sp, int maxIter=50, MDOUBLE epsilon=0.05, MDOUBLE curL = NULL);
	MDOUBLE optimizeBranches(tree& et, const sequenceContainer &sc,
						codonPositionModelLight& sp, int maxIter=50, MDOUBLE epsilon=0.05, MDOUBLE curL = NULL);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}




	//alon:  these function are the ones I maed to implement the branch optimization with gradient descent
	MDOUBLE alon_optimizeBranches(tree& et, const sequenceContainer &sc,
						codonPositionModelLight& sp, int maxIter=50, MDOUBLE epsilon=0.05, MDOUBLE curL = NULL);

	vector<MDOUBLE> alon_aid_find_gradient(tree& et, const sequenceContainer &sc,
						codonPositionModelLight& sp, vector<tree::nodeP> nodesV, MDOUBLE epsilon, MDOUBLE &curL, MDOUBLE min_value = 0, MDOUBLE max_value = 2);

	void alon_aid_gradient_walk_step(tree& leftPoint, tree& rightPoint, const sequenceContainer &sc,codonPositionModelLight& model, MDOUBLE &left_likelihood, MDOUBLE &right_likelihood, string &valid);

	MDOUBLE alon_aid_distance_of_two_points(tree leftPoint, tree rightPoint);
	//end of gradient descent functions




private:
	MDOUBLE _treeLikelihood;
};

class evalBranch{
public:
	explicit evalBranch(tree::nodeP pNode, tree* pTree, const sequenceContainer &scIntron,const sequenceContainer &scExon, exonIntronModel &model)
		:_pNode(pNode),_pTree(pTree), _scIntron(scIntron),_scExon(scExon), _model(model){};
	MDOUBLE operator() (MDOUBLE x);

private:
	tree::nodeP _pNode;
	tree* _pTree;
	const sequenceContainer& _scIntron;
	const sequenceContainer& _scExon;
	exonIntronModel & _model;
};

class evalCodonPositionBranch{
public:
	explicit evalCodonPositionBranch(tree::nodeP pNode, tree* pTree, const sequenceContainer &sc, codonPositionModel &model)
		:_pNode(pNode),_pTree(pTree), _sc(sc), _model(model){};
	MDOUBLE operator() (MDOUBLE x);

private:
	tree::nodeP _pNode;
	tree* _pTree;
	const sequenceContainer& _sc;
	codonPositionModel & _model;
};

class evalCodonPositionBranchLight{
public:
	explicit evalCodonPositionBranchLight(tree::nodeP pNode, tree* pTree, const sequenceContainer &sc, codonPositionModelLight &model)
		:_pNode(pNode),_pTree(pTree), _sc(sc), _model(model){};
	MDOUBLE operator() (MDOUBLE x);

private:
	tree::nodeP _pNode;
	tree* _pTree;
	const sequenceContainer& _sc;
	codonPositionModelLight & _model;
};

#endif
