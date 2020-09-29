#ifndef ___SIMULATE_CHANGES__
#define ___SIMULATE_CHANGES__

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"
#include "sequenceContainer.h"
#include "suffStatComponent.h"
#include "chrNumModel.h"

#include <map>
#include <vector>
using namespace std;

/******************************************************************
This class simulates jumps (events) along  a 
given tree, with the aim of creating a dataset (seqeunceContainer) 
in which we know the exact number of transitions along the tree.
Can also use stochastic mappings to infer the expected number (and type) of transitions that have occured 
*******************************************************************/

class simulateChangesAlongTree  {
public:
	simulateChangesAlongTree(const tree& inTree, const stochasticProcess& sp, finiteIntAlphabet* pAlph);
	virtual ~simulateChangesAlongTree();
	sequenceContainer simulatePosition(); 
	//simulate mutations on the phylogeny a given number of times so the resulting leaves are identical to the input data, 
	//as given in the specified position of the sequenceContainer. 
	//The result is stored in the as the history of mutations for each branch.
	//retruns false if fails to generate mapping for the specified data
	VVdouble generateStochasticMapping(const sequenceContainer& inSc, int iterations, const VVVdouble& jointPost, const suffStatGlobalHomPos& sscUp);

	VVint getChangesForBranch(int nodeID);
	int getNodeContent(int nodeId) {return _nodesContent[nodeId];}
	sequenceContainer toSeqDataWithoutInternalNodes();
	int getTotalJumps(chrNumModel::jumpType type);//0=gain/1=loss/2=dupl/3=demi
	void printEvents(ostream &expFile);

private:
	void init();
	void simulateOnce(tree::nodeP curNode, MDOUBLE disFromNode, int previousContent, int whichSon = 0);
	void calcEventsTypes(); //fill _changesTypes: the number of gain/loss/dupl/demi for each node
	void computeTransitionsFromRootToNodes(VVdouble& nodesTransitions);
	//the following functions are for stochastic mapping
	VVVdouble getSamplingProbs(const suffStatGlobalHomPos& sscUp) const;
	void sampleAncestrals(const Vdouble& rootProbs, const VVVdouble& samplingProbs);
	void sampleMutationsGivenAncestrals();
	vector<chrNumModel::jumpType> sampleMutationsGivenAncestralsPerBranch(int fatherState, int sonState, MDOUBLE branchLength);
	chrNumModel::jumpType getJumpType(int curState, int nextState);
	

private:
	tree _tree;
	stochasticProcess _sp;
	finiteIntAlphabet* _pAlph;

	Vdouble _waitingTimeParams;//each entry is the lambda parameter of the exponential distribution modeling the waiting time for "getting out" of state i
	//_jumpProbs[i][j] is the probability of jumping from state i to state j (given that a change has ocured).
	VVdouble _jumpProbs; 

	VVVint _changesOccurred; // number of times changes from i to j occurred , for each branch
	VVint _changesTypes; //_changesTypes[node][jumpType]; for each node, the number of times for each jump type: gain, loss, dupl, demiDupl, base, or to maxChr.
						//lastPlace is total across tree
	Vint _nodesContent; // the actual state at each node, retrieval according to node id
	sequenceContainer _sc;
};

#endif
