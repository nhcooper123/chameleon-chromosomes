#ifndef __CHR_NUMBER_MNG___H
#define __CHR_NUMBER_MNG___H

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "finiteIntAlphabet.h"
#include "chrNumModel.h"
class computePijHom;
class simulateJumpsLargeAlphabet;
class suffStatGlobalHomPos;

class chrNumberMng {

public:
	explicit chrNumberMng();
	virtual ~chrNumberMng();
    void run();
	void runSimulations();
	void runSimulationsAndInference();
	void runOneSimulation();
	void loopRemoveTaxa();
	void validateData();
	MDOUBLE getBestLL() const {return _bestLL;}
	static MDOUBLE compLogLikelihood(const stochasticProcess* pSp, const tree &tr, const sequenceContainer &sc, chrNumModel::rootFreqType freqType);
	static void scaleTree(tree &tr, MDOUBLE scaleFactor);
	const chrNumModel* getModel() {return static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel());}
	static Vdouble estimateEvents(int fromId, int toId, const stochasticProcess* pSp, const finiteIntAlphabet* pAlph); // a stupid heuristic to estimate the type of events in case it couldn't be computed from the simulations
	Vdouble estimateExpectationsHeuristicly(int nodeID, const VVVdouble& jointPost); //same as estimateEvents but go over all probabilites for the two terminal nodes

private:
	void init(); 
	void clear();
	void printRunInfo();
	void initStartingModel(const finiteIntAlphabet* pAlph, int allowedRange);
	void readFreqsFromFile(const string& inFile, Vdouble& freqs, const finiteIntAlphabet* pAlph);
	pair<int, int> getChrRangeInData(const sequenceContainer& sc); //return <min,max> counts in data
	void getTree();
	void rootTree();
	void printTreeStatesAsBPValues(ostream &out, Vint &states, const tree::nodeP &myNode, VVVdouble *probs, bool printGains, const alphabet* pAlph);
	void initializeStatesVector(Vint& statesOut, const tree& inTree, const sequenceContainer& inSc);
	void ancestralReconstructML(Vint& statesOut);
    void computeAncestralPosterior(const VVVdouble& jointPost); //recieves the precalculated joint post. probs and computes P(N=x|Data) and stores in _ancestralPost[nodeId][letter]
	void computeJointPosterior(VVVdouble &jointPost, suffStatGlobalHomPos& sscUp);//computes P(N=x, father(N)=y|Data), also - return a filled suffStatGlobalHomPos used in later calculations
	MDOUBLE getHeuristicProbOfChange(int nodeID, const VVVdouble& jointPost);////recieves the precalculated joint post and return an estimate of the probability that a change has occured along the branch leading to nodeID.
	//compute LL given the pre-calculated Pij
	static MDOUBLE compLogLikelihood(const computePijHom& pij, const tree &tr, const sequenceContainer &sc, chrNumModel::rootFreqType freqType, const stochasticProcess* pSp);
	static doubleRep compLikelihood(const computePijHom& pij, const tree &tr, const sequenceContainer &sc, chrNumModel::rootFreqType freqType, const stochasticProcess* pSp);
	void traverseUpML(VVdouble &upL, VVint &backtrack, const Vint& states); //states contains the state at the leaves
	MDOUBLE traverseDownML(const VVdouble &upL, const VVint &backtrack, Vint& reconstructStates, VVint& transitionTypeCount);
	//computes the number of transition in which there was a change along the tree of +1, +2,..,+(maxChrSize-1). Similarily for the down changes
	void getTransitionVec(const tree& tr, const Vint& states, Vint& upChangeVec, Vint& downChangeVec);

	MDOUBLE optimizeParams(bool bScaleTree = false); //returns the best log-likelihood
	void inferAncestralStates(const VVVdouble& jointPost);//infer ML and posterior states. receives the joint posterior p(N=x, Father(N)=y|D)
	void computeChangeProbAndExp(const VVVdouble& jointPost, const suffStatGlobalHomPos& sscUp);//infer the expected number of changes between two states AND the probability for each type of change along the tree. receives the joint posterior p(N=x, Father(N)=y|D)
	void computeChangeProbAndExp_SM(const VVVdouble& jointPost, const suffStatGlobalHomPos& sscUp);
	MDOUBLE computeExpectationOfChange(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost, int fromState, int toState);
	MDOUBLE computeExpectationOfChangePerBranch(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost, tree::nodeP node, int fromState, int toState);
	void computeChangeProbAndExpBaseNumber(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost);
	//for each leaf calculate the number of transitions (for each event type) from root to leaf
	void computeExpectationsFromRootToNodes(VVdouble& nodesTransitions);
	MDOUBLE computeProbOfChangePerBranch(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost, tree::nodeP node, int fromState, int toState);
	
	//simulation related
	void getSimulationResults(int iter, Vint& simStates, VVdouble& params, Vint& mlRoot, Vdouble& postRoot, Vdouble& mlAllInternals, Vdouble& postAllInternals);
	void getSimulatedTree(int iter);

	void getStartingData();
	void changesNamesToLeaves(const tree& inTree, const sequenceContainer& inSc); //append the state of the leaves to their names
	// TREE SEARCH PART
	void getStartingTreeFromTreeFile();
	void getStartingTreeNJ_fromDistances(const VVdouble& disTab,const vector<string>& vNames);
	
	void printResults(const VVVdouble& jointPost);
	void printReconstructedTrees(const string& mlTree, const string& postTree); //print the ML and posterior reconstruction
	void printPosterior(ostream &out, const VVdouble& postProbs);
	void printExpectations(const VVVdouble& jointPost, ostream &out);

	bool performSimulationsForExpectation(); //return true if should perform simulations that should be used to compute the expectation of the number of events
	void getBaseTransitions(const string& baseTransitionProbs, vector<pair<int, MDOUBLE> > & compProbs);

private:
	stochasticProcess* _pSp;
	sequenceContainer _sc; 
	tree _tree;
	finiteIntAlphabet* _pAlph;
	MDOUBLE _bestLL;
	VVdouble _ancestralPost;//the posterior at ancestral nodes _ancestralPost[nodeId][letter]
	Vint _ancestralMl; //the ML reconstruction of ancestral nodes _ancestralMl[nodeId]
	VVdouble _expChanges; // for each node, the expected number of changes for [0-4]: gain, loss, dupl, halfDupl, baseNumber
	//VVdouble _probChanges; // for each node, the probability the event ([0-3]: gain, loss, dupl, halfDupl) occured along the branch from nodeId to its father
};
#endif
