#ifndef ___SIMULATE_JUMPS__LARGE
#define ___SIMULATE_JUMPS__LARGE
#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h" 

#include <map>
#include <vector>
using namespace std;

/****************************************************************** 
This class simulates jumps (events) along differing branch lengths (according to a 
given tree), with the aim of giving the expectation of the number of jumps
from state A to state B given that the terminal states at the end of the branch are
x and y.
*******************************************************************/

class simulateJumpsLargeAlphabet  {
public:
	simulateJumpsLargeAlphabet(const tree& inTree, const stochasticProcess& sp, alphabet* pAlph);
	virtual ~simulateJumpsLargeAlphabet();
	//runSimulation: do the actual simulation. iterNum specifies the number of iterations starting from each state
	void runSimulation(int iterNum); 
	
	//for a branch length specified by a nodeName: 
	//give the expected number of jumps (changes) from fromId to toId that occured along the specified branh length, 
	//in which the starting character is terminalStart and the terminal character is terminalEnd
	MDOUBLE getExpectation(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId);
	//same as above, except here we return the probability of a jump from fromId to toId given 
	//terminal states terminalStart, terminalEnd in this branch
	MDOUBLE getProb(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId);
    	
private:
	void init();
	void runOneIter(int state);
	int getCombinedState(int terminalStart, int terminalEnd) const;
	int getCombinedAlphabetSize() const {return _pAlph->size()*_pAlph->size();}
	int getStartId(int combinedState) const;
	int getEndId(int combinedState) const;
	void computeExpectationsAndPosterior();
	

private:
	tree _tree;
	stochasticProcess _sp;
	alphabet* _pAlph;

	Vdouble _waitingTimeParams;//each entry is the lambda parameter of the exponential distribution modeling the waiting time for "getting out" of state i
	//_jumpProbs[i][j] is the probability of jumping from state i to state j (given that a change has ocured).
	VVdouble _jumpProbs; 
	//_node2Jumps: maps a node name (which specify a branch length) to 
	//the expected number of jumps between any two characters along the branch leading from the father to this node
	//given the terminal characters of this branch.
	//The matrix is 2D and not 4D because we use a "combined alphabet" to make access easier. see getCombinedState() for details
	//The first dimension is the combined terminal state and the second dimension is the combined jump state
	//map<string, VVdouble> _nodes2JumpsExp; 
	typedef map<int, map<int, MDOUBLE> > map_i2i2double ;
	map<string, map_i2i2double > _nodes2JumpsExp;
	
	//the number of times we reached a certain combination of terminal states for each branch lengths
	//e.g. the number of times we observed 0,1 at terminal states given branch length 0.03
	//this is used to to afterwards normalize (i.e. compute the expectation) the _nodes2JumpsExp values
	map<string, Vdouble> _totalTerminals; 

	//_node2JumpsProb: maps a node name (which specify a branch length) to 
	//the probability of a jump between any two characters along the branch leading from the father to this node
	//given the terminal characters of this branch.
	//The matrix is 2D and not 4D because we use a "combined alphabet" to make access easier. see getCombinedState() for details
	//The first dimension is the combined terminal state and the second dimension is the combined jump state
	//map<string, VVdouble> _nodes2JumpsProb; 


	vector<tree::nodeP> _orderNodesVec; //internal use: the branch are sorted in ascending order 

};

#endif
