#include "simulateJumpsLargeAlphabet.h"
#include "talRandom.h"
#include "someUtil.h"
#include <algorithm>


simulateJumpsLargeAlphabet::simulateJumpsLargeAlphabet(const tree& inTree, const stochasticProcess& sp, alphabet* pAlph)
: _tree(inTree), _sp(sp), _pAlph(pAlph)	
{
}

simulateJumpsLargeAlphabet::~simulateJumpsLargeAlphabet()
{
}


// a comparison function to be used in sort init
bool compareDist(tree::nodeP node1, tree::nodeP node2)
{
	return (node1->dis2father() < node2->dis2father());
}

void simulateJumpsLargeAlphabet::init()
{
	//init the vector of waiting times. 
	_waitingTimeParams.clear();
	_waitingTimeParams.resize(_pAlph->size());
	int i, j;
	for (i = 0; i < _pAlph->size(); ++i)
	{
		_waitingTimeParams[i] = -_sp.dPij_dt(i, i, 0.0);
		if (_waitingTimeParams[i] < 0.0)
			errorMsg::reportError("problem in simulateJumpsLargeAlphabet::init for state " + int2string(i));
		
	}
	//init _jumpProbs.
	//_jumpProbs[i][j] = Q[i][j] / -Q[i][i]
	_jumpProbs.clear();
	_jumpProbs.resize(_pAlph->size());
	for (i = 0; i < _pAlph->size(); ++i)
	{
		MDOUBLE sum = 0.0;
		_jumpProbs[i].resize(_pAlph->size());
		for (j = 0; j < _pAlph->size(); ++j)
		{
			if (i == j)
				_jumpProbs[i][j] = 0.0;
			else
			{
				_jumpProbs[i][j] = _sp.dPij_dt(i, j, 0.0) / _waitingTimeParams[i];
			}
			sum += _jumpProbs[i][j];
		}
		if (! DEQUAL(sum, 1.0)){
			string err = "error in simulateJumpsLargeAlphabet::init(): sum probabilities is not 1 and equal to ";
			err+=double2string(sum);
			errorMsg::reportError(err);
		}
	}

	//init _orderNodesVec: a vector in which the branch lengths are ordered in ascending order
	_tree.getAllNodes(_orderNodesVec, _tree.getRoot());
	sort(_orderNodesVec.begin(), _orderNodesVec.end(), compareDist); 

	_nodes2JumpsExp.clear();
	//_nodes2JumpsProb.clear();
//	VVdouble zeroMatrix(getCombinedAlphabetSize());
//	for (i = 0; i < getCombinedAlphabetSize(); ++i)
//		zeroMatrix[i].resize(getCombinedAlphabetSize(), 0.0);
	Vdouble zeroVector(getCombinedAlphabetSize(),0.0);
	for (i = 0; i < _orderNodesVec.size(); ++i)
	{
		string nodeName = _orderNodesVec[i]->name();
//		_nodes2JumpsExp[nodeName] = zeroMatrix;
		//_nodes2JumpsProb[nodeName] = zeroMatrix;
		_nodes2JumpsExp[nodeName][0][0] = 0;
		_totalTerminals[nodeName]=zeroVector;
//		for (j=0; j<getCombinedAlphabetSize();++j)
//			_totalTerminals[nodeName]=zeroVector;
	}
}


//runSimulation: do the actual simulation. iterNum specifies the number of iterations starting from each state
void simulateJumpsLargeAlphabet::runSimulation(int iterNum)
{
	init();
	for (int state = 0; state < _pAlph->size(); ++state)
	{
		LOG(5, <<"simulaing state "<<state);
		for (int iter = 0; iter < iterNum; ++iter)
		{
			LOGnOUT(7,<<"state = "<<state<<" iter = "<<iter<<endl);
			runOneIter(state);
		}
	}
	
	computeExpectationsAndPosterior();	
}


//simulate jumps starting from startState. The simulation continue until the maxTime is reached. In each step:
//1. Draw a new waiting time.
//2. Go over all branches shorter than nextJumpTime and update their jumpsNum between the states that were switched 
//	(these branches will not be affected by the current jump): 
//	however they might have been affected by the previous jump
//3. Draw a new state
void simulateJumpsLargeAlphabet::runOneIter(int startState)
{
	MDOUBLE maxTime = _orderNodesVec[_orderNodesVec.size()-1]->dis2father();
	MDOUBLE totalTimeTillJump = 0.0;
	int jumpsNum = 0;
	int curState = startState;
	int smallestBranchNotUpdatedSofar = 0;
	vector<pair<int, int> > jumpsSoFar(0);
	while (totalTimeTillJump < maxTime)
	{
		MDOUBLE avgWaitingTime = 1 / _waitingTimeParams[curState];
		MDOUBLE nextJumpTime = totalTimeTillJump + talRandom::rand_exp(avgWaitingTime);
		//go over all branches that "finished" their simulation (shorter than nextJumpTime) and update with their _nodes2JumpsExp 
		//with the jumps that occured between the terminal Ids: startState-->curState
		for (int b = smallestBranchNotUpdatedSofar; b < _orderNodesVec.size(); ++b)
		{
			if (_orderNodesVec[b]->dis2father() > nextJumpTime)
			{
				smallestBranchNotUpdatedSofar = b;
				break;
			}
			string nodeName = _orderNodesVec[b]->name();
			//update all the jumps that occured along the branch
			int terminalState = getCombinedState(startState, curState);
			_totalTerminals[nodeName][terminalState]++;
			//update all longer branches with all jumps that occurred till now
			vector<bool> jumpsSoFarBool(getCombinedAlphabetSize(),false);
			for (int j = 0; j < jumpsSoFar.size(); ++j)
			{
				int combinedJumpState = getCombinedState(jumpsSoFar[j].first, jumpsSoFar[j].second);
				jumpsSoFarBool[combinedJumpState]=true;
                _nodes2JumpsExp[nodeName][terminalState][combinedJumpState] += 1;
			}
			//for (int combined=0;combined<jumpsSoFarBool.size();++combined)
			//{
			//	if (jumpsSoFarBool[combined])
			//		_nodes2JumpsProb[nodeName][terminalState][combined]+=1;
			//}
		}
		totalTimeTillJump = nextJumpTime;
		int nextState = giveRandomState(_pAlph->size(),curState, _jumpProbs);
		jumpsSoFar.push_back(pair<int,int>(curState, nextState));
		curState = nextState;
		++jumpsNum;
	}
}



void simulateJumpsLargeAlphabet::computeExpectationsAndPosterior(){
	//scale _nodes2JumpsExp so it will represent expectations
	map<string, map_i2i2double >::iterator iterExp = _nodes2JumpsExp.begin();
	for (; iterExp != _nodes2JumpsExp.end(); ++iterExp)
	{
		string nodeName = iterExp->first;
		map_i2i2double::iterator termIter = iterExp->second.begin();
		for (; termIter != iterExp->second.end(); ++termIter)
		//for (int termState = 0; termState < getCombinedAlphabetSize(); ++termState)
		{
			int termState = termIter->first;
			map<int, MDOUBLE>::iterator jumpIter = termIter->second.begin();
			for (; jumpIter != termIter->second.end(); ++jumpIter)
			//for (int jumpState = 0; jumpState < getCombinedAlphabetSize(); ++jumpState)
			{
				int jumpState = jumpIter->first;
				map<string, Vdouble>::iterator iterTerm = _totalTerminals.find(nodeName);
				//map<string, VVdouble>::iterator iterProb = _nodes2JumpsProb.find(nodeName);
				if (iterTerm==_totalTerminals.end())// || (iterProb==_nodes2JumpsProb.end()))
				{
					errorMsg::reportError("error in simulateJumpsLargeAlphabet::runSimulation, unknown reason: cannot find nodeName in map");
				}
				if ((iterTerm->second[termState]==0)){ //never reached these terminal states
					if (iterExp->second[termState][jumpState]==0)// && (iterProb->second[termState][jumpState]==0))
						continue;//leave the value of _nodes2JumpsExp and _nodes2JumpsProb as zero
					else {
						errorMsg::reportError("error in simulateJumpsLargeAlphabet::runSimulation, 0 times reached termState but non-zero for jumpCount");
					}
				}
				(iterExp->second[termState][jumpState]) /= iterTerm->second[termState];
				
				//(iterProb->second[termState][jumpState]) /= iterTerm->second[termState];
				
			}
		}
	}
}



//////////////////////////////////////////////////////////
//combined two characters into a combined state.
//For example. if the alphabet is {0,1,2} then the combined alphabet will be {0,1...8}.
//The states (terminalStart, terminalEnd) = (0,2) then combinedId = 2.
//The states (terminalStart, terminalEnd) = (1,2) then combinedId = 5. etc.
int simulateJumpsLargeAlphabet::getCombinedState(int terminalStart, int terminalEnd) const
{
	return (terminalStart * _pAlph->size() + terminalEnd);
}
int simulateJumpsLargeAlphabet::getStartId(int combinedState) const
{
	return combinedState / _pAlph->size();
}
int simulateJumpsLargeAlphabet::getEndId(int combinedState) const
{
	return combinedState % _pAlph->size();
}
//////////////////////////////////////////////////////////


MDOUBLE simulateJumpsLargeAlphabet::getExpectation(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId) 
{
	map<string, map_i2i2double >::iterator pos;
	if ((pos = _nodes2JumpsExp.find(nodeName)) == _nodes2JumpsExp.end())
	{
		string err="error in simulateJumpsLargeAlphabet::getExpectation: cannot find node "+nodeName;
		LOG(5, <<"names in _nodes2JumpsExp:"<<endl;);
		for (pos = _nodes2JumpsExp.begin(); pos != _nodes2JumpsExp.end(); ++pos)
			LOG(5, <<pos->first<<endl;);
		errorMsg::reportError(err);
	}
	int combinedTerminalState = getCombinedState(terminalStart, terminalEnd);
	int combinedJumpState = getCombinedState(fromId, toId);
	map_i2i2double::iterator iter2 = _nodes2JumpsExp[nodeName].find(combinedTerminalState);
	if (iter2 == pos->second.end())
		return 0;
	map<int, MDOUBLE>::iterator iter3 = _nodes2JumpsExp[nodeName][combinedTerminalState].find(combinedJumpState);
	if (iter3 == iter2->second.end())
		return 0;
	return iter3->second;
}


MDOUBLE simulateJumpsLargeAlphabet::getProb(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId) {
	//map <string, VVdouble>::iterator pos;
	//if ((pos = _nodes2JumpsProb.find(nodeName)) == _nodes2JumpsProb.end())
	//{
	//	string err="error in simulateJumpsLargeAlphabet::getProb: cannot find node "+nodeName;
	//	errorMsg::reportError(err);
	//}
	//int combinedTerminalState = getCombinedState(terminalStart, terminalEnd);
	//int combinedJumpState = getCombinedState(fromId, toId);
	//return (pos->second[combinedTerminalState][combinedJumpState]);
	return -1;
}