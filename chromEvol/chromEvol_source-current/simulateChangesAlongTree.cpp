#include "simulateChangesAlongTree.h"
#include "talRandom.h"
#include "matrixUtils.h"
#include "treeIt.h"
#include "chrNumberMng.h"
#include "chrNumberOptions.h"
#include "treeUtil.h"
#include <algorithm>


simulateChangesAlongTree::simulateChangesAlongTree(const tree& inTree, const stochasticProcess& sp, finiteIntAlphabet* pAlph)
: _tree(inTree), _sp(sp), _pAlph(pAlph)	
{
}

simulateChangesAlongTree::~simulateChangesAlongTree()
{
}

void simulateChangesAlongTree::init()
{
	//init the vector of waiting times. 
	_waitingTimeParams.clear();
	_waitingTimeParams.resize(_pAlph->size());
	int i, j;
	for (i = 0; i < _pAlph->size(); ++i)
	{
		_waitingTimeParams[i] = -_sp.dPij_dt(i, i, 0.0);
		
	}

	//init _jumpProbs. _jumpProbs[i][j] = Q[i][j] / -Q[i][i]
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
			string err = "error in simulateJumps::init(): sum probabilities is not 1 and equal to " + double2string(sum);
			errorMsg::reportError(err);
		}
	}
	int nodesNum = _tree.getNodesNum();
	_changesOccurred.clear();
	_changesOccurred.resize(nodesNum);
	for (int i=0; i<nodesNum; ++i)
		resizeMatrix(_changesOccurred[i], _pAlph->size(), _pAlph->size());
	_nodesContent.clear();
	_nodesContent.resize(nodesNum, 0);
}

sequenceContainer simulateChangesAlongTree::simulatePosition(){
	init();
	Vdouble freqs(_pAlph->size(),0.0);
	for (int i = 0; i< freqs.size(); ++i)
		freqs[i]=_sp.freq(i);
	int rootState = giveRandomState(_pAlph->size(), freqs);
	string name = _tree.getRoot()->name();
	sequence seq(_pAlph->fromInt(rootState),name, "", _tree.getRoot()->id(), _pAlph);
	_sc.add(seq);



	_nodesContent[_tree.getRoot()->id()] = rootState;
	simulateOnce(_tree.getRoot(),0,rootState,0);
	simulateOnce(_tree.getRoot(),0,rootState,1);
	if (_tree.getRoot()->getNumberOfSons() > 2)
		simulateOnce(_tree.getRoot(),0,rootState,2);
	calcEventsTypes();
	return _sc;
}

void simulateChangesAlongTree::simulateOnce(tree::nodeP curNode, 
											MDOUBLE disFromNode, 
											int previousContent, int whichSon){
		tree::nodeP sonNode = curNode->getSon(whichSon);
		MDOUBLE avgWaitingTime = 1.0 / _waitingTimeParams[previousContent];
		MDOUBLE timeTillChange = talRandom::rand_exp(avgWaitingTime);
		disFromNode += timeTillChange;
		int nextContent = giveRandomState(_pAlph->size(), previousContent, _jumpProbs);
		while (disFromNode < sonNode->dis2father()) {
			_changesOccurred[sonNode->id()][previousContent][nextContent]++;
			previousContent=nextContent;
			MDOUBLE avgWaitingTime = 1.0 / _waitingTimeParams[previousContent];
			MDOUBLE timeTillChange = talRandom::rand_exp(avgWaitingTime);
			disFromNode += timeTillChange;
			nextContent = giveRandomState(_pAlph->size(), nextContent, _jumpProbs);
		}
		while (disFromNode >= sonNode->dis2father()) {
			_nodesContent[sonNode->id()] = previousContent;
			string name = sonNode->name();
			sequence seq(_pAlph->fromInt(previousContent),name, "", sonNode->id(), _pAlph);
			_sc.add(seq);
			if (sonNode->isLeaf()) {
				//string name = sonNode->name();
				//sequence seq(int2string(previousContent),name, "", sonNode->id(), _pAlph);
				//_sc.add(seq);
				return;
			}
			for (int s = 1; s < sonNode->getNumberOfSons(); ++s)
                simulateOnce(sonNode, 0, previousContent, s);
			disFromNode-=sonNode->dis2father();
			curNode = sonNode;
			sonNode = curNode->getSon(0);
		}
		_changesOccurred[sonNode->id()][previousContent][nextContent]++;
		simulateOnce(curNode, disFromNode, nextContent, 0);
}

VVint simulateChangesAlongTree::getChangesForBranch(int nodeID){
	if (nodeID>_changesOccurred.size())
		errorMsg::reportError("error in simulateChangesAlongTree::getChangesForBranch, nodeID doesn't exist");
	return _changesOccurred[nodeID];
}

sequenceContainer simulateChangesAlongTree::toSeqDataWithoutInternalNodes() {
	sequenceContainer myseqData;
	int nextId = 0;
	for (sequenceContainer::taxaIterator it =_sc.taxaBegin(); it != _sc.taxaEnd(); ++it)
	{
		string nodeName = it->name();
		tree::nodeP theCurNode = _tree.findNodeByName(nodeName);
		if (theCurNode == NULL)
			errorMsg::reportError("could not find the specified name: " + nodeName);
		if (theCurNode->isInternal()) 
			continue;
		sequence tmpSeq = *it;
		tmpSeq.setID(nextId);
		myseqData.add(tmpSeq);
		++nextId;
	}
	return myseqData;
}

//type = gain/loss/dupl/demi/ or jumps to maxChr
int simulateChangesAlongTree::getTotalJumps(chrNumModel::jumpType type)
{
	int totalPos = _changesTypes.size() -1;
	switch (type)
	{
	case chrNumModel::GAIN_J:
	case chrNumModel::LOSS_J:
	case chrNumModel::DUPL_J:
	case chrNumModel::DEMI_J:
	case chrNumModel::BASE_J:
	case chrNumModel::MAX_CHR_J:
		return _changesTypes[totalPos][type];
		break;
	default:
		errorMsg::reportError("simulateChangesAlongTree::getTotalJumps: unknown type");
		return -1;
	}
}
#define THRESHOLD_R 0.000001
//fill _changesTypes: the number of gain/loss/dupl/demi/to_max_chr for each node
void simulateChangesAlongTree::calcEventsTypes()
{
	int maxChr = _pAlph->size();
	int nodesNum = _tree.getNodesNum();
	_changesTypes.resize(nodesNum+1);
	for (int n = 0; n < _changesTypes.size(); ++n)
		_changesTypes[n].resize(chrNumModel::JUMP_TYPE_MAX, 0);
	treeIterTopDownConst tIt(_tree);
	chrNumModel* pModel = static_cast<chrNumModel*>(_sp.getPijAccelerator()->getReplacementModel());
	int baseNumber = pModel->getBaseNumber(); 
	MDOUBLE treeLength = getSumOfBranchLengths(_tree);
	bool bIgnoreGain = false, bIgnoreLoss = false, bIgnoreBase = false, bIgnoreDupl = false, bIgnoreDemi = false;
	if (pModel->getGainR(0)*treeLength < THRESHOLD_R) //the multiplication by treeLength is to estimate the total expected number of transitions along the tree
		bIgnoreGain = true;
	if (pModel->getLossR(0) * treeLength < THRESHOLD_R)
		bIgnoreLoss = true;
	if (pModel->getBaseNumberR() * treeLength < THRESHOLD_R)
		bIgnoreBase = true;
	if (pModel->getDuplR(0) * treeLength < THRESHOLD_R)
		bIgnoreDupl = true;
	if (pModel->getDemiR(0) * treeLength < THRESHOLD_R)
		bIgnoreDemi = true;

	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		int chr;
		for (chr = 2; chr < maxChr-1; ++chr)
		{
			if (!bIgnoreGain)
				_changesTypes[nodeId][chrNumModel::GAIN_J] += _changesOccurred[nodeId][chr-1][chr];
			if (!bIgnoreLoss)
				_changesTypes[nodeId][chrNumModel::LOSS_J] += _changesOccurred[nodeId][chr-1][chr-2];
			if (!bIgnoreBase) {
				for (int toChr = chr+baseNumber; toChr < maxChr; toChr+=baseNumber) {
					_changesTypes[nodeId][chrNumModel::BASE_J] += _changesOccurred[nodeId][chr-1][toChr-1]; //base number jumps
					if (_changesOccurred[nodeId][chr - 1][toChr - 1])
						--_changesOccurred[nodeId][chr - 1][toChr - 1];
				}
			}
			_changesTypes[nodeId][chrNumModel::MAX_CHR_J] += _changesOccurred[nodeId][chr - 1][maxChr - 1]; //changes to maxChr
			//compute exp of a duplication 
			if (!bIgnoreDupl && chr <= maxChr/2)
			{
				_changesTypes[nodeId][chrNumModel::DUPL_J] += _changesOccurred[nodeId][chr-1][2*chr-1];
			}
			//else 
			//{
			//	_changesTypes[nodeId][2] += _changesOccurred[nodeId][chr-1][maxChr-1];
			//}
			//compute exp of a demi-duplication 
			if (!bIgnoreDemi &&  ((chr <= maxChr * 0.66) && (chr > 2)))
			{
				if (chr == 3)
				{
					_changesTypes[nodeId][chrNumModel::DEMI_J] += _changesOccurred[nodeId][2][4]; //count only 3->5
				}
				else if ((chr % 2) == 0)
				{
                    _changesTypes[nodeId][chrNumModel::DEMI_J] += _changesOccurred[nodeId][chr-1][static_cast<int>(1.5*chr-1)];
				}
				else
				{
					_changesTypes[nodeId][chrNumModel::DEMI_J] += _changesOccurred[nodeId][chr-1][static_cast<int>(1.5*chr-1)];
					_changesTypes[nodeId][chrNumModel::DEMI_J] += _changesOccurred[nodeId][chr-1][static_cast<int>(1.5*chr)];
				}
			}
		}
		//for (maxChr-1) there are only two possibilities : loss and jumps to maxChr
		chr = maxChr-1;
		if (!bIgnoreLoss)
			_changesTypes[nodeId][chrNumModel::LOSS_J] += _changesOccurred[nodeId][chr-1][chr-2];
		_changesTypes[nodeId][chrNumModel::MAX_CHR_J] += _changesOccurred[nodeId][chr-1][maxChr-1];
		//for maxChr there is only one possibility loss 
		chr = maxChr;
		if (!bIgnoreLoss)
			_changesTypes[nodeId][chrNumModel::LOSS_J] += _changesOccurred[nodeId][chr-1][chr-2]; //loss from max to max-1
		//add all results
		_changesTypes[nodesNum][chrNumModel::GAIN_J] += _changesTypes[nodeId][chrNumModel::GAIN_J];
		_changesTypes[nodesNum][chrNumModel::LOSS_J] += _changesTypes[nodeId][chrNumModel::LOSS_J];
		_changesTypes[nodesNum][chrNumModel::DUPL_J] += _changesTypes[nodeId][chrNumModel::DUPL_J];
		_changesTypes[nodesNum][chrNumModel::DEMI_J] += _changesTypes[nodeId][chrNumModel::DEMI_J];
		_changesTypes[nodesNum][chrNumModel::BASE_J] += _changesTypes[nodeId][chrNumModel::BASE_J];
		_changesTypes[nodesNum][chrNumModel::MAX_CHR_J] += _changesTypes[nodeId][chrNumModel::MAX_CHR_J];
	}	
}

void simulateChangesAlongTree::printEvents(ostream &expFile)
{
	treeIterTopDownConst tItr(_tree);
	//print gains
	expFile<<"#Nodes with GAIN events"<<endl;
	MDOUBLE sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _changesTypes[nodeId][chrNumModel::GAIN_J];
		if (_changesTypes[nodeId][chrNumModel::GAIN_J] > 0)
			expFile<<mynode->name()<<": "<<_changesTypes[nodeId][chrNumModel::GAIN_J]<<endl;
	}
	expFile<<"Total number of gain events: "<<sum<<endl;

	//print losses
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with LOSS events"<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _changesTypes[nodeId][chrNumModel::LOSS_J];
		if (_changesTypes[nodeId][chrNumModel::LOSS_J] > 0)
			expFile<<mynode->name()<<": "<<_changesTypes[nodeId][chrNumModel::LOSS_J]<<endl;
	}
	expFile<<"Total number of loss events: "<<sum<<endl;

	//print duplications
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with duplication events"<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _changesTypes[nodeId][chrNumModel::DUPL_J];
		if (_changesTypes[nodeId][chrNumModel::DUPL_J] > 0)
			expFile<<mynode->name()<<": "<<_changesTypes[nodeId][chrNumModel::DUPL_J]<<endl;
	}
	expFile<<"Total number of duplication events: "<<sum<<endl;

	//print demi-duplications
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with demi-duplication events"<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _changesTypes[nodeId][chrNumModel::DEMI_J];
		if (_changesTypes[nodeId][chrNumModel::DEMI_J] > 0)
			expFile<<mynode->name()<<": "<<_changesTypes[nodeId][chrNumModel::DEMI_J]<<endl;
	}
	expFile<<"Total number of demi-dulications events: "<<sum<<endl;
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;

	//print base-number transitions
	bool bBaseNumber = (IGNORE_PARAM != static_cast<chrNumModel*>(_sp.getPijAccelerator()->getReplacementModel())->getBaseNumberR());
	sum = 0.0;
	if (bBaseNumber)
	{
		expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
		expFile<<"#Nodes with base number transitions"<<endl<<endl;
		for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
		{
			if (mynode->isRoot())
				continue;
			int nodeId = mynode->id();
			sum += _changesTypes[nodeId][chrNumModel::BASE_J];
			if (_changesTypes[nodeId][chrNumModel::BASE_J] > 0)
				expFile<<mynode->name()<<": "<<_changesTypes[nodeId][chrNumModel::BASE_J]<<endl;
		}
		expFile<<"Total number of base number transitions: "<<sum<<endl;
		expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	}	


	//print transitions to max chromosome
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with transitions to max chromosome allowed"<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _changesTypes[nodeId][chrNumModel::MAX_CHR_J];
		if (_changesTypes[nodeId][chrNumModel::MAX_CHR_J] > 0)
			expFile<<mynode->name()<<": "<<_changesTypes[nodeId][chrNumModel::MAX_CHR_J]<<endl;
	}
	expFile<<"Total number of transitions to max chromosome: "<<sum<<endl;
	
	//for each leaf print the number of transitions (for each event type) from root to leaf
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	VVdouble nodesTransitions; //leavesTransitions[leafId][chrNumModel::jumpType]
	computeTransitionsFromRootToNodes(nodesTransitions);
	expFile<<"#Expected number of events from root to leaf"<<endl;
	expFile<<"LEAF"<<"\t";
	if (bBaseNumber)
		expFile<<"GAIN"<<"\t"<<"LOSS"<<"\t"<<"DUPLICATION"<<"\t"<<"DEMI-DUPLICATION"<<"\t"<<"BASE-NUMBER"<<endl;
	else
		expFile<<"GAIN"<<"\t"<<"LOSS"<<"\t"<<"DUPLICATION"<<"\t"<<"DEMI-DUPLICATION"<<endl;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (!mynode->isLeaf())
			continue;
		expFile<<mynode->name()<<"\t";
		int nodeId = mynode->id();
		if (bBaseNumber)
			expFile<<nodesTransitions[nodeId][chrNumModel::GAIN_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::LOSS_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::DUPL_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::DEMI_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::BASE_J]<<endl;
		else
			expFile<<nodesTransitions[nodeId][chrNumModel::GAIN_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::LOSS_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::DUPL_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::DEMI_J]<<endl;

		//expFile<<nodesTransitions[nodeId][chrNumModel::GAIN_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::LOSS_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::DUPL_J]<<"\t"<<nodesTransitions[nodeId][chrNumModel::DEMI_J]<<endl;
	}
}

void simulateChangesAlongTree::computeTransitionsFromRootToNodes(VVdouble& nodesTransitions)
{
	nodesTransitions.clear();
	resizeMatrix(nodesTransitions, _tree.getNodesNum(), chrNumModel::JUMP_TYPE_MAX);
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		int fatherId = mynode->father()->id();
		nodesTransitions[nodeId][chrNumModel::GAIN_J] = nodesTransitions[fatherId][chrNumModel::GAIN_J] + _changesTypes[nodeId][chrNumModel::GAIN_J];
		nodesTransitions[nodeId][chrNumModel::LOSS_J] = nodesTransitions[fatherId][chrNumModel::LOSS_J] + _changesTypes[nodeId][chrNumModel::LOSS_J];
		nodesTransitions[nodeId][chrNumModel::DUPL_J] = nodesTransitions[fatherId][chrNumModel::DUPL_J] + _changesTypes[nodeId][chrNumModel::DUPL_J];
		nodesTransitions[nodeId][chrNumModel::DEMI_J] = nodesTransitions[fatherId][chrNumModel::DEMI_J] + _changesTypes[nodeId][chrNumModel::DEMI_J];
		nodesTransitions[nodeId][chrNumModel::BASE_J] = nodesTransitions[fatherId][chrNumModel::BASE_J] + _changesTypes[nodeId][chrNumModel::BASE_J];
		nodesTransitions[nodeId][chrNumModel::MAX_CHR_J] = nodesTransitions[fatherId][chrNumModel::MAX_CHR_J] + _changesTypes[nodeId][chrNumModel::MAX_CHR_J];
	}
}

//return // for each node, the expected number of changes for [0-4]: gain, loss, dupl, halfDupl, baseNumber
VVdouble simulateChangesAlongTree::generateStochasticMapping(const sequenceContainer& inSc, int iterations, const VVVdouble& jointPost, const suffStatGlobalHomPos& sscUp)
{
	init();
	int nodesNum = _tree.getNodesNum();
	_changesTypes.resize(nodesNum+1);
	for (int n = 0; n < _changesTypes.size(); ++n)
		_changesTypes[n].resize(chrNumModel::JUMP_TYPE_MAX, 0);

	Vdouble rootProbs = jointPost[_tree.getRoot()->id()][0];
	VVVdouble samplingProbs = getSamplingProbs(sscUp);

	//store all jumps that occured across all mappings in _changesTypes
	for (int i = 0; i < iterations; ++i)
	{
		//3. sample ancestrals and store in _nodeContants
        sampleAncestrals(rootProbs, samplingProbs);
        //4. simulate mutations based on ancestrals
        sampleMutationsGivenAncestrals();
	}

	//convert from the number of events that occured to expectations
	VVdouble res;
	resizeMatrix(res, _tree.getNodesNum()+1, chrNumModel::JUMP_TYPE_MAX);
	MDOUBLE gainExp = 0.0, lossExp = 0.0, duplExp = 0.0, halfDuplExp = 0.0, maxChrExp = 0.0, baseExp=0.0; 
	for (int n = 0; n < res.size(); ++n) {
		res[n][chrNumModel::GAIN_J] = _changesTypes[n][chrNumModel::GAIN_J]/static_cast<MDOUBLE>(iterations);
		res[n][chrNumModel::LOSS_J] = _changesTypes[n][chrNumModel::LOSS_J]/static_cast<MDOUBLE>(iterations);
		res[n][chrNumModel::DUPL_J] = _changesTypes[n][chrNumModel::DUPL_J]/static_cast<MDOUBLE>(iterations);
		res[n][chrNumModel::DEMI_J] = _changesTypes[n][chrNumModel::DEMI_J]/static_cast<MDOUBLE>(iterations);
		res[n][chrNumModel::BASE_J] = _changesTypes[n][chrNumModel::BASE_J]/static_cast<MDOUBLE>(iterations);
		res[n][chrNumModel::MAX_CHR_J] = _changesTypes[n][chrNumModel::MAX_CHR_J]/static_cast<MDOUBLE>(iterations);
		gainExp += res[n][chrNumModel::GAIN_J];
		lossExp += res[n][chrNumModel::LOSS_J];
		duplExp += res[n][chrNumModel::DUPL_J];
		halfDuplExp += res[n][chrNumModel::DEMI_J];
		baseExp += res[n][chrNumModel::BASE_J];
		maxChrExp += res[n][chrNumModel::MAX_CHR_J];
	}
	res[res.size()-1][chrNumModel::GAIN_J] = gainExp;
	res[res.size()-1][chrNumModel::LOSS_J] = lossExp;
	res[res.size()-1][chrNumModel::DUPL_J] = duplExp;
	res[res.size()-1][chrNumModel::DEMI_J] = halfDuplExp;
	res[res.size()-1][chrNumModel::BASE_J] = baseExp;
	res[res.size()-1][chrNumModel::MAX_CHR_J] = maxChrExp;

	return res;
}

//perform stochastic mapping. The ancestral states are already stored in _nodesContent. 
//The results of each mapping are added to _changesTypes
void simulateChangesAlongTree::sampleMutationsGivenAncestrals() {
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		int nodeId = mynode->id();
		if (mynode->isRoot())
			continue;
		else
		{
			int fatherState = _nodesContent[mynode->father()->id()];
			int sonState = _nodesContent[nodeId];
			vector<chrNumModel::jumpType> res = sampleMutationsGivenAncestralsPerBranch(fatherState, sonState, mynode->dis2father());
			for (int i = 0; i < res.size(); ++i)
				_changesTypes[nodeId][res[i]] += 1;
		}
	}
}



vector<chrNumModel::jumpType> simulateChangesAlongTree::sampleMutationsGivenAncestralsPerBranch(int fatherState, int sonState, MDOUBLE branchLength)
{
	int maxIterNum = 1000;
	for (int i = 0; i < maxIterNum; ++i)
	{
		vector<chrNumModel::jumpType> tryMapping;
		MDOUBLE disFromNode = 0.0;
		int curState = fatherState;
		MDOUBLE avgWaitingTime = 1.0 / _waitingTimeParams[curState];
		MDOUBLE timeTillChange = talRandom::rand_exp(avgWaitingTime);

		if (fatherState != sonState)
		{ //draw jump time given a change has occured
			MDOUBLE u = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
            MDOUBLE tmp = u * (1.0 - exp(branchLength * -_waitingTimeParams[curState]));
            timeTillChange =  -log(1.0 - tmp) / _waitingTimeParams[curState];
			assert (timeTillChange < branchLength);
		}
        
		while (disFromNode + timeTillChange < branchLength)
		{
            //a jump occured but not passed the whole branch. Add the current state and time to branch history and draw next state
			////pair<int, MDOUBLE> oneJump(curState, timeTillChange);
			////tryMapping.push_back(oneJump);
			////curState = giveRandomState(_sp.alphabetSize(), curState, _jumpProbs);

			int nextState = giveRandomState(_sp.alphabetSize(), curState, _jumpProbs);
			tryMapping.push_back(getJumpType(curState, nextState));
			curState = nextState;
			disFromNode += timeTillChange;			
			avgWaitingTime = 1.0 / _waitingTimeParams[curState];
            timeTillChange = talRandom::rand_exp(avgWaitingTime);
       	}
		//the last jump passed the length of the branch
		if (curState != sonState)
		{ //simulation failed
			continue;
		}
		else
		{
			//MDOUBLE timeOfJump = branchLength - disFromNode;
			//pair<int, MDOUBLE> oneJump(curState, timeOfJump);
			//tryMapping.push_back(oneJump);
			return tryMapping;
		}
	}
	//all simulations failed - estimate events heuristically
	string err = "could not produce simulations with father = " +  int2string(_pAlph->id2Count(fatherState)) + " son = " + int2string(_pAlph->id2Count(sonState)) + "branch length = " + double2string(branchLength);
	LOGnOUT(6, <<err<<endl;);
	vector<chrNumModel::jumpType> res;
	Vdouble estEvents = chrNumberMng::estimateEvents(fatherState, sonState, &_sp, _pAlph);
	for (; estEvents[chrNumModel::GAIN_J] > 0.5; --estEvents[chrNumModel::GAIN_J])
		res.push_back(chrNumModel::GAIN_J);
	for (;estEvents[chrNumModel::LOSS_J] > 0.5; --estEvents[chrNumModel::LOSS_J])
		res.push_back(chrNumModel::LOSS_J);
	for (; estEvents[chrNumModel::DUPL_J] > 0.5; --estEvents[chrNumModel::DUPL_J])
		res.push_back(chrNumModel::DUPL_J); 
	for (; estEvents[chrNumModel::DEMI_J] > 0.5; --estEvents[chrNumModel::DEMI_J])
		res.push_back(chrNumModel::DEMI_J);
	for (; estEvents[chrNumModel::BASE_J] > 0.5; --estEvents[chrNumModel::BASE_J])
		res.push_back(chrNumModel::BASE_J);
	return res;
}




VVVdouble simulateChangesAlongTree::getSamplingProbs(const suffStatGlobalHomPos& sscUp) const
{
	//compute sampling probabilities for each node: P(letter at Node is j | Data & letter at Father(Node) is i)
	//P(Node=j | Data, Father=i) = Up[node][j] * Pij(t) / Sum_k {Up[node][k] * Pik(t) }
	VVVdouble samplingProbs(_tree.getNodesNum());//samplingProbs[nodeId][fatherState][sonState]
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		if (mynode->isRoot())
			continue;
		resizeMatrix(samplingProbs[mynode->id()], _sp.alphabetSize(), _sp.alphabetSize());
		for (int fatherState = 0; fatherState < _sp.alphabetSize(); ++fatherState)
		{
			MDOUBLE sum = 0.0; 
			for (int sonState = 0; sonState < _sp.alphabetSize(); ++sonState)
			{
				doubleRep tmp = _sp.Pij_t(fatherState, sonState, mynode->dis2father()) * sscUp.get(mynode->id(),sonState);
				samplingProbs[mynode->id()][fatherState][sonState] = convert(tmp);
				sum += convert(tmp);
			}
			for (int sonState = 0; sonState < _sp.alphabetSize(); ++sonState)
				samplingProbs[mynode->id()][fatherState][sonState] /= sum;
		}
	}
	return samplingProbs;
}


void simulateChangesAlongTree::sampleAncestrals(const Vdouble& rootProbs, const VVVdouble& samplingProbs)
{
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		int nodeId = mynode->id();
		if (mynode->isRoot())
		{
			int rootState = giveRandomState(_sp.alphabetSize(), rootProbs);
			_nodesContent[nodeId] = rootState;
		}
		else
		{
			int fatherState = _nodesContent[mynode->father()->id()];
			int sonState = giveRandomState(_sp.alphabetSize(), fatherState, samplingProbs[nodeId]);
			_nodesContent[nodeId] = sonState;
		}
	}
}


chrNumModel::jumpType simulateChangesAlongTree::getJumpType(int curState, int nextState) {
	int curChrom = _pAlph->id2Count(curState);
	int nextChrom = _pAlph->id2Count(nextState);
	chrNumModel* pModel = static_cast<chrNumModel*>(_sp.getPijAccelerator()->getReplacementModel());
	int baseNumber = pModel->getBaseNumber(); 
	int dif = nextChrom - curChrom;
	if (nextChrom == _pAlph->max())
		return chrNumModel::MAX_CHR_J;
	if ((dif == 1) && (pModel->getGainR(curState) > 0))
		return chrNumModel::GAIN_J;
	else if ((dif == -1) && (pModel->getLossR(curState) > 0))
		return chrNumModel::LOSS_J;
	else if ((curChrom*2 == nextChrom) && (pModel->getDuplR(curState) > 0))
		return chrNumModel::DUPL_J;
	else if ((pModel->getBaseNumber() > chrNumberOptions::_minBaseTransition) && ((dif % baseNumber) == 0))
		return chrNumModel::BASE_J;
	else if (pModel->getDemiR(curState) > 0) {
		if (curChrom % 2 == 0) {
			if (curChrom*1.5 == nextChrom)
				return chrNumModel::DEMI_J;
		}
		else if ((int(curChrom*1.5) == nextChrom) || (int(curChrom*1.5) == nextChrom-1))
			return chrNumModel::DEMI_J;
	}
	string err = "error in simulateChangesAlongTree::getJumpType: could not infer jump type from " + int2string(curChrom) + "to "  + int2string(nextChrom);
	errorMsg::reportError(err);
	return chrNumModel::MAX_CHR_J;

}