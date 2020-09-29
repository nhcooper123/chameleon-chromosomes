#include "evalPPDistMng.h"
#include "chrNumberOptions.h"
#include "treeIt.h"
#include "treeUtil.h"
#include "logFile.h"
#include <cstdlib>

evalPPDistMng::evalPPDistMng() {
	init();
}


evalPPDistMng::~evalPPDistMng() {

}


void evalPPDistMng::init() {
	string logFileName("");
	if (chrNumberOptions::_logFile.size() > 0)
         logFileName = chrNumberOptions::_outDir + "//" + chrNumberOptions::_logFile;
	myLog::setLog(logFileName, chrNumberOptions::_logValue);

}

//the tree height threshold is stored in chrNumberOptions::_optimizePointsNum;
void evalPPDistMng::run()
{
	getTree();
	readExpectations();
	calculateEdgeHeights();
	calcExpectedAndObservedCounts();
	calcChiStat();
}

void evalPPDistMng::calcChiStat() 
{
	MDOUBLE observed1 = _obsChangesSet1[chrNumModel::DUPL_J] + _obsChangesSet1[chrNumModel::DEMI_J];
	MDOUBLE observed2 = _obsChangesSet2[chrNumModel::DUPL_J] + _obsChangesSet2[chrNumModel::DEMI_J];
	MDOUBLE expected1 = _expChangesSet1[chrNumModel::DUPL_J] + _expChangesSet1[chrNumModel::DEMI_J];
	MDOUBLE expected2 = _expChangesSet2[chrNumModel::DUPL_J] + _expChangesSet2[chrNumModel::DEMI_J];
	MDOUBLE diffsq1 = (observed1 - expected1)*(observed1 - expected1)/expected1;
	MDOUBLE diffsq2 = (observed2 - expected2)*(observed2 - expected2)/expected2;
	MDOUBLE chi = diffsq1 + diffsq2;
	LOGnOUT(3,<<chrNumberOptions::_dataFile<<"DUPL+DEMI:"<<endl<<" chi = "<<chi<<" observed1 = "<<observed1<<" exp1 = "<<expected1<<" observed2 = "<<observed2<<" exp2 = "<<expected2<< endl);
	observed1 = _obsChangesSet1[chrNumModel::DUPL_J];
	observed2 = _obsChangesSet2[chrNumModel::DUPL_J];
	expected1 = _expChangesSet1[chrNumModel::DUPL_J];
	expected2 = _expChangesSet2[chrNumModel::DUPL_J];
	diffsq1 = (observed1 - expected1)*(observed1 - expected1)/expected1;
	diffsq2 = (observed2 - expected2)*(observed2 - expected2)/expected2;
	chi = diffsq1 + diffsq2;
	LOGnOUT(3,<<chrNumberOptions::_dataFile<<"DUPL:"<<endl<<" chi = "<<chi<<" observed1 = "<<observed1<<" exp1 = "<<expected1<<" observed2 = "<<observed2<<" exp2 = "<<expected2<< endl);
	observed1 = _obsChangesSet1[chrNumModel::DEMI_J];
	observed2 = _obsChangesSet2[chrNumModel::DEMI_J];
	expected1 = _expChangesSet1[chrNumModel::DEMI_J];
	expected2 = _expChangesSet2[chrNumModel::DEMI_J];
	diffsq1 = (observed1 - expected1)*(observed1 - expected1)/expected1;
	diffsq2 = (observed2 - expected2)*(observed2 - expected2)/expected2;
	chi = diffsq1 + diffsq2;
	LOGnOUT(3,<<chrNumberOptions::_dataFile<<"DEMI:"<<endl<<" chi = "<<chi<<" observed1 = "<<observed1<<" exp1 = "<<expected1<<" observed2 = "<<observed2<<" exp2 = "<<expected2<< endl);
	observed1 = _obsChangesSet1[chrNumModel::GAIN_J];
	observed2 = _obsChangesSet2[chrNumModel::GAIN_J];
	expected1 = _expChangesSet1[chrNumModel::GAIN_J];
	expected2 = _expChangesSet2[chrNumModel::GAIN_J];
	diffsq1 = (observed1 - expected1)*(observed1 - expected1)/expected1;
	diffsq2 = (observed2 - expected2)*(observed2 - expected2)/expected2;
	chi = diffsq1 + diffsq2;
	LOGnOUT(3,<<chrNumberOptions::_dataFile<<"GAIN:"<<endl<<" chi = "<<chi<<" observed1 = "<<observed1<<" exp1 = "<<expected1<<" observed2 = "<<observed2<<" exp2 = "<<expected2<< endl);
	observed1 = _obsChangesSet1[chrNumModel::LOSS_J];
	observed2 = _obsChangesSet2[chrNumModel::LOSS_J];
	expected1 = _expChangesSet1[chrNumModel::LOSS_J];
	expected2 = _expChangesSet2[chrNumModel::LOSS_J];
	diffsq1 = (observed1 - expected1)*(observed1 - expected1)/expected1;
	diffsq2 = (observed2 - expected2)*(observed2 - expected2)/expected2;
	chi = diffsq1 + diffsq2;
	LOGnOUT(3,<<chrNumberOptions::_dataFile<<"LOSS:"<<endl<<" chi = "<<chi<<" observed1 = "<<observed1<<" exp1 = "<<expected1<<" observed2 = "<<observed2<<" exp2 = "<<expected2<< endl);

}

void evalPPDistMng::calcExpectedAndObservedCounts()
{
	if (_cutoffs.size() == 1)
		calcExpectedAndObservedCountsLeaves();
	_expChangesSet1.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	_expChangesSet2.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	_obsChangesSet1.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	_obsChangesSet2.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	int totalId = _expChanges.size()-1;

	MDOUBLE treeLength = getSumOfBranchLengths(_tree);
	//calculate the observed and expected number of events along each edge that belongs to the first or second set
	treeIterDownTopConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isRoot())
			continue;
		int id = mynode->id();
		MDOUBLE lowerBound = _edgeHeights[id].first;
		MDOUBLE upperBound = _edgeHeights[id].second;
		//check if edge belongs to first set
		if (lowerBound < _cutoffs[0])
			lowerBound = _cutoffs[0];
		if (upperBound > _cutoffs[1])
			upperBound = _cutoffs[1];
		if (upperBound > lowerBound) {
			//edge belongs to set1 (but maybe only a fraction of it)
			MDOUBLE edgeLenInSet = upperBound-lowerBound;
			MDOUBLE edgeFracInSet = edgeLenInSet / mynode->dis2father();
			for (int type = 0; type < 4; ++ type) {
				_expChangesSet1[type] += (edgeLenInSet / treeLength) * _expChanges[totalId][type];
				_obsChangesSet1[type] += edgeFracInSet * _expChanges[id][type];
			}
		}
		//check if edge belongs to second set
		lowerBound = _edgeHeights[id].first;
		upperBound = _edgeHeights[id].second;
		if (lowerBound < _cutoffs[2])
			lowerBound = _cutoffs[2];
		if (upperBound > _cutoffs[3])
			upperBound = _cutoffs[3];
		if (upperBound > lowerBound) {
			//edge belongs to set1 (but maybe only a fraction of it)
			MDOUBLE edgeLenInSet = upperBound-lowerBound;
			MDOUBLE edgeFracInSet = edgeLenInSet / mynode->dis2father();
			for (int type = 0; type < 4; ++ type) {
				_expChangesSet2[type] += (edgeLenInSet / treeLength) * _expChanges[totalId][type];
				_obsChangesSet2[type] += edgeFracInSet * _expChanges[id][type];
			}
		}
	}
}

//chech external edges versus internal
void evalPPDistMng::calcExpectedAndObservedCountsLeaves()
{
	_expChangesSet1.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	_expChangesSet2.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	_obsChangesSet1.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	_obsChangesSet2.resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	int totalId = _expChanges.size()-1;

	MDOUBLE treeHeight = _cutoffs[0];
	MDOUBLE maxExternalBranch = treeHeight * 0.5;
	MDOUBLE treeLenthOrig = getSumOfBranchLengths(_tree);
	MDOUBLE treeLength = adjustExpectationsAndTreeLengthExcludingLongExternals(_tree, maxExternalBranch);
	//calculate the observed and expected number of events along each edge that belongs to the first or second set
	treeIterDownTopConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isRoot())
			continue;
		int id = mynode->id();
		//edge belongs to first set (=external) only if is shorter than 1/2 treeHeight
		if ((mynode->isLeaf()) && (mynode->dis2father() < maxExternalBranch))
		{
			for (int type = 0; type < 4; ++type) {
				_expChangesSet1[type] += (mynode->dis2father() / treeLength) * _expChanges[totalId][type];
				_obsChangesSet1[type] += _expChanges[id][type];
			}
		}
		else if (mynode->isInternal())
		{
			for (int type = 0; type < 4; ++ type) {
				_expChangesSet2[type] += (mynode->dis2father() / treeLength) * _expChanges[totalId][type];
				_obsChangesSet2[type] += _expChanges[id][type];
			}
		}
	}
}

void evalPPDistMng::calculateEdgeHeights() {
	_edgeHeights.resize(_tree.getNodesNum()); 
	treeIterDownTopConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isRoot())
			continue;
		else if (mynode->isLeaf()) {
			MDOUBLE lower = 0.0;
			MDOUBLE upper = mynode->dis2father();
			_edgeHeights[mynode->id()] = pair<MDOUBLE, MDOUBLE>(lower, upper);
		}
		else {
			if (mynode->getNumberOfSons() < 2)
				errorMsg::reportError("internal node " + mynode->name() + " must have at least 2 sons");
			int sonId = mynode->getSon(0)->id();
			MDOUBLE lower = _edgeHeights[sonId].second; // get the upper height for the son's edge;
			MDOUBLE lower2 = _edgeHeights[mynode->getSon(1)->id()].second;
			if (! DEQUAL(lower, lower2, 0.001))
				errorMsg::reportError("In an ultrametric tree, the edges from the 2 sons must end up with the same height. father name = " + mynode->name());
			MDOUBLE upper = lower + mynode->dis2father();
			_edgeHeights[mynode->id()] = pair<MDOUBLE, MDOUBLE>(lower, upper);
		}
	}
}


void evalPPDistMng::getTree()
{
	if (!(chrNumberOptions::_treeFile.empty()))
	{
		_tree = tree(chrNumberOptions::_treeFile);
	}

	//if tree is not rooted then user must speicify the rooting position
	if(_tree.getRoot()->getNumberOfSons() > 2) 
	{
		if (chrNumberOptions::_rootAt == "")
            errorMsg::reportError("The input tree is not rooted.");
		else {
			tree::nodeP myroot = _tree.findNodeByName(chrNumberOptions::_rootAt); //returns NULL if not found
			if (myroot){
				_tree.rootAt(myroot);
			}
			else
				errorMsg::reportError("could not find root " + chrNumberOptions::_rootAt);

		}
	}
	_tree.isUltrametric(0.05, true);
	vector<tree::nodeP> nodes;
	_tree.getAllLeaves(nodes, _tree.getRoot());
	MDOUBLE treeHeight = getDistanceFromNode2ROOT(nodes[0]);	
	_cutoffs.resize(chrNumberOptions::_optimizePointsNum.size());
	if (_cutoffs.size() == 1) {
		if (chrNumberOptions::_optimizePointsNum[0] != 100)
            errorMsg::reportError("to test external versus internal nodes _optimizePointsNum should be set to 100");
	}
	else if (chrNumberOptions::_optimizePointsNum.size() != 4)
		errorMsg::reportError("_optimizePointsNum should contain 4 points: lowerHeightSet1,upperHeightSet1,lowerHeightSet2,upperHeightSet2");

	for (int i = 0; i < _cutoffs.size(); ++i) {
		MDOUBLE x = treeHeight * static_cast<MDOUBLE>(chrNumberOptions::_optimizePointsNum[i]) / 100.0;
		_cutoffs[i] = treeHeight * static_cast<MDOUBLE>(chrNumberOptions::_optimizePointsNum[i]) / 100.0;
	}

	string outFileName = chrNumberOptions::_outDir + "//" + "allNodes.tree";
	_tree.output(outFileName, tree::PHYLIP, true);
}

void evalPPDistMng::readExpectations()
{
	int nodesNum = _tree.getNodesNum();
	_expChanges.clear();
	_expChanges.resize(nodesNum+1); //last place is total across tree
	for (int n = 0; n < _expChanges.size(); ++n)
		_expChanges[n].resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
		

	ifstream expFile(chrNumberOptions::_dataFile.c_str());
	vector<string> lines;
	putFileIntoVectorStringArray(expFile, lines);
	expFile.close();
	if (lines.empty()){
		errorMsg::reportError("unable to open file, or file is empty in readExpectations()");
	}
	
	vector<string>::const_iterator it1 = lines.begin();
	int localid=0;
	while ( ( (*it1).find("#ALL EVENTS EXPECTATIONS PER NODE")  == -1) && (it1 != lines.end()))
		++it1;
	if (it1 == lines.end())
		errorMsg::reportError("cannot find expectations per node");
	++it1;
	++it1;
	treeIterTopDownConst tItr(_tree);
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	//for (; ((*it1).find("#+++++")  == -1); ++it1)
	{
		//read tree node
		if (mynode->isRoot())
			continue;
		string treeNodeName = mynode->name();
		int nodeId = mynode->id();
		//read exp file
		vector<string> strVec;
		splitString(*it1, strVec, "\t "); 
		if (strVec.size() != 5)
			errorMsg::reportError("unable to read line");
		string expNodeName = strVec[0];
		if (expNodeName != treeNodeName)
			errorMsg::reportError("tree node is " + treeNodeName + " while exp file name is " + expNodeName);
		_expChanges[nodeId][chrNumModel::GAIN_J] = string2double(strVec[1]);
		_expChanges[nodeId][chrNumModel::LOSS_J] = string2double(strVec[2]);
		_expChanges[nodeId][chrNumModel::DUPL_J] = string2double(strVec[3]);
		_expChanges[nodeId][chrNumModel::DEMI_J] = string2double(strVec[4]);
		_expChanges[nodesNum][chrNumModel::GAIN_J] += _expChanges[nodeId][chrNumModel::GAIN_J];
		_expChanges[nodesNum][chrNumModel::LOSS_J] += _expChanges[nodeId][chrNumModel::LOSS_J];
		_expChanges[nodesNum][chrNumModel::DUPL_J] += _expChanges[nodeId][chrNumModel::DUPL_J];
		_expChanges[nodesNum][chrNumModel::DEMI_J] += _expChanges[nodeId][chrNumModel::DEMI_J];
		++it1;
	}
}

//calculate sum of branch lengths beside those external branches that are longer than the given cutoff
//in addition - subtract the corresponding expectations of the long branches
MDOUBLE evalPPDistMng::adjustExpectationsAndTreeLengthExcludingLongExternals(const tree &tr, MDOUBLE maxExternalBranch) {
	int nodesNum = _tree.getNodesNum();
	treeIterDownTopConst tIt(tr);
	MDOUBLE sum = 0;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!mynode->isRoot()){
			int nodeId = mynode->id();
			if (mynode->isLeaf()) {
				if (mynode->dis2father() > maxExternalBranch) {
					_expChanges[nodesNum][chrNumModel::GAIN_J] -= _expChanges[nodeId][chrNumModel::GAIN_J];
					_expChanges[nodesNum][chrNumModel::LOSS_J] -= _expChanges[nodeId][chrNumModel::LOSS_J];
					_expChanges[nodesNum][chrNumModel::DUPL_J] -= _expChanges[nodeId][chrNumModel::DUPL_J];
					_expChanges[nodesNum][chrNumModel::DEMI_J] -= _expChanges[nodeId][chrNumModel::DEMI_J];
				}
				else 
					sum+=mynode->dis2father();
			}
			else 
				sum+=mynode->dis2father();
		}
	}
	return sum;
}

