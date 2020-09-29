#include "chrNumberMng.h"
#include "chrNumberOptions.h"
#include "finiteIntAlphabet.h"
#include "treeIt.h"
#include "matrixUtils.h"
#include "trivialAccelerator.h"
#include "uniDistribution.h"
#include "sequence.h"
#include "someUtil.h"
#include "chrCountFormat.h"
#include "jcDistance.h"
#include "distanceTable.h"
#include "nj.h"
#include "treeUtil.h"
#include "fastaFormat.h"
#include "generalChrNumModel.h"
#include "chrNumOptimizer.h"
#include "computePijComponent.h"
#include "seqContainerTreeMap.h"
#include "computeUpAlgChrNum.h"
#include "computeDownAlg.h"
#include "simulateTree.h"
#include "evaluateCharacterFreq.h"
#include "simulateJumpsLargeAlphabet.h"
#include "logFile.h"
#include "dataUtils.h"
#include "Parameters.h"
#include "simulateChangesAlongTree.h"
#include "ConversionUtils.h"
#include <time.h>
#include <cstdlib>
#include <algorithm>
#include <cmath>

#define MAX_CHR_FOR_SIMULATIONS 2000
#define MAX_CHR_POSSIBLE 500
#define EXP_THRESHOLD 0.5
#define EXP_THRESHOLD2 0.9
#define EXP_THRESHOLD_MAX_CHR 0.05
#define THRESHOLD_R 0.000001

chrNumberMng::chrNumberMng() {
	init();
}


chrNumberMng::~chrNumberMng() {
	clear();
}

void chrNumberMng::clear() { 
	if (_pSp != NULL)
	{
		delete _pSp;
		_pSp = NULL;
	}
	if (_pAlph != NULL)
	{
		delete _pAlph;
		_pAlph = NULL;
	}
}


void chrNumberMng::init() {
	_pSp = NULL;
	_bestLL = 0;
	createDir("", chrNumberOptions::_outDir);
	string logFileName("");
	if (chrNumberOptions::_logFile.size() > 0)
         logFileName = chrNumberOptions::_outDir + "//" + chrNumberOptions::_logFile;
	myLog::setLog(logFileName, chrNumberOptions::_logValue);
	int max = abs(chrNumberOptions::_maxChrNum);
	if (max < 2)
		max = 2;
	_pAlph = new finiteIntAlphabet(1,max);
	if ((chrNumberOptions::_branchModelType == chrNumModel::COMBINED) || (chrNumberOptions::_branchModelType == chrNumModel::COMBINED_INTERNALS)) {
		if (chrNumberOptions::_scaleBranch < 0)
			errorMsg::reportError("error in chrNumberMng: COMBINED gradual and speciational model with a negative scale factor");
    }
}

void chrNumberMng::validateData()
{
    printRunInfo();
	getStartingData();
	pair<int,int> chrRange = getChrRangeInData(_sc);
	getTree();
	initStartingModel(_pAlph, chrRange.second - chrRange.first);
	//this is just to make sure all names in tree and sc are the same
	compLogLikelihood(_pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
	if ((chrNumberOptions::_smIter > 0) && (chrNumberOptions::_simulationsNum >  0)) {
		errorMsg::reportError("cannot perform both types of simulations to estimate expectations. Set one of the parameter _smIter or _simulationsNum to zero!");
	}
}

void chrNumberMng::run()
{
	time_t t1;
	time(&t1);
	time_t t2;

    printRunInfo();
	getStartingData();
	pair<int,int> minmax = getChrRangeInData(_sc);
	int chrRange;
	if 	(chrNumberOptions::_maxBaseTransition > 0)
		chrRange = chrNumberOptions::_maxBaseTransition;
	else {
		chrRange = minmax.second - minmax.first;
		if (chrNumberOptions::_baseNumber > chrRange)
			chrRange = chrNumberOptions::_baseNumber +1 ;
	}
	initStartingModel(_pAlph, chrRange);
	getTree();
	_bestLL = optimizeParams(true);
	LOGnOUT(5,<<"Inferring ancestral states"<<endl);
	VVVdouble jointPost;
	suffStatGlobalHomPos sscUp;
	computeJointPosterior(jointPost, sscUp);
	inferAncestralStates(jointPost);
	LOG(5,<<"Computing expectations"<<endl);
	computeChangeProbAndExp(jointPost, sscUp);
	LOG(5,<<"Printing results"<<endl);
	printResults(jointPost);
	time(&t2);
	LOGnOUT(5,<<endl<<"TOTAL RUNNING TIME = "<<static_cast<int>(t2-t1)<<endl);
}


void chrNumberMng::runSimulationsAndInference()
{
	printRunInfo();
	int iterNum = chrNumberOptions::_simulationsIter - chrNumberOptions::_startSimulationsIter;
	bool bJumpSimulation = !chrNumberOptions::_simulationsJumpsStats.empty();
	if (bJumpSimulation && (chrNumberOptions::_simulationsNum < 2))
		errorMsg::reportError("performing jumpSimulations but cannot make inference");
	if (!bJumpSimulation)
		chrNumberOptions::_simulationsNum = 0;
	//1a. init model with simulated params - this is actually just for printing
	finiteIntAlphabet simAlph(1,chrNumberOptions::_maxChrNumForSimulations);
	initStartingModel(&simAlph, chrNumberOptions::_maxBaseTransition);
	stochasticProcess* pSimSp = _pSp->clone();
	delete _pSp;
	_pSp = NULL;
	chrNumberOptions::_freqFile = "";
	chrNumModel* pSimModel = static_cast<chrNumModel*>(pSimSp->getPijAccelerator()->getReplacementModel());
	vector<paramObj> paramVec;
	pSimModel->getParams(paramVec);
	//1b. init vectors to store results 
	VVdouble simParams(paramVec.size()), optParams(paramVec.size());
	Vdouble simPostRoot(iterNum), simMlAllInternals(iterNum), simPostAllInternals(iterNum);
	Vint simMlRoot(iterNum), optMlRoot(iterNum);
	Vdouble optPostRoot(iterNum), optMlAllInternals(iterNum), optPostAllInternals(iterNum);  
	Vdouble simLLVec(iterNum), optLLVec(iterNum);
	for (int i = 0; i < simParams.size(); ++i)
	{
		simParams[i].resize(iterNum);
		optParams[i].resize(iterNum);
	}
	VVdouble optTotalJumps(chrNumModel::JUMP_TYPE_MAX), diffTotalJumps(chrNumModel::JUMP_TYPE_MAX); //diff is the difference between the real (simTotalJumps) and inferred number of jumps
	VVint simTotalJumps(chrNumModel::JUMP_TYPE_MAX);
	
	//2. open global output file and print headings 
	string baseOutDirName = chrNumberOptions::_outDir;
	string statsFileName = baseOutDirName + "//stats.txt";
	ofstream statsFile(statsFileName.c_str());
	//statsFile.precision(5);
	statsFile<<"ITER"<<"\t"<<"SIM_LL"<<"\t";
	int p;
	for(p = 0; p < paramVec.size(); ++p) 
		statsFile<<pSimModel->getParamName(paramVec[p]._type)<<"\t";
	statsFile<<"SIM_ROOT"<<"\t"<<"POST_SROOT"<<"\t"<<"SIM_INTERNALS"<<"\t"<<"POST_SINTERNALS"<<"\t";
	statsFile<<"OPT_LL"<<"\t"<<"DIFF_LL"<<"\t";
	for(p = 0; p < paramVec.size(); ++p) 
		statsFile<<pSimModel->getParamName(paramVec[p]._type)<<"\t";
	statsFile<<"OPT_ROOT"<<"\t"<<"POST_ROOT"<<"\t"<<"OPT_INTERNALS"<<"\t"<<"POST_INTERNALS"<<endl;
	//expectation stats file
	string statsExpFileName = baseOutDirName + "//" + chrNumberOptions::_simulationsJumpsStats;
	ofstream statsExpFile(statsExpFileName.c_str());
	statsExpFile<<"ITER"<<"\t"<<"sim_gain"<<"\t"<<"sim_loss"<<"\t"<<"sim_dupl"<<"\t"<<"sim_demi"<<"\t"<<"sim_maxChr"<<"\t";
	statsExpFile<<"opt_gain"<<"\t"<<"opt_loss"<<"\t"<<"opt_dupl"<<"\t"<<"opt_demi"<<"\t";
	statsExpFile<<"diff_gain"<<"\t"<<"diff_loss"<<"\t"<<"diff_dupl"<<"\t"<<"diff_demi"<<endl;
	int expIter = -1;//number of iterations for which expectations were computed
	int simGains=-1, simLoss=-1, simDupl=-1, simDemi=-1, simBaseNum=-1, simMaxChr=-1;  //save the number of events for that simulation - to be stored later
	for (int iter = 0; iter < iterNum; ++iter)
	{
		int treeIter = iter + chrNumberOptions::_startSimulationsIter;
		LOGnOUT(3, <<"simulation iter = "<<treeIter<<endl);
		//3. open simulated directory and out files
		string iterOutDirName = baseOutDirName + "//" + int2string(treeIter);
		chrNumberOptions::_outDir = iterOutDirName;
		createDir("", iterOutDirName);
		string simTreeFileName = chrNumberOptions::_outDir + "//" + "simTree.phr";
		ofstream treeStream(simTreeFileName.c_str());
		string simScFileName = chrNumberOptions::_outDir + "//" + "simCounts.txt";
		ofstream simMsa(simScFileName.c_str());
		string simScAllNodesFileName = chrNumberOptions::_outDir + "//" + "simCountsAllNodes.txt";
		ofstream simMsaAllNodes(simScAllNodesFileName.c_str());
		string trueModelMLTreeFileName = chrNumberOptions::_outDir + "//" + "simModelMlTree.phr";
		string trueModelPostTreeFileName = chrNumberOptions::_outDir + "//" + "simModelPostTree.phr";
		string optModelMLTreeFileName = chrNumberOptions::_outDir + "//" + "optModelMlTree.phr";
		string optModelPostTreeFileName = chrNumberOptions::_outDir + "//" + "optModelPostTree.phr";
		
		//4. run simulation
		getSimulatedTree(treeIter); 
		sequenceContainer allNodesSimSc;
		sequenceContainer leavesSimSc;

		//decide which simulation to perform
		if (bJumpSimulation == false)
		{
			simulateTree sim(_tree, *pSimSp, &simAlph);
			sim.generate_seq(1);
			allNodesSimSc = sim.toSeqData();
			leavesSimSc = sim.toSeqDataWithoutInternalNodes();
		}
		else
		{
			simulateChangesAlongTree simJumps(_tree, *pSimSp, &simAlph);
			allNodesSimSc = simJumps.simulatePosition();
			leavesSimSc = simJumps.toSeqDataWithoutInternalNodes();
            simGains = simJumps.getTotalJumps(chrNumModel::GAIN_J);
			simLoss = simJumps.getTotalJumps(chrNumModel::LOSS_J);
			simDupl = simJumps.getTotalJumps(chrNumModel::DUPL_J);
			simDemi = simJumps.getTotalJumps(chrNumModel::DEMI_J);
			simBaseNum = simJumps.getTotalJumps(chrNumModel::BASE_J);
			simMaxChr = simJumps.getTotalJumps(chrNumModel::MAX_CHR_J);
			string eventsFileName = iterOutDirName + "//simEvents.txt";
            ofstream eventsFile(eventsFileName.c_str());
			simJumps.printEvents(eventsFile);
            eventsFile.close();
		}
		//get simulated tree node states
		Vint simStates(_tree.getNodesNum(),-2);
		treeIterTopDownConst tIt(_tree);
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
			string nodename = mynode->name();
			int seqId = allNodesSimSc.getId(nodename);
			simStates[mynode->id()] = allNodesSimSc[seqId][0];
		}
		//print simulated tree
		printTreeStatesAsBPValues(treeStream, simStates, _tree.getRoot(), NULL, false, &simAlph);
		treeStream.close();
		//print simulated MSA
		fastaFormat::write(simMsaAllNodes, allNodesSimSc);
		simMsaAllNodes.close();
		fastaFormat::write(simMsa, leavesSimSc);
		simMsa.flush();
		simMsa.close();
		chrNumberOptions::_dataFile = simScFileName;
		///end of simulations //////////////////////
		if (chrNumberOptions::_maxOptimizationIterations < 0)
			continue;
		////start of inference
		getStartingData();
		pair<int,int> chrRange = getChrRangeInData(_sc);
		initStartingModel(_pAlph, chrRange.second - chrRange.first);
		simLLVec[iter] = compLogLikelihood(_pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
		//5. infer ancestrals with the simulated model
		chrNumberOptions::_mainType = "Run_Fix_Param";
		VVVdouble simjointPost;
		suffStatGlobalHomPos ssc;
		computeJointPosterior(simjointPost, ssc);
		inferAncestralStates(simjointPost); //put the posterior and ML reconstruction in _ancestralPost[][] and _ancestralMl[]
		printReconstructedTrees(trueModelMLTreeFileName, trueModelPostTreeFileName);
		//store results of simModel
		getSimulationResults(iter, simStates, simParams, simMlRoot, simPostRoot, simMlAllInternals, simPostAllInternals);

		//6. run optimization
		chrNumberOptions::_mainType = "Optimize_Model";
		_bestLL = optLLVec[iter] = optimizeParams();
		LOG(5, <<"opt LL = "<<optLLVec[iter]<<endl;);
		VVVdouble optjointPost;
		suffStatGlobalHomPos ssc1;
		computeJointPosterior(optjointPost, ssc1);
		inferAncestralStates(optjointPost); //put the posterior and ML reconstruction in _ancestralPost[][] and _ancestralMl[]
		if (bJumpSimulation)
		{
			computeChangeProbAndExp(optjointPost, ssc1);
			if (performSimulationsForExpectation())
			{
				++expIter;
				simTotalJumps[chrNumModel::GAIN_J].push_back(simGains);
				simTotalJumps[chrNumModel::LOSS_J].push_back(simLoss);
				simTotalJumps[chrNumModel::DUPL_J].push_back(simDupl);
				simTotalJumps[chrNumModel::DEMI_J].push_back(simDemi);
				simTotalJumps[chrNumModel::BASE_J].push_back(simBaseNum);
				simTotalJumps[chrNumModel::MAX_CHR_J].push_back(simMaxChr);
                optTotalJumps[chrNumModel::GAIN_J].push_back(_expChanges[_expChanges.size()-1][chrNumModel::GAIN_J]);
				optTotalJumps[chrNumModel::LOSS_J].push_back(_expChanges[_expChanges.size()-1][chrNumModel::LOSS_J]);
				optTotalJumps[chrNumModel::DUPL_J].push_back(_expChanges[_expChanges.size()-1][chrNumModel::DUPL_J]);
				optTotalJumps[chrNumModel::DEMI_J].push_back(_expChanges[_expChanges.size()-1][chrNumModel::DEMI_J]);
				optTotalJumps[chrNumModel::BASE_J].push_back(_expChanges[_expChanges.size()-1][chrNumModel::BASE_J]);
				diffTotalJumps[chrNumModel::GAIN_J].push_back(fabs(optTotalJumps[chrNumModel::GAIN_J][expIter] - simTotalJumps[chrNumModel::GAIN_J][expIter]));
				diffTotalJumps[chrNumModel::LOSS_J].push_back(fabs(optTotalJumps[chrNumModel::LOSS_J][expIter] - simTotalJumps[chrNumModel::LOSS_J][expIter]));
				diffTotalJumps[chrNumModel::DUPL_J].push_back(fabs(optTotalJumps[chrNumModel::DUPL_J][expIter] - simTotalJumps[chrNumModel::DUPL_J][expIter]));
				diffTotalJumps[chrNumModel::DEMI_J].push_back(fabs(optTotalJumps[chrNumModel::DEMI_J][expIter] - simTotalJumps[chrNumModel::DEMI_J][expIter]));
				diffTotalJumps[chrNumModel::BASE_J].push_back(fabs(optTotalJumps[chrNumModel::BASE_J][expIter] - simTotalJumps[chrNumModel::BASE_J][expIter]));
			}
		}
		//printResults();
		errorMsg::reportError("If I got here than the version is not updated - printResults is not executed");
		//store results of optModel
		getSimulationResults(iter, simStates, optParams, optMlRoot, optPostRoot, optMlAllInternals, optPostAllInternals);

		statsFile<<treeIter<<"\t"<<simLLVec[iter]<<"\t";
		for(p = 0; p < simParams.size(); ++p) 
			statsFile<<simParams[p][iter]<<"\t";
		statsFile<<simMlRoot[iter]<<"\t"<<simPostRoot[iter]<<"\t"<<simMlAllInternals[iter]<<"\t"<<simPostAllInternals[iter]<<"\t";
		statsFile<<optLLVec[iter]<<"\t"<<optLLVec[iter]-simLLVec[iter]<<"\t";
		for(p = 0; p < optParams.size(); ++p) 
			statsFile<<optParams[p][iter]<<"\t";
		statsFile<<optMlRoot[iter]<<"\t"<<optPostRoot[iter]<<"\t"<<optMlAllInternals[iter]<<"\t"<<optPostAllInternals[iter]<<endl;
		if (bJumpSimulation)
		{
			if (performSimulationsForExpectation())
			{
				statsExpFile<<treeIter<<"\t"<<simTotalJumps[chrNumModel::GAIN_J][expIter]<<"\t"<<simTotalJumps[chrNumModel::LOSS_J][expIter]<<"\t"<<simTotalJumps[chrNumModel::DUPL_J][expIter]<<"\t"<<simTotalJumps[chrNumModel::DEMI_J][expIter]<<"\t"<<simTotalJumps[chrNumModel::MAX_CHR_J][expIter]<<"\t";			
				statsExpFile<<optTotalJumps[chrNumModel::GAIN_J][expIter]<<"\t"<<optTotalJumps[chrNumModel::LOSS_J][expIter]<<"\t"<<optTotalJumps[chrNumModel::DUPL_J][expIter]<<"\t"<<optTotalJumps[chrNumModel::DEMI_J][expIter]<<"\t";
				statsExpFile<<diffTotalJumps[chrNumModel::GAIN_J][expIter]<<"\t"<<diffTotalJumps[chrNumModel::LOSS_J][expIter]<<"\t"<<diffTotalJumps[chrNumModel::DUPL_J][expIter]<<"\t"<<diffTotalJumps[chrNumModel::DEMI_J][expIter]<<endl;
			}
		}
	}
	MDOUBLE avgLLsim = computeAverage(simLLVec);
	MDOUBLE avgLLopt = computeAverage(optLLVec);

	statsFile<<"AVG"<<"\t"<<avgLLsim<<"\t";
	for(p = 0; p < simParams.size(); ++p) 
		statsFile<<computeAverage(simParams[p])<<"\t";
	statsFile<<computeAverage(simMlRoot)<<"\t"<<computeAverage(simPostRoot)<<"\t"<<computeAverage(simMlAllInternals)<<"\t"<<computeAverage(simPostAllInternals)<<"\t";
	statsFile<<avgLLopt<<"\t"<<avgLLopt - avgLLsim<<"\t";
	for(p = 0; p < optParams.size(); ++p) 
		statsFile<<computeAverage(optParams[p])<<"\t";
	statsFile<<computeAverage(optMlRoot)<<"\t"<<computeAverage(optPostRoot)<<"\t"<<computeAverage(optMlAllInternals)<<"\t"<<computeAverage(optPostAllInternals)<<endl;

	statsFile<<"STD"<<"\t"<<computeStd(simLLVec)<<"\t";
	for(p = 0; p < simParams.size(); ++p) 
		statsFile<<computeStd(simParams[p])<<"\t";
	statsFile<<computeStd(simMlRoot)<<"\t"<<computeStd(simPostRoot)<<"\t"<<computeStd(simMlAllInternals)<<"\t"<<computeStd(simPostAllInternals)<<"\t";
	statsFile<<computeStd(optLLVec)<<"\t"<<"\t";
	for(p = 0; p < optParams.size(); ++p) 
		statsFile<<computeStd(optParams[p])<<"\t";
	statsFile<<computeStd(optMlRoot)<<"\t"<<computeStd(optPostRoot)<<"\t"<<computeStd(optMlAllInternals)<<"\t"<<computeStd(optPostAllInternals)<<endl;
	if (bJumpSimulation)
	{
		statsExpFile<<"AVG"<<"\t"<<computeAverage(simTotalJumps[chrNumModel::GAIN_J])<<"\t"<<computeAverage(simTotalJumps[chrNumModel::LOSS_J])<<"\t"<<computeAverage(simTotalJumps[chrNumModel::DUPL_J])<<"\t"<<computeAverage(simTotalJumps[chrNumModel::DEMI_J])<<"\t"<<computeAverage(simTotalJumps[chrNumModel::MAX_CHR_J])<<"\t";			
		statsExpFile<<computeAverage(optTotalJumps[chrNumModel::GAIN_J])<<"\t"<<computeAverage(optTotalJumps[chrNumModel::LOSS_J])<<"\t"<<computeAverage(optTotalJumps[chrNumModel::DUPL_J])<<"\t"<<computeAverage(optTotalJumps[chrNumModel::DEMI_J])<<"\t";
		statsExpFile<<computeAverage(diffTotalJumps[chrNumModel::GAIN_J])<<"\t"<<computeAverage(diffTotalJumps[chrNumModel::LOSS_J])<<"\t"<<computeAverage(diffTotalJumps[chrNumModel::DUPL_J])<<"\t"<<computeAverage(diffTotalJumps[chrNumModel::DEMI_J])<<endl;
		statsExpFile<<"STD"<<"\t"<<computeStd(simTotalJumps[chrNumModel::GAIN_J])<<"\t"<<computeStd(simTotalJumps[chrNumModel::LOSS_J])<<"\t"<<computeStd(simTotalJumps[chrNumModel::DUPL_J])<<"\t"<<computeStd(simTotalJumps[chrNumModel::DEMI_J])<<"\t"<<computeStd(simTotalJumps[chrNumModel::MAX_CHR_J])<<"\t";			
		statsExpFile<<computeStd(optTotalJumps[chrNumModel::GAIN_J])<<"\t"<<computeStd(optTotalJumps[chrNumModel::LOSS_J])<<"\t"<<computeStd(optTotalJumps[chrNumModel::DUPL_J])<<"\t"<<computeStd(optTotalJumps[chrNumModel::DEMI_J])<<"\t";
		statsExpFile<<computeStd(diffTotalJumps[chrNumModel::GAIN_J])<<"\t"<<computeStd(diffTotalJumps[chrNumModel::LOSS_J])<<"\t"<<computeStd(diffTotalJumps[chrNumModel::DUPL_J])<<"\t"<<computeStd(diffTotalJumps[chrNumModel::DEMI_J])<<endl;
	}
	
	statsFile.close();
	statsExpFile.close();
	delete pSimSp;
}



void chrNumberMng::runSimulations()
{
	printRunInfo();
	int iterNum = chrNumberOptions::_simulationsIter - chrNumberOptions::_startSimulationsIter;
	bool bJumpSimulation = !chrNumberOptions::_simulationsJumpsStats.empty();
	if (bJumpSimulation && (chrNumberOptions::_simulationsNum < 2))
		errorMsg::reportError("performing jumpSimulations but cannot make inference");
	if (!bJumpSimulation)
		chrNumberOptions::_simulationsNum = 0;
	//1a. init model with simulated params - this is actually just for printing
	finiteIntAlphabet simAlph(1,chrNumberOptions::_maxChrNumForSimulations);
	initStartingModel(&simAlph, chrNumberOptions::_maxBaseTransition);
	stochasticProcess* pSimSp = _pSp->clone();
	delete _pSp;
	_pSp = NULL;
	chrNumModel* pSimModel = static_cast<chrNumModel*>(pSimSp->getPijAccelerator()->getReplacementModel());
	vector<paramObj> paramVec;
	pSimModel->getParams(paramVec);
		
	//2. open global output file and print headings 
	string baseOutDirName = chrNumberOptions::_outDir;
	//statsFile.precision(5);
	int simGains=-1, simLoss=-1, simDupl=-1, simDemi=-1, simBaseNum=-1, simMaxChr=-1;  //save the number of events for that simulation - to be stored later
	for (int iter = 0; iter < iterNum; ++iter)
	{
		int treeIter = iter + chrNumberOptions::_startSimulationsIter;
		LOGnOUT(3, <<"simulation iter = "<<treeIter<<endl);
		//3. open simulated directory and out files
		string iterOutDirName = baseOutDirName + "//" + int2string(treeIter);
		chrNumberOptions::_outDir = iterOutDirName;
		createDir("", iterOutDirName);
		string simTreeFileName = chrNumberOptions::_outDir + "//" + "simTree.phr";
		ofstream treeStream(simTreeFileName.c_str());
		string simScFileName = chrNumberOptions::_outDir + "//" + "simCounts.txt";
		ofstream simMsa(simScFileName.c_str());
		string simScAllNodesFileName = chrNumberOptions::_outDir + "//" + "simCountsAllNodes.txt";
		ofstream simMsaAllNodes(simScAllNodesFileName.c_str());

		//4. run simulation
		getSimulatedTree(treeIter); 
		sequenceContainer allNodesSimSc;
		sequenceContainer leavesSimSc;

		//decide which simulation to perform
		if (bJumpSimulation == false)
		{
			simulateTree sim(_tree, *pSimSp, &simAlph);
			sim.generate_seq(1);
			allNodesSimSc = sim.toSeqData();
			leavesSimSc = sim.toSeqDataWithoutInternalNodes();
		}
		else
		{
			simulateChangesAlongTree simJumps(_tree, *pSimSp, &simAlph);
			allNodesSimSc = simJumps.simulatePosition();
			leavesSimSc = simJumps.toSeqDataWithoutInternalNodes();
            simGains = simJumps.getTotalJumps(chrNumModel::GAIN_J);
			simLoss = simJumps.getTotalJumps(chrNumModel::LOSS_J);
			simDupl = simJumps.getTotalJumps(chrNumModel::DUPL_J);
			simDemi = simJumps.getTotalJumps(chrNumModel::DEMI_J);
			simBaseNum = simJumps.getTotalJumps(chrNumModel::BASE_J);
			simMaxChr = simJumps.getTotalJumps(chrNumModel::MAX_CHR_J);
			string eventsFileName = iterOutDirName + "//simEvents.txt";
            ofstream eventsFile(eventsFileName.c_str());
			simJumps.printEvents(eventsFile);
            eventsFile.close();
		}
		//get simulated tree node states
		Vint simStates(_tree.getNodesNum(),-2);
		treeIterTopDownConst tIt(_tree);
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
			string nodename = mynode->name();
			int seqId = allNodesSimSc.getId(nodename);
			simStates[mynode->id()] = allNodesSimSc[seqId][0];
		}
		//print simulated tree
		printTreeStatesAsBPValues(treeStream, simStates, _tree.getRoot(), NULL, false, &simAlph);
		treeStream.close();
		//print simulated MSA
		fastaFormat::write(simMsaAllNodes, allNodesSimSc);
		simMsaAllNodes.close();
		fastaFormat::write(simMsa, leavesSimSc);
		simMsa.flush();
		simMsa.close();
		
	}
	delete pSimSp;
}



//store results of simulations
//mlRoot: the error in ML reconstruction: the abs difference between the simulated and inferred root state
//mlAllInternals: same as above for all internal nodes 
//postRoot: the posterior of the simulated root state
//postAllInternals: the average posterior for all internal nodes for the simulated states 
void chrNumberMng::getSimulationResults(int iterNum, Vint& simStates, VVdouble& params, Vint& mlRoot, Vdouble& postRoot, Vdouble& mlAllInternals, Vdouble& postAllInternals)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel());
	vector<paramObj> paramVec;
	pModel->getParams(paramVec);
	if (params.size() != paramVec.size())
		errorMsg::reportError("parameters size is not the same in chrNumberMng::getSimulationResults");
	for(int p = 0; p < paramVec.size(); ++p) 
		params[p][iterNum] = paramVec[p]._value;

	int rootId = _tree.getRoot()->id();
	int simRootState = simStates[rootId];
	int mlRootState = _ancestralMl[rootId];
	MDOUBLE postSimRoot = _ancestralPost[rootId][simRootState];
	mlRoot[iterNum] = abs(mlRootState - simRootState);
	postRoot[iterNum] = postSimRoot;

	//loop over all internal nodes
	vector<tree::nodeP> internalsVec;
    _tree.getAllHTUs(internalsVec, _tree.getRoot());
	MDOUBLE sumML = 0.0, sumPost = 0.0;
	for (int n=0; n< internalsVec.size(); ++n) 
	{
		int nodeId = internalsVec[n]->id();
		int simState = simStates[nodeId];
		int mlState = _ancestralMl[nodeId];
		sumML += abs(mlState - simState);
		MDOUBLE postSim = _ancestralPost[nodeId][simState];
		sumPost += postSim;
	}
	mlAllInternals[iterNum] = sumML / static_cast<MDOUBLE>(internalsVec.size());
	postAllInternals[iterNum] = sumPost / static_cast<MDOUBLE>(internalsVec.size());
}


void chrNumberMng::runOneSimulation()
{
    printRunInfo();
	getTree();
	if (_pAlph)
		delete _pAlph;
	_pAlph = new finiteIntAlphabet(1,chrNumberOptions::_maxChrNumForSimulations);
	initStartingModel(_pAlph, chrNumberOptions::_maxBaseTransition);
	simulateTree sim(_tree, *_pSp, _pAlph);
	sim.generate_seq(1);
	sequenceContainer allNodesSimSc = sim.toSeqData();
	_sc = sim.toSeqDataWithoutInternalNodes();
	
	//print simulated tree
	Vint states(_tree.getNodesNum(),-2);
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
//		if (mynode->isRoot()) continue;
		string nodename = mynode->name();
		int seqId = allNodesSimSc.getId(nodename);
		states[mynode->id()] = allNodesSimSc[seqId][0];
	}

	string simFileName = chrNumberOptions::_outDir + "//" + "simTree.phr";
	ofstream oStream(simFileName.c_str());
	printTreeStatesAsBPValues(oStream, states, _tree.getRoot(), NULL, false, _pAlph);
	oStream.close();

	string simScFileName = chrNumberOptions::_outDir + "//" + "simCounts.txt";
	ofstream simMsa(simScFileName.c_str());
	fastaFormat::write(simMsa, _sc);
	simMsa.close();
	simScFileName = chrNumberOptions::_outDir + "//" + "simCountsAllNodes.txt";
	ofstream simMsaAllNodes(simScFileName.c_str());
	fastaFormat::write(simMsaAllNodes, allNodesSimSc);
	simMsaAllNodes.close();

	_bestLL = compLogLikelihood(_pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
	cerr<<"sim LL = "<<_bestLL<<endl;
	///end of simulations //////////////////////
	
	//reconstruct ML states given the simulated model or an optimized model
	bool bOptimize = false;
	if (bOptimize)
	{
		chrNumberOptions::_mainType = "Optimize_Model";
		MDOUBLE ll = optimizeParams();
		cerr<<"ml LL = "<<ll <<endl;
	}
	
	VVVdouble jointPost;
	suffStatGlobalHomPos ssc;
	computeJointPosterior(jointPost, ssc);
	inferAncestralStates(jointPost);
	computeChangeProbAndExp(jointPost, ssc);
	printResults(jointPost);
}

//infer ancestral states using both ML and posterior probabilities and stores in members _ancestralPost and _ancestralMl.
//returns the joint posterior probabilities in jointPost
void chrNumberMng::inferAncestralStates(const VVVdouble& jointPost)
{
	ancestralReconstructML(_ancestralMl);
	//string ancestralTreeStr = chrNumberOptions::_outDir + "//" + chrNumberOptions::_inferTreeFile;
	//ofstream treeStream(ancestralTreeStr.c_str());
	//printTreeStatesAsBPValues(treeStream, _ancestralMl, _tree.getRoot(), NULL, false, _pAlph);
	//treeStream.close();

	computeAncestralPosterior(jointPost);
	string posteriorFileName = chrNumberOptions::_outDir + "//ancestorsProbs.txt";
	ofstream posteriorFile(posteriorFileName.c_str());
	printPosterior(posteriorFile, _ancestralPost);
	posteriorFile.close();
}

//print the ML and posterior reconstruction
void chrNumberMng::printReconstructedTrees(const string& mlTree, const string& postTree)
{
	Vstring mlPrintData(_tree.getNodesNum());
	Vstring postPrintData(_tree.getNodesNum());
	int alphaSize = _pAlph->size();		
	for (int i=0; i< mlPrintData.size(); ++i) 
	{
		tree::nodeP myNode = _tree.findNodeById(i);
		if (!myNode)
			errorMsg::reportError("error in chrNumberMng::printTrees, cannot find node");
		if (myNode->isLeaf())
			continue;

		mlPrintData[i] = "[" + myNode->name() + "-" + _pAlph->fromInt(_ancestralMl[myNode->id()]) + "]"; //gainX + "//" + gainA + "//" + lossX + "//" + lossA + "[" + myNode->name() + "]";
		//order _ancestralPost so to print only the three highest states
		vector< vecElem<MDOUBLE> > orders;
		orderVec(_ancestralPost[myNode->id()], orders);
		postPrintData[i] = "[" + myNode->name() + "_" + _pAlph->fromInt(orders[alphaSize-1].getPlace()) + "-" + double2string(orders[alphaSize-1].getValue(),2);
		postPrintData[i] += "//" + _pAlph->fromInt(orders[alphaSize-2].getPlace()) + "-" + double2string(orders[alphaSize-2].getValue(),2);
		postPrintData[i] += "]"; 
	}

	tree printTree(_tree);
	changesNamesToLeaves(printTree, _sc);
    ofstream treeStream(mlTree.c_str());
	printDataOnTreeAsBPValues(treeStream, mlPrintData, printTree);
	treeStream.close();

	ofstream postTreeStream(postTree.c_str());
	printDataOnTreeAsBPValues(postTreeStream, postPrintData, printTree);
	postTreeStream.close();
}


void chrNumberMng::printPosterior(ostream &out, const VVdouble& postProbs)
{
	//print heading
	out<<"NODE"<<"\t";
	int letter = _pAlph->min();
	int alphabetSize = _pSp->alphabetSize();
	for (letter = 0; letter < alphabetSize; ++letter)
		out<<_pAlph->id2Count(letter)<<"\t";
	out<<endl;
	
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		out<<mynode->name()<<"\t";
		for (letter = 0; letter < alphabetSize; ++letter)
			out<<postProbs[mynode->id()][letter]<<"\t";
		out<<endl;
	}
}


//should change here - go over jointPost and see the nodes that the node most probably changed
void chrNumberMng::printExpectations(const VVVdouble& jointPost, ostream &expFile)
{
	int maxChr = _pSp->alphabetSize();
	bool bEstimateMissingEvents = true;
	treeIterTopDownConst tItr(_tree);
	if (!performSimulationsForExpectation())
	{//heuristic computation of transitions 
		expFile<<"#cannot perform simulations: Expectations are based on heuristic computations rather than simulations"<<endl;
		for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
		{
			if (mynode->isRoot())
				continue;
			int nodeId = mynode->id();
			int fatherId = mynode->father()->id();
			if (_ancestralMl[nodeId] != _ancestralMl[fatherId])
			{
				if ((_ancestralMl[nodeId] == _pAlph->unknown()) || (_ancestralMl[fatherId] == _pAlph->unknown()))
					continue;
				MDOUBLE sumExp = 0.0;
				for (int e = 0; e < _expChanges[nodeId].size(); ++e)
					sumExp += _expChanges[nodeId][e];
				Vdouble est = estimateExpectationsHeuristicly(nodeId, jointPost);
				_expChanges[nodeId][chrNumModel::GAIN_J] = est[chrNumModel::GAIN_J];
				_expChanges[nodeId][chrNumModel::LOSS_J] = est[chrNumModel::LOSS_J];
				_expChanges[nodeId][chrNumModel::DUPL_J] = est[chrNumModel::DUPL_J];
				_expChanges[nodeId][chrNumModel::DEMI_J] = est[chrNumModel::DEMI_J];
				_expChanges[nodeId][chrNumModel::BASE_J] = est[chrNumModel::BASE_J];
			}
		}
	}
	
	//print gains
	expFile<<"#Nodes with GAIN events with expectation above "<<EXP_THRESHOLD<<endl;
	MDOUBLE sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _expChanges[nodeId][chrNumModel::GAIN_J];
		if (_expChanges[nodeId][chrNumModel::GAIN_J] > EXP_THRESHOLD)
			expFile<<mynode->name()<<": "<<_expChanges[nodeId][chrNumModel::GAIN_J]<<endl;//"   prob = "<<_probChanges[nodeId][chrNumModel::GAIN_J]<<endl;
	}
	expFile<<"Total number of gain events: "<<sum<<endl;

	//print losses
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with LOSS events with expectation above "<<EXP_THRESHOLD<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _expChanges[nodeId][chrNumModel::LOSS_J];
		if (_expChanges[nodeId][chrNumModel::LOSS_J] > EXP_THRESHOLD)
			expFile<<mynode->name()<<": "<<_expChanges[nodeId][chrNumModel::LOSS_J]<<endl;//"   prob = "<<_probChanges[nodeId][chrNumModel::LOSS_J]<<endl;
	}
	expFile<<"Total number of loss events: "<<sum<<endl;

	//print duplications
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with duplication events with expectation above "<<EXP_THRESHOLD<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _expChanges[nodeId][chrNumModel::DUPL_J];
		if (_expChanges[nodeId][chrNumModel::DUPL_J] > EXP_THRESHOLD)
			expFile<<mynode->name()<<": "<<_expChanges[nodeId][chrNumModel::DUPL_J]<<endl;//"   prob = "<<_probChanges[nodeId][chrNumModel::DUPL_J]<<endl;
	}
	expFile<<"Total number of duplication events: "<<sum<<endl;

	//print demi-duplications
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	expFile<<"#Nodes with demi-duplication events with expectation above "<<EXP_THRESHOLD<<endl;
	sum = 0.0;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		sum += _expChanges[nodeId][chrNumModel::DEMI_J];
		if (_expChanges[nodeId][chrNumModel::DEMI_J] > EXP_THRESHOLD)
			expFile<<mynode->name()<<": "<<_expChanges[nodeId][chrNumModel::DEMI_J]<<endl;//"   prob = "<<_probChanges[nodeId][chrNumModel::DEMI_J]<<endl;
	}
	expFile<<"Total number of demi-dulications events: "<<sum<<endl;
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	//print base-number transitions
	bool bBaseNumber = (IGNORE_PARAM != static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel())->getBaseNumberR());
	if (bBaseNumber)
	{
		expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
		expFile<<"#Nodes with transitions in base number with expectation above "<<EXP_THRESHOLD<<endl;
		sum = 0.0;
		for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
		{
			if (mynode->isRoot())
				continue;
			int nodeId = mynode->id();
			sum += _expChanges[nodeId][chrNumModel::BASE_J];
			if (_expChanges[nodeId][chrNumModel::BASE_J] > EXP_THRESHOLD)
				expFile<<mynode->name()<<": "<<_expChanges[nodeId][chrNumModel::BASE_J]<<endl;
		}
		expFile<<"Total number of transitions in base number: "<<sum<<endl;
		expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	}	
	//print events not accounted for in the simulations defined as:
	//the ML reconstruction of node != ML reconstruct of father(node) AND the total expectation over all event types is smaller than EXP_THRESHOLD
	expFile<<"#EVENTS NOT ACCOUNTED FOR IN THE SIMULATIONS: "<<endl;
	sum = 0.0;
	bool bMissing = false;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		int fatherId = mynode->father()->id();
		if ((_ancestralMl[nodeId] != _ancestralMl[fatherId]) && (getHeuristicProbOfChange(nodeId, jointPost) > EXP_THRESHOLD2))
		{
			if ((_ancestralMl[nodeId] == _pAlph->unknown()) || (_ancestralMl[fatherId] == _pAlph->unknown()))
				continue;
			bool bComposite = false;
			if ((_pAlph->isComposite(_ancestralMl[nodeId])) || (_pAlph->isComposite(_ancestralMl[fatherId])))
				bComposite = true;
			MDOUBLE sumExp = 0.0;
			for (int e = 0; e < _expChanges[nodeId].size(); ++e)
				sumExp += _expChanges[nodeId][e];
			if ((sumExp < EXP_THRESHOLD) || (bComposite) || (_expChanges[nodeId][chrNumModel::MAX_CHR_J] > EXP_THRESHOLD_MAX_CHR))
			{
				bMissing = true;
				expFile<<mynode->father()->name()<<":ML="<<_pAlph->fromInt(_ancestralMl[fatherId])<<" to "<<mynode->name()<<":ML="<<_pAlph->fromInt(_ancestralMl[nodeId])<<"  ProbOfChange = "<<getHeuristicProbOfChange(nodeId, jointPost)<<"  sum expectations= "<<sumExp<<endl;
				if (bEstimateMissingEvents )
				{
					Vdouble est = estimateExpectationsHeuristicly(nodeId, jointPost);
					_expChanges[nodeId][chrNumModel::GAIN_J] = est[chrNumModel::GAIN_J];
					_expChanges[nodeId][chrNumModel::LOSS_J] = est[chrNumModel::LOSS_J];
					_expChanges[nodeId][chrNumModel::DUPL_J] = est[chrNumModel::DUPL_J];
					_expChanges[nodeId][chrNumModel::DEMI_J] = est[chrNumModel::DEMI_J];
					_expChanges[nodeId][chrNumModel::BASE_J] = est[chrNumModel::BASE_J];
					if (bBaseNumber)
						expFile<<"\t"<<"HEURISTIC ESTIMATION:\t"<<"Gains = "<<est[chrNumModel::GAIN_J]<<"\t"<<"Losses = "<<est[chrNumModel::LOSS_J]<<"\t"<<"Duplications = "<<est[chrNumModel::DUPL_J]<<"\t"<<"Demi-dupl = "<<est[chrNumModel::DEMI_J]<<"\t"<<"Base-number = "<<est[chrNumModel::BASE_J]<<endl;
					else
						expFile<<"\t"<<"HEURISTIC ESTIMATION:\t"<<"Gains = "<<est[chrNumModel::GAIN_J]<<"\t"<<"Losses = "<<est[chrNumModel::LOSS_J]<<"\t"<<"Duplications = "<<est[chrNumModel::DUPL_J]<<"\t"<<"Demi-dupl = "<<est[chrNumModel::DEMI_J]<<endl;
					}
			}
		}
	}
	if (bMissing == false)
        expFile<<"NONE"<<endl;
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;


	//print all expectaiotns for all nodes 
	expFile<<"#ALL EVENTS EXPECTATIONS PER NODE"<<endl;
	expFile<<"NODE"<<"\t";
	if (bBaseNumber)
		expFile<<"GAIN"<<"\t"<<"LOSS"<<"\t"<<"DUPLICATION"<<"\t"<<"DEMI-DUPLICATION"<<"\t"<<"BASE-NUMBER"<<endl;
	else
		expFile<<"GAIN"<<"\t"<<"LOSS"<<"\t"<<"DUPLICATION"<<"\t"<<"DEMI-DUPLICATION"<<endl;
	for (tree::nodeP mynode = tItr.first(); mynode != tItr.end(); mynode = tItr.next()) 
	{
		if (mynode->isRoot())
			continue;
		expFile<<mynode->name()<<"\t";
		int nodeId = mynode->id();
		if (bBaseNumber)
			expFile<<_expChanges[nodeId][chrNumModel::GAIN_J]<<"\t"<<_expChanges[nodeId][chrNumModel::LOSS_J]<<"\t"<<_expChanges[nodeId][chrNumModel::DUPL_J]<<"\t"<<_expChanges[nodeId][chrNumModel::DEMI_J]<<"\t"<<_expChanges[nodeId][chrNumModel::BASE_J]<<endl;
		else
			expFile<<_expChanges[nodeId][chrNumModel::GAIN_J]<<"\t"<<_expChanges[nodeId][chrNumModel::LOSS_J]<<"\t"<<_expChanges[nodeId][chrNumModel::DUPL_J]<<"\t"<<_expChanges[nodeId][chrNumModel::DEMI_J]<<endl;
	}
	expFile<<"#+++++++++++++++++++++++++++++"<<endl<<endl;
	//for each leaf print the number of transitions (for each event type) from root to leaf
	VVdouble nodesTransitions; //leavesTransitions[leafId][chrNumModel::jumpType]
	computeExpectationsFromRootToNodes(nodesTransitions);
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
	}
}

//recieves the precalculated joint post and return an estimate of the probability that a change has occured along the branch leading to nodeID.
MDOUBLE chrNumberMng::getHeuristicProbOfChange(int nodeID, const VVVdouble& jointPost)
{
	int alphabetSize = _pSp->alphabetSize();
	MDOUBLE sumIdenticalTerminals = 0.0;
	for (int i = 0; i < alphabetSize; ++i){
		sumIdenticalTerminals += jointPost[nodeID][i][i];
	}
	if (sumIdenticalTerminals > 1.0001) {
		string errM = "error in getHeuristicProbOfChange. sumIdenticalTerminals = " + double2string(sumIdenticalTerminals);
		errorMsg::reportError(errM);
	}
	return 1.0-sumIdenticalTerminals;
}
//compute the posterior probabilitie for all states for each ancestral node. 
//The model is irreversible so down computations are conditional on the state at the root.
//--> have to sum also over root states
//prob(state|Data)=sum1{root states}sum2{fatherLetter}[P(root=rootState)*Down(fatherState|rootState)*P(fatherState->state)up(state)]}}
//In practice - we can also use the pre-calculated joint probability P(N=x, father(N)=y|D) and just sum over these probs:
//Prpb(State|Data) = sum{fatherState}[P(N=x, father(N)=y|D)]}
//return the joint posterior probabilities in jointPost
void chrNumberMng::computeAncestralPosterior(const VVVdouble& jointPost)
{
	int numNodes = _tree.getNodesNum();
	int alphabetSize = _pSp->alphabetSize();
	resizeMatrix(_ancestralPost, numNodes, alphabetSize);
    
	treeIterTopDownConst tIt(_tree);
	int letter;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next())
	{
        if (mynode->isRoot()) 
		{
			for(letter = 0; letter<alphabetSize; ++letter) 
                _ancestralPost[mynode->id()][letter] = jointPost[mynode->id()][0][letter];
			continue;
		}
		for(letter = 0; letter < alphabetSize; ++letter) 
		{
			MDOUBLE sum = 0.0;
			for(int fatherLetter = 0; fatherLetter < alphabetSize; ++fatherLetter) 
			{
				sum += jointPost[mynode->id()][fatherLetter][letter];
			}
			_ancestralPost[mynode->id()][letter] = sum;
		}
	}
}
//calculates the joint probability P(N=x, father(N)=y|D)
//jointPost[nodeId][state at father][state here]
//In case of the root there is no meaning to Father=y
//-->put the posterior probability of x at the root in posteriorPerNodePer2States[root_id][0][x]
//also - return a filled suffStatGlobalHomPos used in later calculations
void chrNumberMng::computeJointPosterior(VVVdouble &jointPost, suffStatGlobalHomPos& sscUp)
{
	int numNodes = _tree.getNodesNum();
	int alphabetSize = _pSp->alphabetSize();
	jointPost.resize(numNodes);
	for (int n=0; n<jointPost.size();++n)
		resizeMatrix(jointPost[n], alphabetSize, alphabetSize);
	suffStatGlobalGamPos sscDown;
	sscUp.allocatePlace(numNodes, alphabetSize);
	computePijHom pi;
	pi.fillPij(_tree, *_pSp); 

	computeUpAlgChrNum comp_Up;
	computeDownAlg comp_Down;
	comp_Up.fillComputeUp(_tree, _sc, 0, pi, sscUp);
	comp_Down.fillComputeDownNonReversible(_tree, _sc, 0, pi, sscDown, sscUp);

	treeIterTopDownConst tIt(_tree);
	doubleRep ll = compLikelihood(pi, _tree, _sc, chrNumberOptions::_rootFreqType, _pSp);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		if (mynode->isRoot())
		{
			for(int letter=0; letter < alphabetSize; ++letter)
			{
				doubleRep prob = sscUp.get(mynode->id(), letter) * _pSp->freq(letter) / ll;
				jointPost[mynode->id()][0][letter] = convert(prob);
			}
			continue;
		}
		for (int sonState = 0; sonState<alphabetSize; ++sonState)
		{
			for (int fatherState = 0; fatherState<alphabetSize; ++fatherState)
			{
				doubleRep sum = 0.0; 
				for (int stateAtRoot = 0; stateAtRoot< alphabetSize; ++stateAtRoot)
				{
					sum +=	_pSp->freq(stateAtRoot)*
							sscDown.get(stateAtRoot, mynode->id(), fatherState)*
                            sscUp.get(mynode->id(),sonState)*
                            pi.getPij(mynode->id(), fatherState, sonState);
				}
				MDOUBLE p = convert(sum / ll);
				jointPost[mynode->id()][fatherState][sonState]= p;
			}
		}
	}
}



void chrNumberMng::printRunInfo()
{
	LOGnOUT(1, <<"chromEvol Version: 2.0. Last updated May 2020"<<endl;);
	LOGDO(5, Parameters::dump(myLog::LogFile())); 
	LOG(chrNumberOptions::_logValue, <<"\n ---------------------- THE PARAMETERS ----------------------------"<<endl;);
	if (chrNumberOptions::_mainType.size() > 0) 
		LOG(chrNumberOptions::_logValue, <<"main type: "<<chrNumberOptions::_mainType<<endl;);
	if (chrNumberOptions::_treeFile.size() > 0) 
		LOG(chrNumberOptions::_logValue, <<"tree file: "<<chrNumberOptions::_treeFile<<endl;);
	if (chrNumberOptions::_dataFile.size() > 0) 
		LOG(chrNumberOptions::_logValue, <<"data file: "<<chrNumberOptions::_dataFile<<endl;);
	if (chrNumberOptions::_outFile.size() > 0) 
		LOG(chrNumberOptions::_logValue, <<"output file: "<<chrNumberOptions::_outFile<<endl;);
	
	LOG(chrNumberOptions::_logValue, <<"model Type: "<<chrNumberOptions::getModelTypeStr(chrNumberOptions::_modelType)<<endl;);
	
	if 	(chrNumberOptions::_rootAt.size() > 0) 
		LOG(chrNumberOptions::_logValue, <<"root sequence name: "<<chrNumberOptions::_rootAt<<endl;);

	LOG(chrNumberOptions::_logValue, <<" max chromosome number allowed: "<< chrNumberOptions::_maxChrNum<<endl;);
	LOG(chrNumberOptions::_logValue, <<" _gainConstR: "<< chrNumberOptions::_gainConstR<<endl;);
	LOG(chrNumberOptions::_logValue, <<" _lossConstR: "<< chrNumberOptions::_lossConstR<<endl;);
	LOG(chrNumberOptions::_logValue, <<" _duplConstR: "<< chrNumberOptions::_duplConstR<<endl;);
	LOG(chrNumberOptions::_logValue, <<" _demiPloidyR: "<< chrNumberOptions::_demiPloidyR<<endl;);
	LOG(chrNumberOptions::_logValue, <<" _baseNumber: "<< chrNumberOptions::_baseNumber<<endl;);
	LOG(chrNumberOptions::_logValue, <<" _baseNumberR: "<< chrNumberOptions::_baseNumberR<<endl;);
	
	LOG(chrNumberOptions::_logValue, <<"\n -----------------------------------------------------------------"<<endl;);
}

void chrNumberMng::getStartingData()
{
	if (chrNumberOptions::_dataFile == "") {
		errorMsg::reportError("Please give a data file name in the parameters file");
	}
	 //the maximal chromosome permitted is maxInData + abs(chrNumberOptions::_maxChrNum)
	pair<int,int> min_max = chrCountFormat::getMinMaxCountInFile(chrNumberOptions::_dataFile);
	int maxCount = min_max.second;
	LOG(chrNumberOptions::_logValue, <<"max count = "<<maxCount<<" min count = "<<min_max.first<<endl);
	int maxChrAllowed = (chrNumberOptions::_maxChrNum <= 0)? abs(chrNumberOptions::_maxChrNum)+maxCount: abs(chrNumberOptions::_maxChrNum);
	if (maxChrAllowed < maxCount)
		errorMsg::reportError("the requested maximal chromosome number is lower than the maximal count observed in the input counts file");
	int minChrAllowed = (chrNumberOptions::_minChrNum <= 0)? min_max.first - abs(chrNumberOptions::_minChrNum): chrNumberOptions::_minChrNum;
	minChrAllowed = max(1,minChrAllowed);
	LOG(chrNumberOptions::_logValue, <<"max count allowed= "<<maxChrAllowed<<" min count allowed = "<<minChrAllowed<<endl);
	if (minChrAllowed > min_max.first)
		errorMsg::reportError("the requested minimal chromosome number is higher than the minimum count observed in the counts file");
	if (chrNumberOptions::_minBaseTransition  <= 0) {
		int minBase = min_max.first - abs(chrNumberOptions::_minBaseTransition);
		minBase = max(3,minBase);
		chrNumberOptions::_minBaseTransition = minBase;
		LOG(chrNumberOptions::_logValue, <<"min base transition = "<<minBase<<endl);
	}
	
	if (_pAlph)
		delete _pAlph;
	_pAlph = new finiteIntAlphabet(minChrAllowed,maxChrAllowed);
    ifstream inFile(chrNumberOptions::_dataFile.c_str());
	if (inFile.fail())
        errorMsg::reportError("cannot open file: " + chrNumberOptions::_dataFile);
	_sc = chrCountFormat::read(inFile, _pAlph);
	inFile.close();

	if (_pAlph->getCompositeIdsNum() > 0)
	{
		//the alphabet was changed
		for (sequenceContainer::taxaIterator it =_sc.taxaBegin(); it != _sc.taxaEnd(); ++it)
		{
			it->setAlphabet(_pAlph);
		}
	}
}


void chrNumberMng::getSimulatedTree(int iter)
{
	if (chrNumberOptions::_simulationsTreeDir.empty())
		return getTree();
	string treeFile = chrNumberOptions::_simulationsTreeDir + "/" + int2string(iter) + ".tree";
	chrNumberOptions::_treeFile = treeFile;
	getTree();
}

void chrNumberMng::getTree() {
	if (!(chrNumberOptions::_treeFile.empty()))
	{
		_tree = tree(chrNumberOptions::_treeFile);
        if (!_tree.withBranchLength()) 
            _tree.createFlatLengthMatrix(0.05);
	}
	else
	{
		VVdouble disTab;
		vector<string> vNames;
		jcDistance distanceMethod;
		giveDistanceTable(&distanceMethod, _sc, disTab, vNames);
		NJalg nj1;
        _tree= nj1.computeTree(disTab, vNames);
	}
	rootTree();
	MDOUBLE treeLength = getSumOfBranchLengths(_tree);
	LOG(3, <<"Original total tree length = "<<treeLength<<endl);
	//if (chrNumberOptions::_branchModelType == chrNumModel::SPECIATIONAL)
	//{//set all branches to equal length without changing the total tree length
	//	vector<tree::nodeP> leaves;
	//	_tree.getAllLeaves(leaves, _tree.getRoot());
	//	_tree.createFlatLengthMatrix(treeLength/leaves.size()); 
	//	treeLength = getSumOfBranchLengths(_tree);
	//}

	if (chrNumberOptions::_branchMul != 1.0) {
		if (chrNumberOptions::_branchMul == 999)
		{ //scale tree so that total tree length is equal to the number of different characters (a lower bound for the MP number of transitions+1)
			Vint type4Pos;
            getCharacterType4pos(_sc, type4Pos);
			type4Pos[0] += _pAlph->getCompositeIdsNum();
            int newTotalLength = type4Pos[0];
			MDOUBLE mul = newTotalLength / treeLength;
			_tree.multipleAllBranchesByFactor(mul);
			LOG(6, <<"rescaling tree by "<<mul<<" so that total tree length is "<<newTotalLength<<endl);
			chrNumberOptions::_branchMul = mul;
		}
		else
            _tree.multipleAllBranchesByFactor(chrNumberOptions::_branchMul);
		treeLength = getSumOfBranchLengths(_tree);
	}
	if (!DEQUAL(0.0, chrNumberOptions::_simulationsTreeLength, 0.0001))
	{
		MDOUBLE requestedLength = fabs(chrNumberOptions::_simulationsTreeLength);
		if (chrNumberOptions::_simulationsTreeLength < 0.0)
		{//_simulationsTreeLength indicates the average length from the root the the leaves
			vector<tree::nodeP> leaves;
			_tree.getAllLeaves(leaves, _tree.getRoot());
			MDOUBLE sum = 0.0;
			for (int i=0;i<leaves.size();++i) 
				sum += _tree.findLengthBetweenAnyTwoNodes(leaves[i], _tree.getRoot());
			MDOUBLE avg = sum / leaves.size();
			_tree.multipleAllBranchesByFactor(requestedLength / avg);
		}
		else // the finalTreeLength should be _simulationsTreeLength
            _tree.multipleAllBranchesByFactor(requestedLength / treeLength);
		treeLength = getSumOfBranchLengths(_tree);
	}
	if (chrNumberOptions::_branchModelType == chrNumModel::SPECIATIONAL)
	{//set all branches to equal length without changing the total tree length
		int edgesNum = _tree.getNodesNum() - 1;
		_tree.createFlatLengthMatrix(treeLength/edgesNum); 
		treeLength = getSumOfBranchLengths(_tree);
	}

	LOG(3, <<"total tree length = "<<treeLength<<endl);
	//print tree with internal nodes
	string outFileName = chrNumberOptions::_outDir + "//" + "allNodes.tree";
	//ofstream treeStream(outFileName.c_str());
	//printTreeStatesAsBPValues(treeStream, states, _tree.getRoot(), NULL, false, _pAlph);
	//treeStream.close();
	_tree.output(outFileName, tree::PHYLIP, true);
}

void chrNumberMng::initializeStatesVector(Vint& statesOut, const tree& inTree, const sequenceContainer& inSc)
{
	statesOut.resize(inTree.getNodesNum(),-2);
	bool bLeavesOnly = true;
	if (inTree.getNodesNum() == inSc.numberOfSeqs())
		bLeavesOnly = false;
	checkThatNamesInTreeAreSameAsNamesInSequenceContainer(inTree, inSc, bLeavesOnly);
	seqContainerTreeMap scTreeMap(inSc, inTree);	
	vector <tree::nodeP> leaves;
	inTree.getAllLeaves(leaves, inTree.getRoot());
	for (int i=0; i< leaves.size(); ++i){
		int myleafId = (leaves[i])->id();
		int mySeqId = scTreeMap.seqIdOfNodeI(myleafId);
		string name = (leaves[i])->name();
		statesOut[myleafId] = inSc[mySeqId][0];
	}
}

void chrNumberMng::changesNamesToLeaves(const tree& inTree, const sequenceContainer& inSc)
{
	seqContainerTreeMap scTreeMap(inSc, inTree);	
	vector <tree::nodeP> leaves;
	inTree.getAllLeaves(leaves, inTree.getRoot());
	for (int i=0; i< leaves.size(); ++i){
		int myleafId = (leaves[i])->id();
		int mySeqId = scTreeMap.seqIdOfNodeI(myleafId);
		int state = inSc[mySeqId][0];
		string newName = leaves[i]->name() + "-" + _pAlph->fromInt(state);
		leaves[i]->setName(newName);
	}
}

void chrNumberMng::printTreeStatesAsBPValues(ostream &out, Vint &states, const tree::nodeP &myNode, 
							   VVVdouble *probs,bool printGains, const alphabet* pAlph)  {
	if (myNode->isLeaf()) {
		//out << myNode->name()<< ":"<<myNode->dis2father();
		out << myNode->name()<<"-"<<pAlph->fromInt(states[myNode->id()])<< ":"<<myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (int i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			printTreeStatesAsBPValues(out,states,myNode->getSon(i),probs, printGains, pAlph);
		}
		out <<")";
		//if (myNode->isRoot()==false) {
			//out<<states[myNode->id()]<<"--";
			//out<<myNode->name();
			out.precision(3);
			if (probs){
				if (printGains)
					out<<(*probs)[myNode->id()][0][2]<<"//"<<(*probs)[myNode->id()][1][2]; 
				else //print losses
					out<<(*probs)[myNode->id()][2][0]<<"//"<<(*probs)[myNode->id()][2][1]; 
			}
			//out << "["<<myNode->name()<<"]"; 
			//out<<":"<<myNode->dis2father();
			out << "["<<myNode->name()<<"-"<<pAlph->fromInt(states[myNode->id()])<<"]"; 
			out<<":"<<myNode->dis2father();

		//}
	}
}

void chrNumberMng::rootTree() {

	tree::nodeP myroot;
	if (!(chrNumberOptions::_rootAt =="")){
		myroot = _tree.findNodeByName(chrNumberOptions::_rootAt); //returns NULL if not found
		if (myroot){
			_tree.rootAt(myroot);
		}
	}
	myroot = _tree.getRoot();
	LOG(chrNumberOptions::_logValue, <<"tree rooted at "<<myroot->name()<<" id, "<<myroot->id()<<endl);
	LOG(chrNumberOptions::_logValue, <<"sons of root are: "<<endl);
	for (int i=0;i<myroot->getNumberOfSons();++i) 
		LOG(chrNumberOptions::_logValue, <<_tree.getRoot()->getSon(i)->name()<<endl);
}

void chrNumberMng::initStartingModel(const finiteIntAlphabet* pAlph, int allowedRange){
	Vdouble freqs(pAlph->size(), 1.0 / pAlph->size());
	if (!(chrNumberOptions::_freqFile.empty()))
	{
		readFreqsFromFile(chrNumberOptions::_freqFile, freqs, pAlph);
	}
	
	chrNumModel* pModel = NULL;
	switch (chrNumberOptions::_modelType)
	{
	case chrNumModel::GENERAL_CHR_MODEL:
		{
		vector<paramObj> initParams;
		chrNumberOptions::getInitParams(initParams);
		if (chrNumberOptions::_baseTransitionProbs != "") {
			vector<pair<int, MDOUBLE> > compProbs;
			getBaseTransitions(chrNumberOptions::_baseTransitionProbs, compProbs);
			pModel = new generalChrNumModel(pAlph, freqs, chrNumberOptions::_rootFreqType, initParams, compProbs);
		}
		else {
			pModel = new generalChrNumModel(pAlph, freqs, chrNumberOptions::_rootFreqType, initParams, allowedRange);
		}
		break;
		}
	default:
		errorMsg::reportError("unknown model type in chrNumberMng::initStartingModel");
	}

	trivialAccelerator pijAcc(pModel);
    uniDistribution uniDistr;
	if (_pSp)
	{
		delete _pSp;
		_pSp = NULL;
	}
    _pSp = new stochasticProcess(&uniDistr, &pijAcc, false);
	delete pModel;
}


class valCmp {
public:
	bool operator()(const pair<int, MDOUBLE> & elem1, const pair<int, MDOUBLE> & elem2) {
		return (elem1.first < elem2.first);
	}
};



void chrNumberMng::getBaseTransitions(const string& baseTransitionProbs, vector<pair<int, MDOUBLE> > & compProbs) {
	if ((chrNumberOptions::_baseNumber == IGNORE_PARAM) || (chrNumberOptions::_baseNumberR == IGNORE_PARAM)) 
		errorMsg::reportError("initialization of possible base transitions is only possible when the baseNumber parameter is included in the model");
	vector<string> strVec;
	MDOUBLE sum = 0.0;
    splitString(baseTransitionProbs, strVec, "_"); 
	for (int i = 0; i < strVec.size(); ++i)
	{
		string transitionStr, probStr;
		splitString2(strVec[i], "=", transitionStr, probStr);
		if (probStr == "NULL")
			probStr = "1.0";
		MDOUBLE prob = atof(probStr.c_str());
		int transition = atoi(transitionStr.c_str());
		if (transition % chrNumberOptions::_baseNumber != 0)
			errorMsg::reportError("possible transitions must be multiplications of _baseNumber");
		compProbs.push_back(pair<int, MDOUBLE>(transition , prob));
		sum += prob;
	}
	if (!DEQUAL(sum, 1.0))
		errorMsg::reportError("The possible base transition probabilities should sum to 1.0");
	sort(compProbs.begin(), compProbs.end(), valCmp());
}

pair<int, int> chrNumberMng::getChrRangeInData(const sequenceContainer& sc)
{
	int maxInData = 0;
	int minInData = MAX_CHR_POSSIBLE;
	sequenceContainer::constTaxaIterator seqItr = sc.constTaxaBegin();
    for (; seqItr != sc.constTaxaEnd(); ++seqItr)
	{
		for (int pos = 0; pos < sc.seqLen(); ++pos)
		{
			int id = (*seqItr)[pos];
			//should check if this is the corrected ID
			if (id == _pAlph->unknown()) 
				continue;
			if (_pAlph->isSpecific(id))
			{
				if (maxInData < id)
					maxInData = id;
				if (minInData > id)
					minInData = id;
			}
			else
			{//composite
				Vint compIds;
				Vdouble compProbs;
				_pAlph->getCompositeParts(id, compIds, compProbs);
				for (int c = 0; c < compIds.size(); ++c)
				{
					if (compIds[c] > maxInData)
						maxInData = compIds[c];
					if (compIds[c] < minInData)
						minInData = compIds[c];
				}
			}
		}
	}
	pair<int, int> res(_pAlph->id2Count(minInData), _pAlph->id2Count(maxInData));
	return res;
}


/*
Note that not all counts must be in the input file. Missing (but allowed) counts will be set to zero.
input file looks like this:
F[1]=0.1
F[2]=0.15
...
F[MAX]=0.002
*/
void chrNumberMng::readFreqsFromFile(const string& inFileName, Vdouble& freqs, const finiteIntAlphabet* pAlph)
{
	int countNum = pAlph->size();
	freqs.clear();
	freqs.resize(countNum, 0.0);
	ifstream inF(inFileName.c_str());
	vector<string> lines;
	putFileIntoVectorStringArray(inF, lines);
//	if (lines.size() != maxChrNum)
//		errorMsg::reportError("error in chrNumberMng::readFreqsFromFile the input frequency vector is not equal to the maximal chromosome number allowed!");
	MDOUBLE sumF = 0.0;
	int n = 0;
	for (int l = 0; l < lines.size(); ++l)
	{
		if (lines[l].empty())
			continue;
		if ((lines[l][0] != 'F') || (lines[l].find("=")  == -1))
			continue;
		string number, rest;
		splitString2(lines[l], "=", rest, number);
		number = trim(number);
		int startNum = 1+rest.find_first_of("[");
		int endNum = rest.find_first_of("]");
		string numStr = rest.substr(startNum,endNum - startNum);
		int count = atoi(numStr.c_str());
		if ((count > pAlph->max()) || (count < pAlph->min()))
			errorMsg::reportError("cannot set root frequency to count = " + int2string(count));
		int placeInAlph = count - pAlph->min();
		freqs[placeInAlph] = string2double(number.c_str());
		sumF += freqs[placeInAlph];
	}

	if (!DEQUAL(1.0, sumF, 0.0001))
		errorMsg::reportError("error in chrNumberMng::readFreqsFromFile. sum of input frequencies should sum to 1.0");
	if (!DEQUAL(1.0, sumF))
		scaleVec(freqs, 1.0/freqs.size());
	makeSureNoZeroFreqs(freqs);
	inF.close();
}


MDOUBLE chrNumberMng::optimizeParams(bool bScaleTree)
{
	MDOUBLE res = 0.0;
	//chrNumModel* pModel = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel());
	if (chrNumberOptions::_mainType == "Run_Fix_Param")
	{
		if (bScaleTree == true) {
			if ((chrNumberOptions::_branchModelType == chrNumModel::COMBINED) || (chrNumberOptions::_branchModelType == chrNumModel::COMBINED_INTERNALS))
			{
				MDOUBLE scaleF = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel())->getScaleFactor();
				scaleTree(_tree, scaleF);
			}
		}

		res = compLogLikelihood(_pSp, _tree, _sc, chrNumberOptions::_rootFreqType);
		LOG(3, <<"LL = "<<res<<endl;);
		return res;
	}
	else if (chrNumberOptions::_mainType == "Optimize_Model")
	{
		LOGnOUT(5,<<"Optimizing parameters"<<endl);
		chrNumOptimizer opt(_tree, _sc);
		Vdouble tols(chrNumberOptions::_optimizeIterNum.size(), chrNumberOptions::_epsilonLLimprovement);
		if (tols.size() > 2)
            tols[1] = chrNumberOptions::_epsilonLLimprovement * 2;
		res = opt.optimizeModelManyStart(_pSp, chrNumberOptions::_optimizePointsNum, chrNumberOptions::_optimizeIterNum, tols);
		LOG(3, <<"after optmizations"<<endl;);
		if (bScaleTree == true) {
			if ((chrNumberOptions::_branchModelType == chrNumModel::COMBINED) || (chrNumberOptions::_branchModelType == chrNumModel::COMBINED_INTERNALS))
			{
				MDOUBLE scaleF = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel())->getScaleFactor();
				scaleTree(_tree, scaleF);
			}
		}
		//LOGDO(3, pModel->printModelParams(myLog::LogFile())); 
	}
	return res;
}

void chrNumberMng::printResults(const VVVdouble& jointPost)
{
	chrNumModel* pModel = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel());
    string fileName = chrNumberOptions::_outDir + "//" + chrNumberOptions::_outFile;
	ofstream outF(fileName.c_str());
	pair<int,int> chrRange = getChrRangeInData(_sc);
	outF<<"#Model: "<<chrNumberOptions::getModelTypeStr(chrNumberOptions::_modelType)<<endl;
	outF<<"#Branching Model: "<<chrNumberOptions::getBranchModelTypeStr(chrNumberOptions::_branchModelType)<<endl;
	outF<<"#Input data file: "<<chrNumberOptions::_dataFile<<endl;
	outF<<"#Input tree file: "<<chrNumberOptions::_treeFile<<endl;
	tree::nodeP myroot = _tree.getRoot();
	outF<<"#Tree rooted at "<<myroot->name()<<"\t";
	outF<<"sons of root are: ";
	for (int i=0;i<myroot->getNumberOfSons();++i) 
		outF<<_tree.getRoot()->getSon(i)->name()<<", ";
	outF<<endl;
	outF<<"#total tree length = "<<getSumOfBranchLengths(_tree)<<endl;
	outF<<"#min chromosome in data: "<<chrRange.first<<endl;
	outF<<"#min chromosome allowed: "<<_pAlph->min()<<endl;
	outF<<"#max chromosome in data: "<<chrRange.second<<endl;
	outF<<"#max chromosome allowed: "<<_pAlph->max()<<endl;
	outF<<"#Half_duplication rate is ";
	if (chrNumberOptions::_demiPloidyR > 0)
		outF<<"estimated"<<endl;
	else if (chrNumberOptions::_demiPloidyR == DEMI_EQUAL_DUPL)
		outF<<"same as duplication rate"<<endl;
	else
		outF<<"zero"<<endl;
	outF<<"#Root frequencies: "<<chrNumberOptions::getRootFreqTypeStr(chrNumberOptions::_rootFreqType)<<endl;
	MDOUBLE branchMul = chrNumberOptions::_branchMul;
	if (chrNumberOptions::_branchMul == 999)
	{
		Vint type4Pos;
        getCharacterType4pos(_sc, type4Pos);
        int newTotalLength = type4Pos[0];
		outF<<"#Tree was rescaled so that the total tree length is "<< newTotalLength<<endl;
		outF<<"#The model parameters were computed according to the scaled tree!"<<endl;
		branchMul  = 1.0;
	}
	else if (chrNumberOptions::_branchMul != 1.0) {
		outF<<"#For optimization issues the tree branches were multiplied by "<<chrNumberOptions::_branchMul<<endl;
		outF<<"#To preserve the original time unit the model parameters should be multiplied by the same factor !!!"<<endl;
	}
	//print params
	outF<<"#Final Model parameters"<<endl;
	vector<paramObj> paramVec;
	pModel->getParams(paramVec);
	for(int p = 0; p < paramVec.size(); ++p) {
		if (paramVec[p]._type == SCALE_BRANCH)
            outF<<pModel->getParamName(paramVec[p]._type)<<"\t"<<paramVec[p]._value<<endl;
		else
			outF<<pModel->getParamName(paramVec[p]._type)<<"\t"<<paramVec[p]._value<<endl;
			//outF<<pModel->getParamName(paramVec[p]._type)<<"\t"<<paramVec[p]._value * branchMul<<endl;
	}
	
	//print frequencies:
	//outF.precision(3);
	MDOUBLE aic = (-2 * _bestLL) + 2* paramVec.size();
	for (int c = 0; c < pModel->alphabetSize(); ++c)
		outF<<"F["<<_pAlph->id2Count(c)<<"]="<<pModel->freq(c)<<"   ";
	outF<<endl;
	outF << "LogLikelihood = "<<_bestLL<<endl;
	outF << "Likelihood = "<<exp(_bestLL)<<endl;
	outF << "AIC (Akaike information criterion) = "<<aic<<endl;
	outF.close();

	string expFileName = chrNumberOptions::_outDir + "//expectations.txt";
	ofstream expFile(expFileName.c_str());
	printExpectations(jointPost, expFile);
	expFile.close();
	string mlTreeStr = chrNumberOptions::_outDir + "//" + chrNumberOptions::_inferTreeFile;
	string postTreeStr = chrNumberOptions::_outDir + "//" + "posteriorAncestors.tree";
	printReconstructedTrees(mlTreeStr, postTreeStr);
}


//this is similar to likelihoodComputation::getLofPos(int pos,tree& et,sequenceContainer& sc, const computePijHom& pi,const stochasticProcess& sp)
//but here the frequencies are allowed to change after the computation is done
//so that freq[i] = Prob(data| i at root) / Prob(data)
MDOUBLE chrNumberMng::compLogLikelihood(const stochasticProcess* pSp, const tree &tr, const sequenceContainer &sc, chrNumModel::rootFreqType freqType)
{
	if (sc.seqLen() != 1)
		errorMsg::reportError("are you sure seqLen != 1. have to change LL computation");
	computePijHom pij;
	pij.fillPij(tr, *pSp, 0, pSp->isReversible()); 
	MDOUBLE res = compLogLikelihood(pij, tr, sc, freqType, pSp);
	return res;
}

doubleRep chrNumberMng::compLikelihood(const computePijHom& pij, const tree &tr, const sequenceContainer &sc, chrNumModel::rootFreqType freqType, const stochasticProcess* pSp)
{
	computeUpAlgChrNum cup;
	suffStatGlobalHomPos ssc;
	cup.fillComputeUp(tr, sc, 0, pij, ssc);
	if (freqType == chrNumModel::ROOT_LL)
	{
		doubleRep sum = 0.0;
		int let = 0;
		for (; let < pSp->alphabetSize(); ++let) 
		{
			sum += ssc.get(tr.getRoot()->id(),let);
		}
		if (!(sum>0.0))
		{
			LOGnOUT(3, <<"sum = "<<sum<<endl);
			LOGDO(3, static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel())->printModelParams(myLog::LogFile())); 
			errorMsg::reportError("in chrNumberMng::compLogLikelihood: sum of likelihood in root over all letters is ZERO");
		}
		chrNumModel* pModel = static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel());
		Vdouble freqs(pSp->alphabetSize());

		MDOUBLE sumF = 0.0;
		for (let = 0; let < pSp->alphabetSize(); ++let)
		{
			freqs[let] = convert(ssc.get(tr.getRoot()->id(),let) / sum);
			sumF += freqs[let];
			//MDOUBLE val = convert(ssc.get(tr.getRoot()->id(),let) / sum);
			//pModel->setFreq(let, val);
		}
		if (!DEQUAL(1.0, sumF))
		{
            cerr<<"SUMFFFF != 1.0!!! sumF = "<<sumF<<endl;
			scaleVec(freqs, 1.0/freqs.size());
		}
        makeSureNoZeroFreqs(freqs);
		for (let = 0; let < pSp->alphabetSize(); ++let)
			pModel->setFreq(let, freqs[let]);
	}

	doubleRep res = 0.0;

	for (int let = 0; let < pSp->alphabetSize(); ++let) {
        doubleRep tmpLcat = ssc.get(tr.getRoot()->id(),let)* pSp->freq(let);
		if (!DBIG_EQUAL(convert(tmpLcat), 0.0))
		{
			cerr<<"tmpLcat = "<<convert(tmpLcat)<<" letter = "<< let<<" freq(let) = "<<pSp->freq(let)<<endl;
			errorMsg::reportError("error in likelihoodComputation::getLofPos. likelihood is smaller than zero");
		}
		res +=tmpLcat;
	}

	if (!(res>0.0)){
		LOG(5,<<"likelihood of pos was zero!!!"<<endl;);
		LOG(5,<<"likelihoodComputation::getLofPos: "<< res<<endl;);
	}
	return res;
}

MDOUBLE chrNumberMng::compLogLikelihood(const computePijHom& pij, const tree &tr, const sequenceContainer &sc, chrNumModel::rootFreqType freqType, const stochasticProcess* pSp)
{
	doubleRep l = compLikelihood(pij, tr, sc, freqType, pSp);
	return convert(log(l));
}

void chrNumberMng::ancestralReconstructML(Vint& statesOut)
{
	initializeStatesVector(statesOut, _tree, _sc);
	VVdouble upL;
	VVint backtrack;
	traverseUpML(upL, backtrack, statesOut);
	VVint transitionTypeCount;
	MDOUBLE LofJoint = traverseDownML(upL, backtrack, statesOut, transitionTypeCount);
}


// upL[node][father_letter] = best reconstruction of subtree with node as the root given that the state at the node's father is father_letter
// -->upL[node][father_letter]=max(letter_here){P(father_letter->letter_here)*upL[son1][letter_here]*upL[son2][letter_here]} 
// backtrack[node][father_letter] = argmax of above 
void chrNumberMng::traverseUpML(VVdouble &upL, VVint &backtrack ,const Vint& states)
{
	resizeMatrix<MDOUBLE>(upL, _tree.getNodesNum(), _pSp->alphabetSize());
	resizeMatrix<int>(backtrack, _tree.getNodesNum(), _pSp->alphabetSize());
	
	computePijHom pij;
	pij.fillPij(_tree, *_pSp, 0, _pSp->isReversible()); 
	
	treeIterDownTopConst tIt(_tree);
	int father_state;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int myState = states[mynode->id()];
		if (mynode->isLeaf()) {
			for (father_state=0; father_state< _pSp->alphabetSize(); ++father_state)
			{ // looping over states at father
				MDOUBLE totalVal = 0.0;
				for (int let=0; let< _pSp->alphabetSize();++let) 
				{
					MDOUBLE val = static_cast<const finiteIntAlphabet*>(_pAlph)->relationsProbs(myState, let);
					//MDOUBLE val = _pAlph->relations(myState, let);
					if (val>0) 
					{
						if (val < 0.99)
							MDOUBLE x = 1;
						val *= pij.getPij(mynode->id(), father_state, let);
						totalVal +=val;
					}
				}
				upL[mynode->id()][father_state]= totalVal;
				backtrack[mynode->id()][father_state] = myState;
			}
		}
		else if (!(mynode->isRoot())) { //internal node
			for (father_state=0; father_state<_pSp->alphabetSize(); ++father_state)
			{ // looping over states at father. For each has to compute Up[myNode][fatherState]
				MDOUBLE myMax = -1;
				int myArgMax=-1;
				for (int my_state=0; my_state<_pSp->alphabetSize(); ++my_state)
				{ // loop to find max over current node
					MDOUBLE val=_pSp->Pij_t(father_state, my_state, mynode->dis2father());
					for (int son=0;son<mynode->getNumberOfSons();son++)
						val *= upL[mynode->getSon(son)->id()][my_state];
					if (val>myMax){
						myMax = val;
						myArgMax = my_state;
					}
				}
				if ((myMax<0) || (myArgMax<0))
					errorMsg::reportError("Error in traverseUpML: cannot find maximum");
				upL[mynode->id()][father_state] = myMax;
				backtrack[mynode->id()][father_state] = myArgMax;
			}
		}
		else {// root
			for (int root_state=0; root_state<_pSp->alphabetSize(); ++root_state){ 
				MDOUBLE val = _pSp->freq(root_state);
				for (int son=0; son<mynode->getNumberOfSons(); ++son)
	                val *= upL[mynode->getSon(son)->id()][root_state];
				upL[mynode->id()][root_state]=val;
			}
		}
	}
}


//fill the vectors reconstructStates (the states at the leaves are already filled) and transitionTypeCount
//returns likelihood of max joint reconstruction. 
//Assumes that upL and backtrack are already filled.
MDOUBLE chrNumberMng::traverseDownML(const VVdouble &upL, const VVint &backtrack, Vint& reconstructStates, VVint& transitionTypeCount) { 
	
	if (backtrack.size() != _tree.getNodesNum()) 
		errorMsg::reportError("error in chrNumberMng::traverseDownML, input vector backtrack must be filled (call traverseUpML() first)");
	transitionTypeCount.resize(_pSp->alphabetSize());
	for (int i = 0; i < transitionTypeCount.size(); i++) 
		transitionTypeCount[i].resize(_pSp->alphabetSize(),0);

	MDOUBLE LofJoint;
	int stateOfRoot;
	findMaxInVector(upL[(_tree.getRoot())->id()], LofJoint, stateOfRoot);
	reconstructStates[(_tree.getRoot())->id()] = stateOfRoot;
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isRoot()) 
			continue;
		int myId = mynode->id();
		int stateAtFather = reconstructStates[mynode->father()->id()];
		if (mynode->isLeaf()) {
			if (!_pAlph->isComposite(reconstructStates[mynode->id()])&& (reconstructStates[mynode->id()] != _pAlph->unknown()))
                ++transitionTypeCount[stateAtFather][reconstructStates[mynode->id()]];
			if ((reconstructStates[mynode->id()] != stateAtFather))
				LOG(10, <<"switch from "<<mynode->father()->name()<<"("<<stateAtFather<<") to "<<mynode->name()<<"("<<reconstructStates[mynode->id()]<<")"<<endl;);
			continue;
		}
		reconstructStates[myId] = backtrack[myId][stateAtFather];
		++transitionTypeCount[stateAtFather][reconstructStates[myId]];
	}
	return log(LofJoint);
}

void chrNumberMng::getTransitionVec(const tree& tr, const Vint& treeStates, Vint& upChangeVec, Vint& downChangeVec)
{
	upChangeVec.resize(_pSp->alphabetSize(), 0); //upChangeVec[i] = number of branches in which i upward changes occured (i=0--> no change)
	downChangeVec.resize(_pSp->alphabetSize(), 0); //downChangeVec[i] = number of branches in which i downward changes occured (i=0--> no change)
	treeIterTopDownConst tIt(tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
        if (mynode->isRoot()) 
			continue;
		int stateAtFather = treeStates[mynode->father()->id()];
		int myState = treeStates[mynode->id()];
		int diff = myState - stateAtFather;
		if (diff > 0)
			upChangeVec[diff]++;
		else if (diff < 0)
			downChangeVec[-diff]++;
		else
		{ //no change
			upChangeVec[diff]++;
			downChangeVec[diff]++;
		}
	}
	LOGDO(3, printVec(upChangeVec, myLog::LogFile())); 
}


//infer the expected number of changes between two states along the whole tree 
//AND the probability for each type of change for each branch. 
//receives the joint posterior p(N=x, Father(N)=y|D)
//results are stored in members _expChanges[nodeId][0-4] 
//0 = gain along the branch from nodeId to its father, 1 = lost, 2 = dupl, 3 = demi-dupl 4=base_number
//NOTE: THE COMPUTATION OF THE EVENTS PROBABILITIES IS NOT CORRECT AS IT SUMS OVER NON-MUTUALLY-EXCLUSIVE EVENTS:
//-->PROB[GAIN] = SUM{PROB X->X+1} BUT IT IS POSSIBLE THAT OVER THE SAME BRANCH IN THE SIMULATIONS WE WILL GET A SERIES OF EVENTS X->X+1->X+2
//WHEN WE SUM THE PROBS, THE PROB CAN BE > 1 AS BOTH X->X+1 AND X+1->X+2 WILL BE COUNTED.
void chrNumberMng::computeChangeProbAndExp(const VVVdouble& jointPost, const suffStatGlobalHomPos& sscUp)
{
	if (chrNumberOptions::_smIter > 0)
		return computeChangeProbAndExp_SM(jointPost, sscUp);
	int maxCount = _pAlph->max();
	int minCount = _pAlph->min();
	Vstring expPrintData(_tree.getNodesNum());
	MDOUBLE gainExp = 0.0, lossExp = 0.0, duplExp = 0.0, halfDuplExp = 0.0, maxChrExp = 0.0, baseExp=0.0;
	_expChanges.clear();
	//_probChanges.clear();
	_expChanges.resize(_tree.getNodesNum()+1); //last place is total across tree
	//_probChanges.resize(_tree.getNodesNum());
	int n;
	for (n = 0; n < _expChanges.size(); ++n)
	{
		_expChanges[n].resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
		//_probChanges[n].resize(chrNumModel::JUMP_TYPE_MAX, 0.0);
	}
	
	//this is done after the resize since we need this matrix later on even of simulations weren't performed
	if (!performSimulationsForExpectation())
	{
		if (chrNumberOptions::_simulationsNum < 1) {
			LOGnOUT(1, <<"Cannot perform simulations. The input number of simulations is lower than 1"<<endl);
		}
		else 
			LOGnOUT(1, <<"Cannot perform simulations. The alpahabet size of "<<_pAlph->size()<<" is too high"<<endl);
		return;
	}

	LOGnOUT(5,<<endl<<"running "<<chrNumberOptions::_simulationsNum<<" simulations"<<endl);
	simulateJumpsLargeAlphabet sim(_tree, *_pSp, _pAlph);
	sim.runSimulation(chrNumberOptions::_simulationsNum);
	LOGnOUT(5, <<"finished simulations"<<endl);
	if ((chrNumberOptions::_baseNumber != IGNORE_PARAM) && (chrNumberOptions::_baseNumberR != IGNORE_PARAM)) {
		computeChangeProbAndExpBaseNumber(sim, jointPost);
	}

	const chrNumModel* pModel = getModel();
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		int chr;
		for (chr = minCount; chr < maxCount-1; ++chr)
		{
			int state = _pAlph->count2Id(chr);
			if (pModel->getGainR(state) > THRESHOLD_R)
				_expChanges[nodeId][chrNumModel::GAIN_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, state+1);
			//_probChanges[nodeId][chrNumModel::GAIN_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, chr);
			if ((chr != minCount) && (pModel->getLossR(state) > THRESHOLD_R))
                _expChanges[nodeId][chrNumModel::LOSS_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, state-1);
			//_probChanges[nodeId][chrNumModel::LOSS_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, chr-2);
			_expChanges[nodeId][chrNumModel::MAX_CHR_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(maxCount));
			//_probChanges[nodeId][chrNumModel::MAX_CHR_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, maxChr-1);
			//compute exp of a duplication 
			if (chr <= maxCount/2) 
			{
				_expChanges[nodeId][chrNumModel::DUPL_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(2*chr));
				//_probChanges[nodeId][chrNumModel::DUPL_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, 2*chr-1);
			}
			//compute exp of a demi-duplication 
			if ((chr <= maxCount * 0.66) && (chr > 2) && (pModel->getDemiR(state) > THRESHOLD_R))
			{
				if (chr == 3)
				{
					_expChanges[nodeId][chrNumModel::DEMI_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(5)); //count only 3->5
					//_probChanges[nodeId][chrNumModel::DEMI_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, 2, 4);
				}
				else if ((chr % 2) == 0)
				{
					_expChanges[nodeId][chrNumModel::DEMI_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(static_cast<int>(1.5*chr)));
					//_probChanges[nodeId][chrNumModel::DEMI_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, static_cast<int>(1.5*chr-1));
				}
				else
				{
					_expChanges[nodeId][chrNumModel::DEMI_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(static_cast<int>(1.5*chr)));
					_expChanges[nodeId][chrNumModel::DEMI_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(static_cast<int>(1.5*chr) + 1));
					//_probChanges[nodeId][chrNumModel::DEMI_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, static_cast<int>(1.5*chr-1));
					//_probChanges[nodeId][chrNumModel::DEMI_J] += computeProbOfChangePerBranch(sim, jointPost, mynode, chr-1, static_cast<int>(1.5*chr));
				}
			}
		}
		//for (maxCount -1) there are only two possibilities : loss and jumps to maxChr
		int state = _pAlph->count2Id(maxCount-1);
		if (pModel->getLossR(state) > THRESHOLD_R)
			_expChanges[nodeId][chrNumModel::LOSS_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(maxCount-2));
		_expChanges[nodeId][chrNumModel::MAX_CHR_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(maxCount));
		//for maxChr there is only one possibility loss 
		state = _pAlph->count2Id(maxCount);
		if (pModel->getLossR(state) > THRESHOLD_R)
			_expChanges[nodeId][chrNumModel::LOSS_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(maxCount-1));
		//add all results
		gainExp += _expChanges[nodeId][chrNumModel::GAIN_J];
		lossExp += _expChanges[nodeId][chrNumModel::LOSS_J];
		duplExp += _expChanges[nodeId][chrNumModel::DUPL_J];
		halfDuplExp += _expChanges[nodeId][chrNumModel::DEMI_J];
		maxChrExp += _expChanges[nodeId][chrNumModel::MAX_CHR_J];
		if (mynode->isLeaf())
			expPrintData[nodeId] = "[" + double2string(_expChanges[mynode->id()][chrNumModel::GAIN_J]) + "//" + double2string(_expChanges[mynode->id()][chrNumModel::LOSS_J]) + "//" +  double2string(_expChanges[mynode->id()][chrNumModel::DUPL_J]) + "//" +  double2string(_expChanges[mynode->id()][chrNumModel::DEMI_J]) + "//"+ "]";
		else
			expPrintData[nodeId] = "[" + mynode->name() + "-" + double2string(_expChanges[mynode->id()][chrNumModel::GAIN_J]) + "//" + double2string(_expChanges[mynode->id()][chrNumModel::LOSS_J]) + "//" +  double2string(_expChanges[mynode->id()][chrNumModel::DUPL_J]) + "//" +  double2string(_expChanges[mynode->id()][chrNumModel::DEMI_J]) + "//"+ "]";
	}	
	LOGnOUT(3,<<endl<<"total expectations"<<endl);
	LOGnOUT(3,<<"dupl="<<duplExp<<endl);
	LOGnOUT(3,<<"gain="<<gainExp<<endl);
	LOGnOUT(3,<<"loss="<<lossExp<<endl);
	LOGnOUT(3,<<"halFDupl="<<halfDuplExp<<endl);
	LOGnOUT(3,<<"baseNumber="<<_expChanges[_expChanges.size()-1][chrNumModel::BASE_J]<<endl);
	LOGnOUT(3,<<"toMaxChr="<<maxChrExp<<endl);
	_expChanges[_expChanges.size()-1][chrNumModel::GAIN_J] = gainExp;
	_expChanges[_expChanges.size()-1][chrNumModel::LOSS_J] = lossExp;
	_expChanges[_expChanges.size()-1][chrNumModel::DUPL_J] = duplExp;
	_expChanges[_expChanges.size()-1][chrNumModel::DEMI_J] = halfDuplExp;
	
	//print tree with expectations as BP values
	//have to change the name of the leaves because can't print the expectations as bootstrap
	string expTeeStr = chrNumberOptions::_outDir + "//" + "exp.tree";
	ofstream expTreeStream(expTeeStr .c_str());
	tree printTree(_tree);
	seqContainerTreeMap scTreeMap(_sc, printTree);	
	vector <tree::nodeP> leaves;
	printTree.getAllLeaves(leaves, printTree.getRoot());
	for (int i=0; i< leaves.size(); ++i){
		int myleafId = (leaves[i])->id();
		int mySeqId = scTreeMap.seqIdOfNodeI(myleafId);
		string newName = leaves[i]->name() + "-" + expPrintData[myleafId];
		leaves[i]->setName(newName);
	}
	printDataOnTreeAsBPValues(expTreeStream, expPrintData, printTree);
	expTreeStream.close();
}

void chrNumberMng::computeChangeProbAndExp_SM(const VVVdouble& jointPost, const suffStatGlobalHomPos& sscUp)
{
	simulateChangesAlongTree simJumps(_tree, *_pSp, _pAlph);
	_expChanges = simJumps.generateStochasticMapping(_sc, chrNumberOptions::_smIter, jointPost, sscUp);

	LOGnOUT(3,<<endl<<"total expectations"<<endl);
	LOGnOUT(3,<<"dupl="<<_expChanges[_expChanges.size()-1][chrNumModel::DUPL_J]<<endl);
	LOGnOUT(3,<<"gain="<<_expChanges[_expChanges.size()-1][chrNumModel::GAIN_J]<<endl);
	LOGnOUT(3,<<"loss="<<_expChanges[_expChanges.size()-1][chrNumModel::LOSS_J]<<endl);
	LOGnOUT(3,<<"halFDupl="<<_expChanges[_expChanges.size()-1][chrNumModel::DEMI_J]<<endl);
	LOGnOUT(3,<<"baseNumber="<<_expChanges[_expChanges.size()-1][chrNumModel::BASE_J]<<endl);
	LOGnOUT(3,<<"toMaxChr="<<_expChanges[_expChanges.size()-1][chrNumModel::MAX_CHR_J]<<endl);
}

void chrNumberMng::computeChangeProbAndExpBaseNumber(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost)
{
	int maxCount = _pAlph->max();
	int minCount = _pAlph->min();
	int baseNumber = static_cast<chrNumModel*>(_pSp->getPijAccelerator()->getReplacementModel())->getBaseNumber();

	MDOUBLE baseExp = 0.0;
	treeIterTopDownConst tIt(_tree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
	{
		if (mynode->isRoot())
			continue;
		int nodeId = mynode->id();
		int chr;
		for (chr = minCount; chr < maxCount-1; ++chr)
		{
			//for each chr go over all possible baseNumber transitions
			int state = _pAlph->count2Id(chr);
			int toChrom = chr + baseNumber;
			for (; toChrom < maxCount; toChrom += baseNumber)
			{
				if (chr*2 == toChrom )//this will go into the duplication count. Note that a similar modification should be made for DEMI, but at present we do not allow for both demi-ploidy and base-number transitions 
					continue;
				_expChanges[nodeId][chrNumModel::BASE_J] += computeExpectationOfChangePerBranch(sim, jointPost, mynode, state, _pAlph->count2Id(toChrom));
			}
		}

		baseExp += _expChanges[nodeId][chrNumModel::BASE_J];
	}	
	_expChanges[_expChanges.size()-1][chrNumModel::BASE_J] = baseExp;
}


bool chrNumberMng::performSimulationsForExpectation()
{
	if ((chrNumberOptions::_simulationsNum < 1) && (chrNumberOptions::_smIter < 1))
		return false;
	int maxChr = _pSp->alphabetSize();
    if (maxChr > MAX_CHR_FOR_SIMULATIONS)
	{
		return false;
	}
	return true;

}

//Expectation of number of changes from fromState to toState along the whole tree 
//NF=sum over all changes x,y:
//Posterior(Node=x,Father=y|D)*Exp(changes u to v|Node=x,Father=y)
//The second term is given to the function as input
//MDOUBLE chrNumberMng::computeExpectationOfChange(simulateJumpsLargeAlphabet &sim,const	VVVdouble &jointPost, int fromState, int toState)
//{
//	MDOUBLE res = 0;
//	treeIterTopDownConst tIt(_tree);
//	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
//		MDOUBLE exp = computeExpectationOfChangePerBranch(sim, jointPost, mynode, fromState, toState);
//		res += exp;
//	}
//	return res;
//}

MDOUBLE chrNumberMng::computeExpectationOfChangePerBranch(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost, tree::nodeP node, int fromState, int toState)
{
	int alphabetSize = _pSp->alphabetSize();
	MDOUBLE nodeExpectation = 0;
	for (int x = 0; x<alphabetSize; ++x)
	{
		for (int y = 0; y<alphabetSize; ++y)
		{
			nodeExpectation += (jointPost[node->id()][x][y]* sim.getExpectation(node->name(), x, y, fromState, toState));
		}
	}
	LOGnOUT(10,<<"node:"<<node->name()<<" from,to="<<_pAlph->id2Count(fromState)<<","<<_pAlph->id2Count(toState)<<" exp="<<nodeExpectation<<endl);
	return nodeExpectation;
}


MDOUBLE chrNumberMng::computeProbOfChangePerBranch(simulateJumpsLargeAlphabet &sim, const VVVdouble &jointPost, tree::nodeP node, int fromState, int toState)
{
	int alphabetSize = _pSp->alphabetSize();
	MDOUBLE nodeProb = 0;
	for (int x = 0; x<alphabetSize; ++x)
	{
		for (int y = 0; y<alphabetSize; ++y)
		{
			nodeProb += (jointPost[node->id()][x][y]* sim.getProb(node->name(), x, y, fromState, toState));
		}
	}
	LOGnOUT(10,<<"node:"<<node->name()<<" from,to="<<fromState+1<<","<<toState+1<<" prob="<<nodeProb<<endl);
	return nodeProb;
}


//models: states which model to perform inference on
void chrNumberMng::loopRemoveTaxa()
{
	//should remove taxa from sc and tree
	//print the original LL diff and for each node the LL diff after removal
	printRunInfo();
	getStartingData();
	getTree();
	initStartingModel(_pAlph, chrNumberOptions::_maxBaseTransition);
	sequenceContainer origSc = _sc;
	map<string, int> name2idSc;
	for (int j = 0; j<origSc.numberOfSeqs(); ++j){
		name2idSc[origSc[j].name()] = origSc[j].id();
	}


	//print heading
	string llFileName = chrNumberOptions::_outDir + "//" + "llDif.txt";
	ofstream out(llFileName.c_str());
	out<<"NODE"<<"\t";
	out<<chrNumberOptions::getModelTypeStr(chrNumberOptions::_modelType)<<"\t";
	out<<endl;
	out<<"ALL"<<"\t";
	MDOUBLE ll = optimizeParams();
	out<<ll<<"\t";
	out<<endl;

	string tmpTreeFileName = chrNumberOptions::_outDir + "//" + "tmpTree.txt";
    tree originalTree = _tree;
	vector <tree::nodeP> leaves;
	originalTree.getAllLeaves(leaves, originalTree.getRoot());
	for (int i=0; i< leaves.size(); ++i) 
	{
		//remove leaf
		tree tmpTree(originalTree);
		string leafName = (leaves[i])->name();
		out<<leafName<<"\t";
		tree::nodeP thisNode = tmpTree.findNodeByName(leafName);
		if (thisNode) 
			tmpTree.removeLeaf(thisNode);
		else
            errorMsg::reportError("cannot find: " + leafName);

		tmpTree.output(tmpTreeFileName);
		_tree = tree(tmpTreeFileName);
		sequenceContainer blanko;
		buildScFromTree(_tree, origSc, name2idSc, blanko);
		_sc = blanko;

		MDOUBLE ll = optimizeParams();
		out<<ll<<"\t";
		out<<endl;
	}
}

//scalefactor: speciationalProportion (SP). 0<=SP<1
//In this case the total tree length stays the same after the scaling.
//if SP==0 then the process is unrelated to spection and the original tree is used.
//If SP-->1: (internal) branch lengths are nearly equal. 
//SP = x*internal_branch_num/(origTreeLength + x*internal_branch_num). 
//Thus, x, the amount to add to each branch is 
// origTreeLength*SP/(internal_branch_num*(1-SP))
//after adding x to all branches: have to scale tree to preserve original tree length
//in practice: if only internal branches are scaled then have to preserve only the length of the internal branches.
void chrNumberMng::scaleTree(tree &inTree, MDOUBLE scaleFactor)
{
	vector <tree::nodeP> nodes;
	if (chrNumberOptions::_branchModelType == chrNumModel::COMBINED)
        inTree.getAllNodes(nodes, inTree.getRoot());
	else if (chrNumberOptions::_branchModelType == chrNumModel::COMBINED_INTERNALS)
		inTree.getAllHTUs(nodes, inTree.getRoot());
	else 
		errorMsg::reportError("cannot scale tree when _branchModelType is not COMBINED*");

	MDOUBLE origTreeLength = 0.0;
	int i;
	int nodesNum = 0;
	for (i = 0; i < nodes.size(); ++i)
	{
		if (nodes[i]->isRoot())
			continue;
        origTreeLength += nodes[i]->dis2father();
		nodesNum++;
	}

	MDOUBLE add = scaleFactor * origTreeLength / ((1-scaleFactor) * nodesNum);
	MDOUBLE lenAfterScale = origTreeLength + nodesNum*add;
	for (i = 0; i < nodes.size(); ++i)
	{
		if (nodes[i]->isRoot())
			continue;
		MDOUBLE oldDist = nodes[i]->dis2father();
		//nodes[i]->setDisToFather(oldDist + add);
		MDOUBLE newDist = (oldDist + add)*origTreeLength / lenAfterScale;
		nodes[i]->setDisToFather(newDist);

	}
	//rescale tree to original length
	MDOUBLE newTreeLength = getSumOfBranchLengths(inTree);
	inTree.multipleAllBranchesByFactor(origTreeLength / newTreeLength);
	MDOUBLE finalTreeLength = getSumOfBranchLengths(inTree);
	//cerr<<"treeLength0="<<origTreeLength<<" final = "<<finalTreeLength<<endl;

}


Vdouble chrNumberMng::estimateExpectationsHeuristicly(int nodeID, const VVVdouble& jointPost)
{
	int eventsNum = 5; //should be JUMP_TYPE_MAX
	Vdouble res(eventsNum, 0.0);
	int alphabetSize = _pSp->alphabetSize();

	for (int fromState = 0; fromState<alphabetSize; ++fromState)
	{
		for (int toState = 0; toState<alphabetSize; ++toState)
		{
			Vdouble est = estimateEvents(fromState, toState, _pSp, _pAlph);
			for (int e = 0; e < eventsNum; ++e)
					res[e] += jointPost[nodeID][fromState][toState] * est[e];
		}
	}
	return res;
}

// a stupid heuristic to estimate the type of events in case it couldn't be computed from the simulations
//because of the possibility of a composite id - we have to loop over all possible 
Vdouble chrNumberMng::estimateEvents(int fromId, int toId, const stochasticProcess* pSp, const finiteIntAlphabet* pAlph)
{
	//find the to and from chromosomes (not trivial in case of composite0
	Vint fromChroms, toChroms;
	Vdouble fromProbs, toProbs;
	Vint compIds1, compIds2;
	pAlph->getCompositeParts(fromId, compIds1, fromProbs);
	for (int i = 0; i < compIds1.size(); ++i) 
		fromChroms.push_back(pAlph->id2Count(compIds1[i]));
	pAlph->getCompositeParts(toId, compIds2, toProbs);
	for (int i = 0; i < compIds2.size(); ++i) 
		toChroms.push_back(pAlph->id2Count(compIds2[i]));
	

	MDOUBLE thresholdR = THRESHOLD_R; //rates of events types below threshold will not be considered
	Vdouble res(5, 0.0);
	chrNumModel* pModel = static_cast<chrNumModel*>(pSp->getPijAccelerator()->getReplacementModel());
	for (int f = 0 ; f < fromChroms.size(); ++f) {
		for (int t = 0 ; t < toChroms.size(); ++t) {
			int fromChrom = fromChroms[f];
			int toChrom = toChroms[t];
			MDOUBLE pairProb = toProbs[t] * fromProbs[f];
            int diff = toChrom - fromChrom;
            if (diff == 0)
                continue;
			if (diff < 0) {
				res[chrNumModel::LOSS_J] = -diff * pairProb;
			}
			MDOUBLE gainR, lossR, duplR, demiR, baseR;
			while (fromChrom < toChrom)
			{
				gainR = pModel->getGainR(pAlph->count2Id(fromChrom));
				lossR = pModel->getLossR(pAlph->count2Id(fromChrom));
				duplR = pModel->getDuplR(pAlph->count2Id(fromChrom));
				demiR = pModel->getDemiR(pAlph->count2Id(fromChrom));
				baseR = pModel->getBaseNumberR();

				//divide all rates by the max rate so that the most probable transition type is the one with the max rate 
				//otherwise - if the rate is very high then the probability of one transition is low since there is a higher prob for a higher number of transitions
				MDOUBLE maxR = max(max(baseR, max(gainR, lossR)), max(duplR, demiR));
				gainR/= maxR;
				lossR/= maxR;
				duplR/= maxR;
				demiR/= maxR;
				baseR/= maxR;
				thresholdR *=maxR;

				MDOUBLE ratio = toChrom/MDOUBLE(fromChrom);
				MDOUBLE baseProb = -1;
				int gainNumberB = -1;
				int lossNumberB = -1;
				if (baseR > thresholdR) {
					//a baseNumber move can either go (1) directly to toChrom, (2) to a number above toChrom (in which case we'll need additional loses) or (3) to a number below toChrom (in which case we'll need additional gains)
					int diff = toChrom - fromChrom;
					int baseNumber = pModel->getBaseNumber(); 
					int multiplier = (int) diff / baseNumber;
					int reminder = diff % baseNumber;
					MDOUBLE base_loss_prob = 0.0; //the prob of base move + losses
					MDOUBLE base_gain_prob = 0.0;//the prob of base move + gains
					gainNumberB = reminder;
					lossNumberB = baseNumber - reminder;
					if (reminder > 0) {//need for losses
						base_loss_prob = copmutePoissonProbability(1, baseR) * copmutePoissonProbability(lossNumberB, lossR);
						base_gain_prob = copmutePoissonProbability(1, baseR) * copmutePoissonProbability(gainNumberB, gainR);
					}
					else {
						base_gain_prob = copmutePoissonProbability(1, baseR);
					}
					baseProb = max(base_gain_prob, base_loss_prob);
					if (base_gain_prob > base_loss_prob)
						lossNumberB = 0;
					else 
						gainNumberB = 0;
				}
				if ((ratio >= 2) && (duplR > thresholdR))
				{	//here we should first check if demi_PP is in the model. 
					//If yes - check the prob for a single duplications (don't worry about additional transitions)
					//If no - we can calculate the number of duplications (+ gains or losses) needed 
					if (demiR <= thresholdR) {
						int lowerDuplNum = log(ratio)/log(2.0); //calculate log_base_2(ratio)
						int gainNumber = toChrom - fromChrom * pow(2.0, lowerDuplNum);
						int lossNumber = fromChrom * pow(2.0,lowerDuplNum+1) - toChrom;
						MDOUBLE duplGainProb = copmutePoissonProbability(lowerDuplNum, duplR) * copmutePoissonProbability(gainNumber, gainR);
						MDOUBLE duplLossProb = copmutePoissonProbability(lowerDuplNum+1, duplR) * copmutePoissonProbability(lossNumber, lossR);
						
						if ((baseProb > duplGainProb) && (baseProb > duplLossProb)) {
							res[chrNumModel::BASE_J] += (1 * pairProb);	
							res[chrNumModel::GAIN_J] += (gainNumberB*pairProb);
							res[chrNumModel::LOSS_J] += (lossNumberB*pairProb);
							fromChrom = toChrom;
							continue;
						}
						else if (duplGainProb > duplLossProb) {
							res[chrNumModel::DUPL_J] += (lowerDuplNum * pairProb);	
							res[chrNumModel::GAIN_J] += (gainNumber*pairProb);
							fromChrom = toChrom;
							continue;
						}
						else { // the best is duplication + losses
							res[chrNumModel::DUPL_J] += ((1+lowerDuplNum) * pairProb);	
							res[chrNumModel::LOSS_J] += (lossNumber*pairProb);
						}

					} // turning to the case where there is a demi model - so we go one duplicaiton at a time
					MDOUBLE duplProb = copmutePoissonProbability(1, duplR); //assuming only one duplication
					if (duplProb > baseProb) {
						res[chrNumModel::DUPL_J] += (1 * pairProb);		
						fromChrom *=2;
						continue;
					}
					else { //base_number move
						res[chrNumModel::BASE_J] += (1 * pairProb);	
						res[chrNumModel::GAIN_J] += (gainNumberB*pairProb);
						res[chrNumModel::LOSS_J] += (lossNumberB*pairProb);
						fromChrom = toChrom;
					}
				}
				else if ((ratio >= 1.5) && ((duplR > thresholdR) || (demiR > thresholdR)))
				{ //check if there was higher probability for demi+gains or dupl+losses
					int lossNumber =  ((toChrom % 2)==0)? (2*fromChrom - toChrom)/2 : (1+(2*fromChrom - (toChrom+1))/2); //first loss as many as you can and then dupl
					if (duplR <= thresholdR)
						lossNumber = 0;
					//int gainNumber = toChrom - static_cast<int>(1.5*fromChrom);
					int gainNumber = static_cast<int>(2.0/3.0 * (toChrom - toChrom%3) - fromChrom + toChrom%3);  //first gain as many as you can and then demi-dupl. The demi-dupl is always made from an even number
					if (ceil(fromChrom * 1.5) == toChrom)
						gainNumber = 0;
					MDOUBLE duplProb = copmutePoissonProbability(1, duplR) * copmutePoissonProbability(lossNumber, lossR);
					MDOUBLE demiProb = copmutePoissonProbability(1, demiR) * copmutePoissonProbability(gainNumber, gainR);
					if (baseProb > max(duplProb, demiProb)) {
						res[chrNumModel::BASE_J] += (1 * pairProb);	
						res[chrNumModel::GAIN_J] += (gainNumberB*pairProb);
						res[chrNumModel::LOSS_J] += (lossNumberB*pairProb);
						fromChrom = toChrom;
						continue;
					}
		            else if (duplProb >= demiProb)
					{
						res[chrNumModel::DUPL_J] +=(1*pairProb);
						res[chrNumModel::LOSS_J] += (lossNumber*pairProb);
						fromChrom = toChrom;
						continue;
					}
					else
					{
						res[chrNumModel::DEMI_J] +=(1*pairProb);
						res[chrNumModel::GAIN_J] += (gainNumber*pairProb);
						fromChrom = toChrom;
						continue;
					}
				}
				else 
				{	//should be here if the ratio is less than 1.5
					//decide between demi-losses, dupl-losses, gains
					int lossNumberForDupl = ((toChrom % 2)==0)? (2*fromChrom - toChrom)/2 : (1+(2*fromChrom - (toChrom+1))/2); //first loss as many as you can and then dupl
					int lossNumberForDemi = fromChrom - 2*static_cast<int>((toChrom+2)/3);
					if (floor(fromChrom * 1.5) == toChrom)
						lossNumberForDemi = 0;
					int gainNumber = toChrom - fromChrom;
					MDOUBLE duplProb = 0.0;
					MDOUBLE demiProb = 0.0;
					MDOUBLE gainProb = 0.0;
					if ((lossNumberForDupl > 0) && (duplR > thresholdR))
						duplProb = copmutePoissonProbability(1, duplR) * copmutePoissonProbability(lossNumberForDupl, lossR);
					if ((lossNumberForDemi >= 0) && (demiR > thresholdR)) {
						if (floor(fromChrom * 1.5) == toChrom) {
							lossNumberForDemi = 0; }
						else if (toChrom % 3 ==1) {
							lossNumberForDemi +=2;}
						else if (toChrom % 3 ==2) 
							lossNumberForDemi +=1;
						if (lossNumberForDemi < 0) 
							LOGnOUT(3, <<"error in  chrNumberMng::estimateEvents. lossNumberForDemi < 0 fromChrom = "<<fromChrom<<" toChrom="<<toChrom);
						demiProb = copmutePoissonProbability(1, demiR) * copmutePoissonProbability(lossNumberForDemi, lossR);
					}
					if ((gainNumber > 0) && (gainR > thresholdR))
                        gainProb = copmutePoissonProbability(gainNumber, gainR);
					MDOUBLE maxProb = max(max(baseProb, gainProb), max(duplProb, demiProb));
					if (maxProb == 0)
						LOGnOUT(6,<<"the probability of gains, demi-PP and duplications is zero in chrNumberMng::estimateEvents("<<int2string(fromChrom)<<","<<int2string(toChrom)<<")"<<endl);
					if (maxProb == gainProb)
					{
						res[chrNumModel::GAIN_J] += (gainNumber*pairProb);
						fromChrom += gainNumber;
						continue;
					}
					else if (maxProb == duplProb)
					{
						res[chrNumModel::DUPL_J] +=(1*pairProb);
						res[chrNumModel::LOSS_J] += (lossNumberForDupl*pairProb);
						fromChrom = toChrom;
						continue;
					}
					else if (maxProb == demiProb)
					{
						res[chrNumModel::DEMI_J] +=(1*pairProb);
						res[chrNumModel::LOSS_J] += (lossNumberForDemi*pairProb);
						fromChrom = toChrom;
						continue;
					}
					else if (maxProb == baseProb) 
					{
						res[chrNumModel::BASE_J] += (1 * pairProb);	
						res[chrNumModel::GAIN_J] += (gainNumberB*pairProb);
						res[chrNumModel::LOSS_J] += (lossNumberB*pairProb);
						fromChrom = toChrom;
						continue;
					}
				}
			}
		}
	}
	return res;
}

void chrNumberMng::computeExpectationsFromRootToNodes(VVdouble& nodesTransitions)
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
		nodesTransitions[nodeId][chrNumModel::GAIN_J] = nodesTransitions[fatherId][chrNumModel::GAIN_J] + _expChanges[nodeId][chrNumModel::GAIN_J];
		nodesTransitions[nodeId][chrNumModel::LOSS_J] = nodesTransitions[fatherId][chrNumModel::LOSS_J] + _expChanges[nodeId][chrNumModel::LOSS_J];
		nodesTransitions[nodeId][chrNumModel::DUPL_J] = nodesTransitions[fatherId][chrNumModel::DUPL_J] + _expChanges[nodeId][chrNumModel::DUPL_J];
		nodesTransitions[nodeId][chrNumModel::DEMI_J] = nodesTransitions[fatherId][chrNumModel::DEMI_J] + _expChanges[nodeId][chrNumModel::DEMI_J];
		nodesTransitions[nodeId][chrNumModel::BASE_J] = nodesTransitions[fatherId][chrNumModel::BASE_J] + _expChanges[nodeId][chrNumModel::BASE_J];
		nodesTransitions[nodeId][chrNumModel::MAX_CHR_J] = nodesTransitions[fatherId][chrNumModel::MAX_CHR_J] + _expChanges[nodeId][chrNumModel::MAX_CHR_J];
	}
}
