#include "likelihoodComputation.h"
#include "stochasticProcess.h"
#include "pijAccelerator.h"
#include "tree.h"
#include "someUtil.h"
#include "treeUtil.h"
#include "recognizeFormat.h"
#include "chrNumberMng.h"
#include "chrNumberOptions.h"
#include "chrNumModel.h"
#include "talRandom.h"
#include "matrixUtils.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "dataUtils.h"
#include "chrCountFormat.h"
#include "fastaFormat.h"
#include "amino.h"
#include "codon.h"
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;



void mainSimulate();
void mainDataPrep();
void translateToPhylip(int argc, char **argv);
void translatePhylipTree(int argc, char **argv);
void convertNexusFormat(int argc, char **argv);
void removeTaxa();
void getTreeHeight(int argc, char **argv);
//void concatenateAlignments(int argc, char **argv);
void mainRunAllModels();
void mainRunAllModelsScaleBranch();
void mainCompareLLAfterRemoveLeaves();
void mainGenerateLRTDistribution();
void evaluatePPDist();
void simulatePP();


int main(int argc, char **argv) {
	//talRandom::setSeed(2);
	
	if (argc < 2)
		errorMsg::reportError("No input arguements. Please specify path to parameter file");
	string paramStr = argv[1];
	chrNumberOptions::initOptions(paramStr);
	

	if (chrNumberOptions::_mainType == "mainSimulate")
	{
		mainSimulate();
		return 0;
	}
	else if (chrNumberOptions::_mainType == "mainValidateData")
	{
		chrNumberMng mng;
		mng.validateData();
		return 0;
	}

	else if (chrNumberOptions::_mainType == "All_Models")
	{
		mainRunAllModels();
		return 0;
	}
	else if (chrNumberOptions::_mainType == "All_Models_Scale_Branch")
	{
		mainRunAllModelsScaleBranch();
		return 0;
	}
	else if (chrNumberOptions::_mainType == "mainGenerateLRTDistribution")
	{
		mainGenerateLRTDistribution();
		return 0;
	}
	else if (chrNumberOptions::_mainType == "evaluatePPDist")
	{
		evaluatePPDist();
		return 0;
	}
	else if (chrNumberOptions::_mainType == "simulatePPInPopulation") {
        	simulatePP();
		return 0;
	}
	else if (chrNumberOptions::_mainType == "removeTaxa") {
		removeTaxa();
		return 1;
	}

	chrNumberMng mng;
	mng.run();
	return 0;

}

void mainRunAllModelsScaleBranch()
{
	chrNumberOptions::_mainType = "Optimize_Model";
	string baseOutDirName = chrNumberOptions::_outDir;
	string statsFileName = baseOutDirName + "//models_summary.txt";
	vector<chrNumModel::branchModelType> models;
	models.push_back(chrNumModel::GRADUAL);
	models.push_back(chrNumModel::SPECIATIONAL);
	models.push_back(chrNumModel::COMBINED);
	models.push_back(chrNumModel::COMBINED_INTERNALS);
	

	ofstream statsFile(statsFileName.c_str());
	statsFile.precision(4);
	statsFile<<"MODEL"<<"\t"<<"Log-likelihood"<<"\t"<<"AIC"<<endl;

	for (int m = 0; m < models.size(); ++m)
	{
		string modelName = chrNumberOptions::getBranchModelTypeStr(models[m]);
		chrNumberOptions::_branchModelType = models[m];
		if ((models[m] == chrNumModel::GRADUAL) || (models[m] == chrNumModel::SPECIATIONAL))
			chrNumberOptions::_scaleBranch = IGNORE_PARAM;
		else
			chrNumberOptions::_scaleBranch = 0.5;
		//run chromEvol
		cerr<<endl<<endl<<"START OF OPTIMIZATION: "<<modelName<<endl<<endl;
		string modelOutDirName = baseOutDirName + "//" + modelName;
		chrNumberOptions::_outDir = modelOutDirName;
		createDir("", modelOutDirName);
		chrNumberMng mng;
		mng.run();
		MDOUBLE ll = mng.getBestLL();
		const chrNumModel* pModel = mng.getModel();
		vector<paramObj> paramVec;
		pModel->getParams(paramVec);
		MDOUBLE aic = (-2 * ll) + 2* paramVec.size();
		statsFile<<modelName<<"\t"<<ll<<"\t"<<aic<<endl;
	}
	statsFile.close();
}

void mainRunAllModels()
{
	chrNumberOptions::_mainType = "Optimize_Model";
	string baseOutDirName = chrNumberOptions::_outDir;
	createDir("", chrNumberOptions::_outDir);
	string statsFileName = baseOutDirName + "//models_summary.txt";
	vector<chrNumModel::modelType> models;
	models.push_back(chrNumModel::CONST_RATE_TRUNCATE);
	models.push_back(chrNumModel::CONST_RATE_NO_DUPL);
	models.push_back(chrNumModel::LINEAR_RATE);
	models.push_back(chrNumModel::LINEAR_RATE_NO_DUPL);
	if (chrNumberOptions::_baseNumber != IGNORE_PARAM) {
		models.push_back(chrNumModel::CONST_RATE_BASE_NUMBER);
		models.push_back(chrNumModel::BASE_NUMBER_NO_DUPL);
	}
	

	//vector<chrNumModel::rootFreqType> freqTypes;
	//freqTypes.push_back(chrNumModel::ROOT_LL);
	Vdouble halfDuplType;
	halfDuplType.push_back(IGNORE_PARAM);
	halfDuplType.push_back(DEMI_EQUAL_DUPL);
	halfDuplType.push_back(1);

	int origBase = chrNumberOptions::_baseNumber;
	ofstream statsFile(statsFileName.c_str());
	statsFile.precision(4);
	statsFile<<"MODEL"<<"\t"<<"Log-likelihood"<<"\t"<<"AIC"<<endl;

	for (int m = 0; m < models.size(); ++m)
	{
		for (int f = 0; f < halfDuplType.size(); ++f)
		{
			string fullModelName = chrNumberOptions::getModelTypeStr(models[m]);
			chrNumberOptions::_duplConstR= IGNORE_PARAM;
			chrNumberOptions::_baseNumberR = IGNORE_PARAM;
			chrNumberOptions::_baseNumber = IGNORE_PARAM;
			switch (models[m])
			{
			case chrNumModel::CONST_RATE_TRUNCATE:
				chrNumberOptions::_duplConstR = 0.5;
			case chrNumModel::CONST_RATE_NO_DUPL:
				chrNumberOptions::_gainConstR = 0.5;
				chrNumberOptions::_lossConstR = 0.5;
				chrNumberOptions::_lossLinearR = IGNORE_PARAM; 
				chrNumberOptions::_gainLinearR = IGNORE_PARAM;
				break;
			case chrNumModel::LINEAR_RATE:
				chrNumberOptions::_duplConstR = 0.5;
			case chrNumModel::LINEAR_RATE_NO_DUPL:
				chrNumberOptions::_gainConstR = 0.5;
				chrNumberOptions::_lossConstR = 0.5;
				chrNumberOptions::_lossLinearR = 0.5; 
				chrNumberOptions::_gainLinearR = 0.5;
				break;
			case chrNumModel::CONST_RATE_BASE_NUMBER:
				chrNumberOptions::_duplConstR = 0.5;
			case chrNumModel::BASE_NUMBER_NO_DUPL:
				chrNumberOptions::_baseNumberR = 0.5;
				chrNumberOptions::_baseNumber = origBase;
				chrNumberOptions::_gainConstR = 0.5;
				chrNumberOptions::_lossConstR = 0.5;
				chrNumberOptions::_lossLinearR = IGNORE_PARAM; 
				chrNumberOptions::_gainLinearR = IGNORE_PARAM;
				break;
			default:
				errorMsg::reportError("unknown model type in mainRunAllModels()");
				break;
			}
			
			//fix the demi-ploidy params
            chrNumberOptions::_demiPloidyR = halfDuplType[f];
			string halfDuplName = "";
			if (halfDuplType[f] == DEMI_EQUAL_DUPL)
				halfDuplName = "_DEMI";
			else if (halfDuplType[f] > 0)
				halfDuplName = "_DEMI_EST";
			fullModelName += halfDuplName;
            if ((models[m] == chrNumModel::LINEAR_RATE_NO_DUPL) || (models[m] == chrNumModel::CONST_RATE_NO_DUPL) || (models[m] == chrNumModel::CONST_RATE_BASE_NUMBER) || (models[m] == chrNumModel::BASE_NUMBER_NO_DUPL))
			{ // should perform the loop over demi-ploidy types only once: for the IGNORE_PARAM
				if (halfDuplType[f] != IGNORE_PARAM)
					errorMsg::reportError("cannot run model with no polyploidy but with demi-ploidy");
				else
                    f +=100; //perform the loop only once
			}
			
			//run chromEvol
			cerr<<endl<<endl<<"START OF OPTIMIZATION: "<<fullModelName<<endl<<endl;
			string modelOutDirName = baseOutDirName + "//" + fullModelName;
			chrNumberOptions::_outDir = modelOutDirName;
			chrNumberOptions::_modelType = chrNumModel::GENERAL_CHR_MODEL;
			createDir("", modelOutDirName);
			chrNumberMng mng;
			mng.run();
			MDOUBLE ll = mng.getBestLL();
			const chrNumModel* pModel = mng.getModel();
			vector<paramObj> paramVec;
			pModel->getParams(paramVec);
			MDOUBLE aic = (-2 * ll) + 2* paramVec.size();
			statsFile<<fullModelName<<"\t"<<ll<<"\t"<<aic<<endl;
		}
	}
	statsFile.close();
}



//gets a fasta MSA file. (1) shorten names to less than 10 characters (2) translate to phylip format
//outputs:
void translateToPhylip(int argc, char **argv)
{
	string inMsa = argv[1];
	string outPhylip= argv[2];
	string outConvertFile =  argv[3];
	convertToPhylipFormat(inMsa, outPhylip, outConvertFile);
}
void translatePhylipTree(int argc, char **argv)
{
	string inTreeFileName = argv[1];
	string convertFile = argv[2];
	string outTranslateTree =  argv[3];
	translatePhylipTree(inTreeFileName, convertFile, outTranslateTree);
}



void removeTaxa()
{
	string inTreeFileName = chrNumberOptions::_treeFile;
	string inMissingListFileName =  chrNumberOptions::_dataFile;
	string outTreeFileName =  chrNumberOptions::_outFile;
	ofstream outTreesFile(outTreeFileName.c_str());
	vector<tree> allTrees = getStartingTreeVecFromFile(inTreeFileName);
	for (int t = 0; t < allTrees.size(); ++t) {
		removeTaxa(allTrees[t], inMissingListFileName);
		allTrees[t].output(outTreesFile);
	}
	outTreesFile.close();
}


void convertNexusFormat(int argc, char **argv)
{
	string inFileName = argv[1];
	string outFileName =  argv[2];

	nucleotide ab;
	ifstream inFile(inFileName.c_str());
	sequenceContainer sc = recognizeFormat::read(inFile, &ab);
	inFile.close();
	ofstream outFile(outFileName.c_str());
	fastaFormat::write(outFile, sc);
	outFile.close();
}

void mainDataPrep()
{
	string inMsa = "Data/Agrodiaetus/142_gbLessStringent.fas";
	string outPhylip= "Data/Agrodiaetus/142_gbLessStringent.phy";
	string outConvertFile =  "Data/Agrodiaetus/phyml/142_gbLessStringent_translateTable.txt";
	string inTreeFileName =  "Data/Agrodiaetus/phyml/142_gbLessStringent_phy_phyml_tree.txt";
	string outTranslateTree =  "Data/Agrodiaetus/phyml/142_gbLessStringent_phymlTranslate.tree";
	string inNexusTreeFileName =  "Data/SimulateTrees/30Taxa/100Trees_30Taxa.nex";
	string outDir =  "Data/SimulateTrees/30Taxa/";

	string inNumberedTreeFileName =  "Data/Passiflora/katie1.numbered.tree";
	string outNewickTree =  "Data/Passiflora/12cp_bayes.tree";
	//string outTreeFileName = chrNumberOptions::_outDir + "//" + "outTree.phr";
	string missingListFileName = "Data/helianthus/removeTaxa.txt";
    string logFileName = chrNumberOptions::_outDir + "//" + chrNumberOptions::_logFile;
	myLog::setLog(logFileName, chrNumberOptions::_logValue);
	
	string inTranslateFileName = "Data/Aristolochiaceae/fromTsue/matK_fromWord.txt";
	string numberedTree = "Data/Aristolochiaceae/fromTsue/matK_106_ML_numbers.tree";
	
	
	//string tre1 =  "Data/SimulateTrees/20Taxa/0.tree";
	//string tre2 =  "Data/SimulateTrees/20Taxa/0_noLength.tree";
	//tree tr2(tre2);

	//tree tr1(tre1);
	
	//seperateNexusTrees(inNexusTreeFileName, outDir);
	//return;



	string countsFile = "Data/Aristolochiaceae/fromTsue/final/matK_106_counts.fasta";
	
	string outTreeFileName = "Data/Primula/phyml/primula_filterNoOg.tree";
	string missingFile = "Data/Primula/phyml/removeOutgroup.txt";
	string treeFile =  "Data/Primula/phyml/primula_filter.tree";
	tree inTree(treeFile);

	removeTaxa(inTree, missingFile);
	inTree.output(outTreeFileName);

	//LOG(chrNumberOptions::_logValue, <<"tree rooted at "<<myTree.getRoot()->name()<<" id, "<<myTree.getRoot()->id()<<endl);
	//LOG(chrNumberOptions::_logValue, <<"sons of root are: "<<endl);
	//for (int s = 0; s < myTree.getRoot()->getNumberOfSons(); ++s)
	//{
 //       LOG(chrNumberOptions::_logValue, <<myTree.getRoot()->getSon(s)->name()<<endl;);
	//}

	//rootTree(myTree, chrNumberOptions::_rootAt);
	//removeTaxa(myTree, missingListFileName);
	//myTree.output(outTreeFileName, tree::PHYLIP, true);
	//myTree.output(outTreeFileName);
	//convertNexusTreeToNewickFormat(inNexusTreeFileName, outNewickTree, inNumberedTreeFileName);
	//translatePhylipTree(inTreeFileName, outConvertFile, outTranslateTree);
	//convertToPhylipFormat(inMsa, outPhylip, outConvertFile);

}


void mainSimulate()
{
	chrNumberMng mng;
    mng.runSimulations();
}

void mainCompareLLAfterRemoveLeaves()
{
	//chrNumberOptions::_mainType = "Run_Fix_Param";
	chrNumberMng mng;

	vector<chrNumModel::modelType> models;
	models.push_back(chrNumModel::CONST_RATE_TRUNCATE);
	//models.push_back(chrNumModel::CONST_RATE_EPS);
	mng.loopRemoveTaxa();
}


//if needed see Whealan & Goldman: Mol. Biol. Evol. 16(9):1292–1299. 1999 for details 
void mainGenerateLRTDistribution()
{
	time_t ltime;
	time(&ltime);
	long seed = static_cast<long>(ltime);
	talRandom::setSeed(seed);
	cout<<"the seed = " <<seed<<endl;
	
	chrNumberOptions::_mainType = "Optimize_Model";
	
	string baseOutDirName = chrNumberOptions::_outDir;
	string statsFileName = baseOutDirName + "//stats.txt";
	vector<chrNumModel::modelType> models =  chrNumberOptions::_simModels;
	if (models.size() != chrNumberOptions::_simDemiTypes.size())
		errorMsg::reportError("each model type should also have the demi model type");
	//print headings
	ofstream statsFile(statsFileName.c_str());
	statsFile.precision(6);
	statsFile<<"ITER"<<"\t"<<"SIM"<<"\t";

	int m = 0;
	for (m = 0; m < models.size(); ++m)
	{
		string fullModelName = chrNumberOptions::getModelTypeStr(models[m]);
		if (chrNumberOptions::_simDemiTypes[m] == DEMI_EQUAL_DUPL)
			fullModelName += "_demi";
        else if (chrNumberOptions::_simDemiTypes[m] > 0)
			fullModelName += "_EstDemi";
		statsFile<<fullModelName<<"\t";
	}
	statsFile<<endl;
	string originalFreqFile = chrNumberOptions::_freqFile;
	string originalDataFile = chrNumberOptions::_dataFile;	
	for (int iter = 0; iter < chrNumberOptions::_simulationsIter; ++iter)
	{
		cerr<<"simulation iter = "<<iter<<endl;
		//create dir for the iteration
		string iterOutDirName = baseOutDirName + "//" + int2string(iter);
		chrNumberOptions::_outDir = iterOutDirName;
		createDir("", iterOutDirName);

		//run simulation with first model
		chrNumberOptions::_dataFile = originalDataFile;
		chrNumberOptions::_freqFile = originalFreqFile;
		chrNumberOptions::_outDir = iterOutDirName + "//sim";
        chrNumberOptions::_modelType = models[0];
		chrNumberOptions::_demiPloidyR = static_cast<MDOUBLE>(chrNumberOptions::_simDemiTypes[0]);
        chrNumberMng mng;
        mng.runOneSimulation();
        MDOUBLE ll = mng.getBestLL();
		statsFile<<iter<<"\t"<<ll<<"\t";
		string simScFileName = chrNumberOptions::_outDir + "//" + "simCounts.txt";
		chrNumberOptions::_dataFile = simScFileName;
		chrNumberOptions::_freqFile = "";
        //run inference with first model and all other models
        for (int m = 0; m < models.size(); ++m)
		{
			string fullModelName = chrNumberOptions::getModelTypeStr(models[m]);
			if (chrNumberOptions::_simDemiTypes[m] == DEMI_EQUAL_DUPL)
				fullModelName += "_demi";
			else if (chrNumberOptions::_simDemiTypes[m] > 0)
				fullModelName += "_EstDemi";

			string modelOutDirName = iterOutDirName + "//" + fullModelName;
			chrNumberOptions::_outDir = modelOutDirName;
			chrNumberOptions::_modelType = models[m];
			chrNumberOptions::_demiPloidyR = static_cast<MDOUBLE>(chrNumberOptions::_simDemiTypes[m]);
			createDir("", modelOutDirName);
			chrNumberMng mng;
			mng.run();
			MDOUBLE ll = mng.getBestLL();
			statsFile<<ll<<"\t";
		}
		statsFile<<endl;
	}
	statsFile.close();
}

#include "evalPPDistMng.h"
//the expectations file name is given in chrNumberOptions::_dataFile
void evaluatePPDist() {
	evalPPDistMng mng;
	mng.run();
}

#include "simulatePPInPopulation.h"
#include <algorithm>
void simulatePP() {
	MDOUBLE specRate = chrNumberOptions::_gainConstR;
	MDOUBLE extRate = chrNumberOptions::_lossConstR;
	MDOUBLE ppSpec = chrNumberOptions::_duplConstR;
	int startDiploidsNum = chrNumberOptions::_minChrNum;
	int finalTaxaNum = chrNumberOptions::_maxChrNum;
	string diversificationDist = chrNumberOptions::_dataFile;
	
	VVint allRes(chrNumberOptions::_simulationsNum);
	Vdouble avgs;
	for (int iter = 0; iter < chrNumberOptions::_simulationsNum; ++iter)
	{
        simulatePPInPopulation mng(specRate, extRate, ppSpec, finalTaxaNum);
        string outFile = chrNumberOptions::_outDir + "/sim" + int2string(iter) + ".txt";
        mng.run(startDiploidsNum, diversificationDist, outFile);
		allRes[iter] = mng._ppDist; 
		avgs.push_back(mng._avg);
	}
	ofstream outF(chrNumberOptions::_outFile.c_str());
	outF<<"#speciation rate = "<<specRate<<endl;
	outF<<"#extinction rate = "<<extRate<<endl;
	outF<<"#polyploid speciation = "<<ppSpec<<endl;
	outF<<"#taxa num = "<<finalTaxaNum<<endl;
	outF<<"#start diploids num = "<<startDiploidsNum<<endl;
	outF<<"#Avg number of PP events = "<<computeAverage(avgs)<<endl;
	outF<<"#Std number of PP events = "<<computeStd(avgs)<<endl;
	outF<<"#####"<<endl;
	//transpose matrix
	VVint matT;
    int simNum = allRes.size();
	int levels = allRes[0].size();
	resizeMatrix(matT, levels, simNum);
	for (int i=0; i< levels; i++){
		for (int j=0; j<simNum;j++) {
			matT[i][j]=allRes[j][i];
		}
	}
	outF<<"#pp_level"<<"\t"<<"AVG"<<"\t"<<"STD"<<endl;
	for (int pp = 0; pp < matT.size(); ++pp) {
		outF<<pp+1<<"\t"<<computeAverage(matT[pp])<<"\t"<<computeStd(matT[pp])<<endl;
	}
	outF.close();
}



//void concatenateAlignments(int argc, char **argv)
//{
//	if (argc < 5)
//		errorMsg::reportError("concatenateMsa USAGE: <outFile> <alphabetSize (4, 20, or 61)> <MSA_FILE_1> <MSA_FILE_2> ...");
//	int msaNum = argc-3;
//	string outFile = argv[1];
//	int alphSize = atoi(argv[2]);
//	vector<string> msaNames;
//	for (int i = 0; i < msaNum; ++i)
//		msaNames.push_back(argv[i+3]);
//	alphabet* pAlph;
//	if (alphSize == 4)
//		pAlph = new nucleotide();
//	else if (alphSize == 20)
//		pAlph = new amino();
//	else if (alphSize == 61)
//		pAlph = new codon();
//
//	sequenceContainer outSc = concatenateMsa(msaNames, pAlph);
//	ofstream outF(outFile.c_str());
//	fastaFormat::write(outF, outSc);
//	outF.close();
//	delete pAlph;
//}

#include "treeIt.h"
void getTreeHeight(int argc, char **argv)
{
	string inList = argv[1];
	string inDir= argv[2];
	string outFile= argv[3];
	string treeExt = argv[4];
	ifstream inListF(inList.c_str());
	ofstream outF(outFile.c_str());
	outF<<"ds"<<"\t"<<"Height"<<endl;
	vector<string> inFileData;
	putFileIntoVectorStringArray(inListF, inFileData);
	if (inFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in codonUtility::readCodonUsage");
	}
	vector<string>::const_iterator it1;
	for (it1 = inFileData.begin(); it1!= inFileData.end(); ++it1) {
		if (it1->empty()) 
			continue; // empty line continue
		string ds = *it1;
		string treeFileName = inDir + "//" + ds + treeExt;
		tree inTree(treeFileName);
		Vdouble heights(inTree.getNodesNum(), -1);
		treeIterDownTopConst tit(inTree);
		for (tree::nodeP mynode = tit.first(); mynode != tit.end(); mynode = tit.next()) {
			if (mynode->isLeaf()) {
				heights[mynode->id()] = 0.0;
				continue;
			}
			MDOUBLE sumHeights = 0.0;
			for (int i = 0; i < mynode->getNumberOfSons(); ++i) {
				tree::nodeP sonP = mynode->getSon(i);
				MDOUBLE sonHeight = heights[sonP->id()];
				if (sonHeight < 0)
					errorMsg::reportError("height of son" + sonP->name() + "must be positive");

				sumHeights += sonHeight + sonP->dis2father();
			}
			MDOUBLE nodeHeight = sumHeights / mynode->getNumberOfSons();
			heights[mynode->id()] = nodeHeight;
		}
		MDOUBLE treeHeight = heights[inTree.getRoot()->id()];
		outF<<ds<<"\t"<<treeHeight<<endl;


	}

	inListF.close();
	outF.close();
}
