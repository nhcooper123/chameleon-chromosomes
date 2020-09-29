#ifndef __CHR_NUMBER_OPTION
#define __CHR_NUMBER_OPTION

#include "definitions.h"
#include "chrNumModel.h"
#include <string>
#include <fstream>
using namespace std;

class chrNumberOptions {
		
public:
    static void initOptions(const string& paramFileName);
	virtual ~chrNumberOptions ();

	static string getModelTypeStr(chrNumModel::modelType type);
    static chrNumModel::modelType getModelTypeFromStr(const string& str);
	static vector<chrNumModel::modelType> getModelsFromStr(const string& inStr);
	static string getStrFromModels(const vector<chrNumModel::modelType>& inVec);

	static string getRootFreqTypeStr(chrNumModel::rootFreqType type);
    static chrNumModel::rootFreqType getRootFreqTypeFromStr(const string& str);
	static string getBranchModelTypeStr(chrNumModel::branchModelType type);
    static chrNumModel::branchModelType getBranchModelTypeFromStr(const string& str);

	static void getInitParams(vector<paramObj>& initParams);
private:
	static void initDefault();
	static void getParamsFromFile(const string& paramFileName);

public:
	static string _mainType;
	static string _treeFile;
	static string _dataFile;
	static string _logFile;
    static int _logValue; 
	static string _outDir;
	static string _outFile;
	static string _inferTreeFile;
	static string _freqFile; //input frequencies file. each freq is in a seperat line and looks like: F[1]=0.1
	static chrNumModel::modelType _modelType;
	static chrNumModel::rootFreqType _rootFreqType;
	
	static string _rootAt;
	//the maximal and minimal chromose number allowed . 
	static int _maxChrNum; //if negative: add abs(value) to the maximal count observed in the data.
	static int _minChrNum; //if negative: subtract abs(value) to the minimal count observed in the data.
	static MDOUBLE _branchMul; // if > 0 then multiply all branches by this scalar. if = 999 then total tree length is equal to the number of different character types
	static int _pow2Scale;//when calculating e^(Qt): scale the Q so that its norm*t will be smaller than 2^x
	
	//rate parameters:
	//constant = does not depend on the current number of chromosomes
	static MDOUBLE _gainConstR; 
	static MDOUBLE _gainLinearR; 
	static MDOUBLE _lossConstR; 
	static MDOUBLE _lossLinearR; 
	static MDOUBLE _duplConstR; 
	static MDOUBLE _demiPloidyR; //specify i->1.5i transitions. If ==-1 then ignore. If ==DEMI_EQUAL_DUPL then equal to _duplConstR
	static MDOUBLE _epsR; 
	//specify the basic unit of duplication. so that if _baseNumber=4 then trnasitions from 7 to 11 or 15 are allowed
	static int _baseNumber;
	static MDOUBLE _baseNumberR; 

	//params for root frequencies distribution
	//speciational change model versus gradual
	static chrNumModel::branchModelType _branchModelType;
	static MDOUBLE _scaleBranch; //adds to all branches of the tree a fixed length
	
	//optimization params
	static int _maxOptimizationIterations; //the total number optimization iterations. Will not be changed by users
	static MDOUBLE _epsilonLLimprovement; //if the log-likelihood after optimization is lower than this threshold - then optimize again. Will not be changed by users
	static MDOUBLE _tolParamOptimization; 
	//parameters for optimize many starts
	//the input to these paraqmeters in the param fule is a string of a list seperated by commas ("10,2,1")
	static Vint _optimizePointsNum; //a vector with the number of points to peformed the current cycle of optimization.
	static Vint _optimizeIterNum; //the number of iterations to perform in each cycle.
	static bool _bOptBaseNumber;
	
	//param to calculate the expected number of events
	static int _simulationsNum;//number of simulations to compute the expectation # of changes of certain type along each branch
	static int _smIter;//number of stochastic mapping iterations to be used. Note that the user cannot enter a positive number for both _smIter and _simulationsNum since onkly one kind of simulations should be performed
	//simulation parameters
	static string _simulationsTreeDir;//specify the folder where all the simulated trees are
	static MDOUBLE _simulationsTreeLength;//the total length of the simulated tree. 0.0: don't scale. positive(1.5): indicate total length. Negative(-1.5): the average length from each leaf to the root. We assume the tree is ultrametric
	static string _simulationsJumpsStats;//out statistics file to compare simulated jumps versus inferred expected events. If empty the run regular simulations
	static int _simulationsIter;//number of simulation iterations: this is used to compute the deltaLL distribution or when checking the accuracy
	static int _startSimulationsIter; //when simulating random trees - simulate tree numbers _startSimulationsIter to  _simulationsIter.
	static int _maxChrNumForSimulations; //specify the maximal chromosome number when performing simulations
	static int _maxBaseTransition; //used only when using the base number parameter - the largest base number transition allowed
	static int _minBaseTransition; //used only when using the base number parameter - the minimum base number transition allowed
	static string _baseTransitionProbs; //used only when using the base number parameter - a string representing the possible transition and its possible probabilities. 9=0.4_18=0.6 means that only 2 base transitions are possible for j=i+9 woth probability 0.4 and to j=i+18 with probability 0.6
	static vector<chrNumModel::modelType> _simModels;//_simModels[0] is the null model and all others are the models used for inference
	static Vint _simDemiTypes;//specify for each _simModels which demiPP type to use: DEMI_EQUAL_DUPL= demi=dupl / -1= no demi / 1 = demi estimated
};
#endif
