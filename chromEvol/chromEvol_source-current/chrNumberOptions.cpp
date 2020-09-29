#include "chrNumberOptions.h"
#include "errorMsg.h"
#include "someUtil.h"
#include "Parameters.h"
#include <iostream>
using namespace std;

string chrNumberOptions::_treeFile;
string chrNumberOptions::_dataFile;
string chrNumberOptions::_rootAt;
string chrNumberOptions::_mainType;
string chrNumberOptions::_logFile;
int chrNumberOptions::_logValue;
string chrNumberOptions::_inferTreeFile;
string chrNumberOptions::_outDir;
string chrNumberOptions::_outFile;
string chrNumberOptions::_freqFile;
int chrNumberOptions::_maxChrNum;
int chrNumberOptions::_minChrNum;
int chrNumberOptions::_baseNumber;
bool chrNumberOptions::_bOptBaseNumber; 
MDOUBLE chrNumberOptions::_baseNumberR;
MDOUBLE chrNumberOptions::_gainConstR;
MDOUBLE chrNumberOptions::_gainLinearR;
MDOUBLE chrNumberOptions::_lossConstR; 
MDOUBLE chrNumberOptions::_lossLinearR; 
MDOUBLE chrNumberOptions::_duplConstR; 
MDOUBLE chrNumberOptions::_epsR;
int chrNumberOptions::_maxOptimizationIterations; 
MDOUBLE chrNumberOptions::_epsilonLLimprovement; 
MDOUBLE chrNumberOptions::_tolParamOptimization; 
Vint chrNumberOptions::_optimizePointsNum; 
Vint chrNumberOptions::_optimizeIterNum; 
chrNumModel::modelType chrNumberOptions::_modelType;
chrNumModel::rootFreqType chrNumberOptions::_rootFreqType;
MDOUBLE chrNumberOptions::_demiPloidyR;
int chrNumberOptions::_simulationsNum;
int chrNumberOptions::_smIter;
int chrNumberOptions::_simulationsIter;
int chrNumberOptions::_startSimulationsIter;
vector<chrNumModel::modelType> chrNumberOptions::_simModels;
Vint chrNumberOptions::_simDemiTypes;
MDOUBLE chrNumberOptions::_branchMul;
int chrNumberOptions::_pow2Scale;
string chrNumberOptions::_simulationsTreeDir;
MDOUBLE chrNumberOptions::_simulationsTreeLength;
string chrNumberOptions::_simulationsJumpsStats;
int chrNumberOptions::_maxChrNumForSimulations;
chrNumModel::branchModelType chrNumberOptions::_branchModelType;
MDOUBLE chrNumberOptions::_scaleBranch; 
int chrNumberOptions::_maxBaseTransition;
int chrNumberOptions::_minBaseTransition;
string chrNumberOptions::_baseTransitionProbs;

chrNumberOptions::~chrNumberOptions() 
{}

void chrNumberOptions::initOptions(const string& paramFileName)
{
	initDefault();
	getParamsFromFile(paramFileName);
}

void chrNumberOptions::initDefault()
{
// DEFAULTS VALUES:
    _treeFile = "";
	_dataFile = "";
	_rootAt = "";
    _mainType = "Optimize_Model";
	_logFile = "log.txt";
	_logValue = 5;
	_outDir = "RESULTS";
	_outFile= "chromEvol.res";
	_inferTreeFile = "mlAncestors.tree";
	_freqFile= "";
	_maxChrNum = -10;
	_minChrNum = 1;
	_baseNumber = IGNORE_PARAM;
	_baseNumberR = IGNORE_PARAM;
	_gainConstR = IGNORE_PARAM;
	_gainLinearR = IGNORE_PARAM;
	_lossConstR = IGNORE_PARAM; 
	_lossLinearR = IGNORE_PARAM; 
	_duplConstR = IGNORE_PARAM; 
	_demiPloidyR = IGNORE_PARAM;
	_bOptBaseNumber = false;
	_epsR = IGNORE_PARAM; 
	_maxOptimizationIterations = 5;
	_tolParamOptimization = 0.01; 
	_epsilonLLimprovement = 0.1; 
	_modelType = chrNumModel::GENERAL_CHR_MODEL;
	_rootFreqType = chrNumModel::ROOT_LL;
	_optimizePointsNum.push_back(10); 
	_optimizePointsNum.push_back(3); 
	_optimizePointsNum.push_back(1); 
	_optimizeIterNum.push_back(0); 
	_optimizeIterNum.push_back(2); 
	_optimizeIterNum.push_back(_maxOptimizationIterations); 
	_simulationsNum = 10000;
	_smIter = 0;
	int _simulationsIter = 50;
	int _startSimulationsIter = 0;
	_simModels.push_back(chrNumModel::CONST_RATE_NO_DUPL);
	_simModels.push_back(chrNumModel::CONST_RATE_TRUNCATE);
	_simModels.push_back(chrNumModel::CONST_RATE_TRUNCATE);
	_simDemiTypes.push_back(IGNORE_PARAM);
	_simDemiTypes.push_back(IGNORE_PARAM);
	_simDemiTypes.push_back(DEMI_EQUAL_DUPL);
	_branchMul = 999;
	_simulationsTreeDir = "";
	_simulationsTreeLength = 0.0;
	_pow2Scale = 1;
	_simulationsJumpsStats = "";
	_maxChrNumForSimulations = 0;
	_branchModelType = chrNumModel::GRADUAL;
	_scaleBranch = IGNORE_PARAM;
	_maxBaseTransition = 0;
	_minBaseTransition = 3;
	_baseTransitionProbs = "";
	
	Parameters::addParameter("_treeFile", _treeFile);
	Parameters::addParameter("_dataFile", _dataFile);
	Parameters::addParameter("_rootAt", _rootAt);
	Parameters::addParameter("_mainType", _mainType);
	Parameters::addParameter("_logFile", _logFile);
	Parameters::addParameter("_logValue", _logValue);
	Parameters::addParameter("_outDir", _outDir);
	Parameters::addParameter("_outFile", _outFile);
	Parameters::addParameter("_inferTreeFile", _inferTreeFile);
	Parameters::addParameter("_freqFile", _freqFile);
	Parameters::addParameter("_maxChrNum", _maxChrNum);
	Parameters::addParameter("_minChrNum", _minChrNum);
	Parameters::addParameter("_baseNumber", _baseNumber);
	Parameters::addParameter("_baseNumberR", _baseNumberR);
	Parameters::addParameter("_gainConstR", _gainConstR);
	Parameters::addParameter("_gainLinearR", _gainLinearR);
	Parameters::addParameter("_lossConstR", _lossConstR);
	Parameters::addParameter("_lossLinearR", _lossLinearR);
	Parameters::addParameter("_modelType", getModelTypeStr(_modelType));
	Parameters::addParameter("_rootFreqType", getRootFreqTypeStr(_rootFreqType));
	Parameters::addParameter("_branchModelType", getBranchModelTypeStr(_branchModelType));
	Parameters::addParameter("_duplConstR", _duplConstR);
	Parameters::addParameter("_demiPloidyR", _demiPloidyR);
	Parameters::addParameter("_epsR", _epsR);
	Parameters::addParameter("_maxOptimizationIterations", _maxOptimizationIterations);
	Parameters::addParameter("_tolParamOptimization", _tolParamOptimization); 
	Parameters::addParameter("_epsilonLLimprovement", _epsilonLLimprovement); 
	Parameters::addParameter("_optimizePointsNum", getStrFromVint(_optimizePointsNum)); 
	Parameters::addParameter("_optimizeIterNum", getStrFromVint(_optimizeIterNum)); 
	Parameters::addParameter("_simulationsNum", _simulationsNum);
	Parameters::addParameter("_smIter", _smIter);
	Parameters::addParameter("_simulationsIter", _simulationsIter);
	Parameters::addParameter("_startSimulationsIter", _startSimulationsIter);
	Parameters::addParameter("_simModels", getStrFromModels(_simModels));
	Parameters::addParameter("_simDemiTypes", getStrFromVint(_simDemiTypes)); 
	Parameters::addParameter("_branchMul", _branchMul); 
	Parameters::addParameter("_simulationsTreeDir", _simulationsTreeDir); 
	Parameters::addParameter("_simulationsTreeLength", _simulationsTreeLength); 
	Parameters::addParameter("_pow2Scale", _pow2Scale); 
	Parameters::addParameter("_simulationsJumpsStats", _simulationsJumpsStats); 
	Parameters::addParameter("_maxChrNumForSimulations", _maxChrNumForSimulations); 
	Parameters::addParameter("_scaleBranch", _scaleBranch); 
	Parameters::addParameter("_bOptBaseNumber",(_bOptBaseNumber == true) ? 1 : 0);
	Parameters::addParameter("_maxBaseTransition",_maxBaseTransition);
	Parameters::addParameter("_minBaseTransition",_minBaseTransition);
	Parameters::addParameter("_baseTransitionProbs", _baseTransitionProbs);
}

void chrNumberOptions::getParamsFromFile(const string& paramFileName)
{
	ifstream params(paramFileName.c_str());
	if(params.fail())
        errorMsg::reportError("cannot open parameter file: " + paramFileName);

	if(params.good())
        Parameters::readParameters(params);
	params.close();

	_treeFile = Parameters::getString("_treeFile");
	_dataFile = Parameters::getString("_dataFile");
	_rootAt = Parameters::getString("_rootAt");
	_mainType = Parameters::getString("_mainType");
	_logFile= Parameters::getString("_logFile");
	_logValue = Parameters::getInt("_logValue");
	_outDir = Parameters::getString("_outDir");
	_outFile = Parameters::getString("_outFile");
	_inferTreeFile = Parameters::getString("_inferTreeFile");
	_freqFile = Parameters::getString("_freqFile");
	_maxChrNum = Parameters::getInt("_maxChrNum");
	_minChrNum = Parameters::getInt("_minChrNum");
	_baseNumber = Parameters::getInt("_baseNumber");
	_baseNumberR = Parameters::getFloat("_baseNumberR");
	_gainConstR= Parameters::getFloat("_gainConstR");
	_gainLinearR= Parameters::getFloat("_gainLinearR");
	_lossConstR = Parameters::getFloat("_lossConstR");
	_lossLinearR = Parameters::getFloat("_lossLinearR");
	_duplConstR = Parameters::getFloat("_duplConstR");
	_demiPloidyR= Parameters::getFloat("_demiPloidyR");
	_epsR = Parameters::getFloat("_epsR");
	_modelType = getModelTypeFromStr(Parameters::getString("_modelType"));
	_rootFreqType = getRootFreqTypeFromStr(Parameters::getString("_rootFreqType"));
	_maxOptimizationIterations = Parameters::getInt("_maxOptimizationIterations");
	_optimizePointsNum = getVintFromStr(Parameters::getString("_optimizePointsNum"));
	_optimizeIterNum = getVintFromStr(Parameters::getString("_optimizeIterNum"));
	_epsilonLLimprovement = Parameters::getFloat("_epsilonLLimprovement");
	_simulationsNum = Parameters::getInt("_simulationsNum");
	_smIter = Parameters::getInt("_smIter");
	_simulationsIter = Parameters::getInt("_simulationsIter");
	_startSimulationsIter = Parameters::getInt("_startSimulationsIter");
	_simModels = getModelsFromStr(Parameters::getString("_simModels"));
	_simDemiTypes = getVintFromStr(Parameters::getString("_simDemiTypes"));
	_branchMul = Parameters::getFloat("_branchMul");
	_simulationsTreeDir = Parameters::getString("_simulationsTreeDir");
	_simulationsTreeLength = Parameters::getFloat("_simulationsTreeLength");
	_pow2Scale = Parameters::getInt("_pow2Scale");
	_simulationsJumpsStats = Parameters::getString("_simulationsJumpsStats");	
	_maxChrNumForSimulations = Parameters::getInt("_maxChrNumForSimulations");
	_branchModelType = getBranchModelTypeFromStr(Parameters::getString("_branchModelType"));
	_scaleBranch = Parameters::getFloat("_scaleBranch");
	_bOptBaseNumber = (Parameters::getInt("_bOptBaseNumber") == 1) ? true : false;
	_maxBaseTransition = Parameters::getInt("_maxBaseTransition");
	_minBaseTransition = Parameters::getInt("_minBaseTransition");
	_baseTransitionProbs = Parameters::getString("_baseTransitionProbs");
}

chrNumModel::modelType chrNumberOptions::getModelTypeFromStr(const string& str)
{
	if (str == "CONST_RATE_TRUNCATE")
		return chrNumModel::CONST_RATE_TRUNCATE;
	else if (str == "CONST_RATE")
		return chrNumModel::CONST_RATE_TRUNCATE;
	else if (str == "CONST_RATE_NO_DUPL")
		return chrNumModel::CONST_RATE_NO_DUPL;
	else if (str == "LINEAR_RATE_NO_DUPL")
		return chrNumModel::LINEAR_RATE_NO_DUPL;
	else if (str == "LINEAR_RATE")
		return chrNumModel::LINEAR_RATE;
	else if (str == "BASE_NUMBER")
		return chrNumModel::CONST_RATE_BASE_NUMBER;
	else if (str == "BASE_NUMBER_NO_DUPL")
		return chrNumModel::BASE_NUMBER_NO_DUPL;
	else if (str == "CONST_RATE_EPS")
		return chrNumModel::CONST_RATE_EPS;
	else if (str == "CONST_RATE_EPS_NO_DUPL")
		return chrNumModel::CONST_RATE_EPS_NO_DUPL;
	else if (str == "GENERAL_CHR_MODEL")
		return chrNumModel::GENERAL_CHR_MODEL;
	else
        errorMsg::reportError("unknown type in chrNumberOptions::getModelTypeFromStr");
	return chrNumModel::CONST_RATE_TRUNCATE;
}


string chrNumberOptions::getModelTypeStr(chrNumModel::modelType type)
{
	string res = "";
	switch (type)
	{
	case chrNumModel::CONST_RATE_TRUNCATE:
		res = "CONST_RATE";
		break;
	case chrNumModel::CONST_RATE_NO_DUPL:
		res = "CONST_RATE_NO_DUPL";
		break;
	case chrNumModel::LINEAR_RATE_NO_DUPL:
		res = "LINEAR_RATE_NO_DUPL";
		break;
	case chrNumModel::LINEAR_RATE:
		res = "LINEAR_RATE";
		break;
	case chrNumModel::CONST_RATE_BASE_NUMBER:
		res = "BASE_NUMBER";
		break;
	case chrNumModel::BASE_NUMBER_NO_DUPL:
		res = "BASE_NUMBER_NO_DUPL";
		break;
	case chrNumModel::CONST_RATE_EPS:
		res = "CONST_RATE_EPS";
		break;
	case chrNumModel::CONST_RATE_EPS_NO_DUPL:
		res = "CONST_RATE_EPS_NO_DUPL";
		break;
	case chrNumModel::GENERAL_CHR_MODEL:
		res = "GENERAL_CHR_MODEL";
		break;
	default:
		errorMsg::reportError("unknown type in chrNumberOptions::getModelTypeStr");
	}
	return res;
}

chrNumModel::rootFreqType chrNumberOptions::getRootFreqTypeFromStr(const string& str)
{
	if (str == "UNIFORM")
		return chrNumModel::UNIFORM;
	else if (str == "ROOT_LL")
		return chrNumModel::ROOT_LL;
	else if (str == "STATIONARY")
		return chrNumModel::STATIONARY;
	else if (str == "FIXED")
		return chrNumModel::FIXED;
	else
        errorMsg::reportError("unknown type in chrNumberOptions::getRootFreqTypeFromStr:" + str);
	return chrNumModel::ROOT_LL;
}

string chrNumberOptions::getRootFreqTypeStr(chrNumModel::rootFreqType type)
{
	string res = "";
	switch (type)
	{
	case chrNumModel::UNIFORM:
		res = "UNIFORM";
		break;
	case chrNumModel::ROOT_LL:
		res = "ROOT_LL";
		break;
	case chrNumModel::FIXED:
		res = "FIXED";
		break;
	case chrNumModel::STATIONARY:
		res = "STATIONARY";
		break;
	default:
		errorMsg::reportError("unknown type in chrNumberOptions::getRootFreqTypeStr");
	}
	return res;
}

chrNumModel::branchModelType chrNumberOptions::getBranchModelTypeFromStr(const string& str)
{
	if (str == "SPECIATIONAL")
		return chrNumModel::SPECIATIONAL;
	else if (str == "GRADUAL")
		return chrNumModel::GRADUAL;
	else if (str == "COMBINED")
		return chrNumModel::COMBINED;
	else if (str == "COMBINED_INTERNALS")
		return chrNumModel::COMBINED_INTERNALS;
	else
        errorMsg::reportError("unknown type in chrNumberOptions::getBranchModelTypeFromStr");
	return chrNumModel::GRADUAL;
}

string chrNumberOptions::getBranchModelTypeStr(chrNumModel::branchModelType type)
{
	string res = "";
	switch (type)
	{
	case chrNumModel::SPECIATIONAL:
		res = "SPECIATIONAL";
		break;
	case chrNumModel::GRADUAL:
		res = "GRADUAL";
		break;
	case chrNumModel::COMBINED:
		res = "COMBINED";
		break;
	case chrNumModel::COMBINED_INTERNALS:
		res = "COMBINED_INTERNALS";
		break;
	default:
		errorMsg::reportError("unknown type in chrNumberOptions::getBranchModelTypeStr");
	}
	return res;
}


vector<chrNumModel::modelType> chrNumberOptions::getModelsFromStr(const string& inStr)
{
	vector<chrNumModel::modelType> res;
	vector<string> outStr;
	splitString(inStr, outStr, ",");
	for (int i = 0; i < outStr.size(); ++i)
	{
		chrNumModel::modelType type = chrNumberOptions::getModelTypeFromStr(outStr[i]);
		res.push_back(type);
	}
	return res;
}

string chrNumberOptions::getStrFromModels(const vector<chrNumModel::modelType>& inVec)
{
	string res("");
	for (int i = 0; i < inVec.size(); ++i)
	{
		if (i > 0)
			res += ",";
		res += getModelTypeStr(inVec[i]);
	}
	return res;
}

void chrNumberOptions::getInitParams(vector<paramObj>& initParams)
{
	initParams.clear();
	initParams.push_back(paramObj(LOSS_CONST, _lossConstR));
	initParams.push_back(paramObj(GAIN_CONST, _gainConstR));
	initParams.push_back(paramObj(DUPL, _duplConstR));
	initParams.push_back(paramObj(LOSS_LINEAR, _lossLinearR));
	initParams.push_back(paramObj(GAIN_LINEAR, _gainLinearR));
	//initParams.push_back(paramObj(DUPL_LINEAR, _duplLinearR));
	initParams.push_back(paramObj(HALF_DUPL, _demiPloidyR));
	initParams.push_back(paramObj(SCALE_BRANCH, _scaleBranch));
	initParams.push_back(paramObj(BASE_NUMBER, _baseNumber));
	initParams.push_back(paramObj(BASE_NUMBER_R, _baseNumberR));
}
