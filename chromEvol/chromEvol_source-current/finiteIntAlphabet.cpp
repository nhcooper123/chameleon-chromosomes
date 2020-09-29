#include "finiteIntAlphabet.h"
#include "someUtil.h"
#include "matrixUtils.h"

#include <algorithm>
using namespace std;

finiteIntAlphabet::finiteIntAlphabet(int minInt, int maxInt) 
: _max(maxInt), _min(minInt)
{
	if ((_max - _min) < 1)
		errorMsg::reportError("error in finiteIntAlphabet. The size of the data cannot be lower than 1");
	_size = maxInt-minInt+1;
	//unitMatrix(_relations, _size);
	unitMatrix(_relationsProbs, _size);
}

finiteIntAlphabet::finiteIntAlphabet(const finiteIntAlphabet& other) {
	copy(other);
}

void finiteIntAlphabet::copy(const finiteIntAlphabet& other) {
	_min= other._min;
	_max= other._max;
	_size = other._size;
	//_relations = other._relations;
	_relationsProbs = other._relationsProbs;
	_id2str = other._id2str;
	_str2id = other._str2id;
}

finiteIntAlphabet& finiteIntAlphabet::operator=(const finiteIntAlphabet &other) {
	copy(other);
	return *this;
}

finiteIntAlphabet::~finiteIntAlphabet() {
}


//should correct this function to allow for uncertainties (1,2--> the count is 1 or 2)
int finiteIntAlphabet::fromChar(const string& str, const int pos) const{
	if (pos != 0)
		errorMsg::reportError("finiteIntAlphabet currently support only one position");

	vector<string> strVec;
	splitString(str, strVec, "_"); 
	Vint countVec;
    if (strVec.size() < 1)
		errorMsg::reportError("unable to read line:" + str);
	else if (strVec.size() == 1)
	{
		int correctedId = count2Id(atoi(strVec[0].c_str()));
		if (isSpecific(correctedId))
			return correctedId;
		if ((str[0] == 'X') || (str[0] == 'x'))
			return unknown();
		vector<string> err;
		err.push_back(" The chromosome number sequences contained the character: ");
		err.push_back(str);
		err.push_back(" while the maximum chromosome number allowed is: ");
		err.push_back(int2string(_size));
		errorMsg::reportError(err);
		return -99; // never suppose to be here.	
	}
	//composite ID
	//for (int i = 0; i < strVec.size(); ++i)
	//{
	//	int correctedId = atoi(strVec[i].c_str()) - 1;
	//	countVec.push_back(correctedId);
	//}
	int compId = getCompositeId(str);
	if (compId == -1)
		errorMsg::reportError("cannot find composite ID: " + str);
	return compId;
}

vector<int> finiteIntAlphabet::fromString(const string &str) const {
	vector<int> vec;
	for (int i=0;i<str.size();i++)
		vec.push_back(fromChar(str, 0));
	return vec;
}

string finiteIntAlphabet::fromInt(const int in_id) const{
	string res;
	if (in_id == unknown())
		res = "X";
	else if (isSpecific(in_id)) 
	{
		res = int2string(_min +in_id);
	}
	else if (isComposite(in_id)) 
	{
		res = _id2str[in_id-_size];
	}
	else
	{
		vector<string> err;
		err.push_back("unable to print id: " + int2string(in_id));
		errorMsg::reportError(err);
	}
	return res;
}


int finiteIntAlphabet::relations(const int charInSeq, const int charToCheck) const {
	errorMsg::reportError("relations is implemented as relationsProbs in finiteIntAlphabet::relations");
	return -1;
	//if (charInSeq == unknown())
	//	return 1;
	//if ((charInSeq < 0) || (charInSeq >= _relations.size()))
 //       return 0;
 //   int xx = _relations[charInSeq][charToCheck];
	//return _relations[charInSeq][charToCheck];
}

MDOUBLE finiteIntAlphabet::relationsProbs(const int charInSeq, const int charToCheck) const {
	if (charInSeq == unknown())
		return 1;
	if ((charInSeq < 0) || (charInSeq >= _relationsProbs.size()))
        return 0;
    return _relationsProbs[charInSeq][charToCheck];
}


int finiteIntAlphabet::getCompositeId(const string& str) const
{
	string compName = getCompositeStr(str);
	map<string, int>::const_iterator itr;
	if ((itr = _str2id.find(compName)) ==  _str2id.end())
	{
		return -1;
	}
	else
		return itr->second;
}


string finiteIntAlphabet::setCompositeId(const string& taxaName, const string& compNameOrig) 
{
	vector<pair<int, MDOUBLE> > compId;
	vector<string> strVec;
    splitString(compNameOrig, strVec, "_"); 
	Vint idsVec;
	if (strVec.size() < 2)
        errorMsg::reportError("there is no composite id here");
	 //composite id ("17_34") or with probabilities ("17=0.6_34=0.4")
	Vdouble probs(_size, 0);
	MDOUBLE sum = 0.0;
	bool bProbVec = false;
	for (int i = 0; i < strVec.size(); ++i)
	{
		string s = strVec[i];
		string count, probStr;
		splitString2(strVec[i], "=", count, probStr);
		MDOUBLE prob;
		if (probStr == "NULL")
			prob = 1.0;
		else
		{
			prob = atof(probStr.c_str());
			sum += prob;
			bProbVec = true;
		}
		int correctedId = count2Id(atoi(count.c_str()));
		if ((prob <0) || (prob > 1))
			errorMsg::reportError("not a probability value in: " +compNameOrig); 

		compId.push_back(pair<int, MDOUBLE>(correctedId, prob));
		probs[correctedId] = prob;
	}

	//string compName = getCompositeStr(compIds);
	//Vint ids(_size, 0);
	//MDOUBLE sum = 0.0;
	//for (int i = 0; i < compIds.size(); ++i)
	//{
	//	ids[compIds[i]] = 1;
	//	if (i == 0)
	//		probs[compIds[i]] = 0.7;
	//	else
	//		probs[compIds[i]] = 0.3;
	//}
	if (!DEQUAL(sum, 1.0) && bProbVec)
		errorMsg::reportError("The count probabilities for taxa " + taxaName + "do not sum to 1.0");
	string compName = getCompositeStr(compNameOrig);
	map<string, int>::const_iterator itr;
	if ((itr = _str2id.find(compName)) ==  _str2id.end())
	{
		int nextId = _relationsProbs.size();
        _str2id[compName] = nextId;
	//	_relations.push_back(ids);
		_relationsProbs.push_back(probs);
		_id2str.push_back(compName);
	}
	return compName;
}

string finiteIntAlphabet::getCompositeStr(vector<pair<int, MDOUBLE> > compId) const
{
	if (compId.size() < 2)
		errorMsg::reportError("there is no composite id here");
	string compName("");
	for (int i = 0; i < compId.size(); ++i)
	{
		if (i > 0)
			compName += "_";
		if (!isSpecific(compId[i].first))
			errorMsg::reportError("unknown count:" + int2string(compId[i].first));
		int count = id2Count(compId[i].first);
        compName += int2string(count);
		compName += "=";
		compName += double2string(compId[i].second, 2, true);
	}
	return compName;
}

bool finiteIntAlphabet::isComposite(int id) const
{
	int compId = id-_size;
	if ((compId >= _id2str.size()) || (compId < 0))
		return false;
	return true;
}
 //return the ids that constitute the composite id
void finiteIntAlphabet::getCompositeParts(int compId, Vint& compIds, Vdouble& compProbs) const
{
	if (!isComposite(compId)) {
		compIds.push_back(compId);
		compProbs.push_back(1.0);
        return;
	}
	string compStr = fromInt(compId);
	vector<string> strVec;
	splitString(compStr, strVec, "_"); 
	if (strVec.size() < 2)
		errorMsg::reportError("composite name contains only one id. in finiteIntAlphabet::getCompositeParts"); 
	for (int i = 0; i < strVec.size(); ++i)
	{
		string count, probStr;
		splitString2(strVec[i], "=", count, probStr);
		if (probStr == "NULL")
			probStr = "1.0";
		int correctedId = count2Id(atoi(count.c_str()));


		//int correctedId = count2Id(atoi(strVec[i].c_str()));
		compIds.push_back(correctedId);
		MDOUBLE prob = atof(probStr.c_str());
		compProbs.push_back(prob);
	}
	return;
}

//convert from the corrected id to the real count (e.g., id = 0 represents count of 1)
int finiteIntAlphabet::id2Count(int id) const
{
	if (id == unknown())
		return unknown(); 
	if (isComposite(id))
		return _max*10;
	string intStr = fromInt(id);
	return atoi(intStr.c_str());
}

int finiteIntAlphabet::count2Id(int count) const
{
	return count - _min;
}

//get the composite name as represented in _id2Str
string finiteIntAlphabet::getCompositeStr(const string& name) const 
{
	vector<pair<int, MDOUBLE> > compId;
	vector<string> strVec;
    splitString(name, strVec, "_"); 
	Vint idsVec;
	if (strVec.size() < 2)
        errorMsg::reportError("there is no composite id here");
	 //composite id ("17_34") or with probabilities ("17=0.6_34=0.4")
	Vdouble probs(_size, 0);
	for (int i = 0; i < strVec.size(); ++i)
	{
		string s = strVec[i];
		string count, probStr;
		splitString2(strVec[i], "=", count, probStr);
		if (probStr == "NULL")
			probStr = "1.0";
		int correctedId = count2Id(atoi(count.c_str()));
		MDOUBLE prob = atof(probStr.c_str());
		compId.push_back(pair<int, MDOUBLE>(correctedId, prob));
		probs[correctedId] = prob;
	}
	string compName = getCompositeStr(compId);
	return compName;
}
