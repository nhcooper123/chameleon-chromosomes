#ifndef _FINITE_INTEGER_ALPH
#define _FINITE_INTEGER_ALPH

#include "matrixUtils.h"
#include "alphabet.h"
#include "definitions.h"
#include "errorMsg.h"
#include <map>
/* 
this class represent an alphabet of positive integers. Zero is not allowed.
Note that the id of the characters are ("input char" -1).
That is: if the sequence file has a character 1 its id will be zero
*/

class finiteIntAlphabet : public alphabet {
public:
	explicit finiteIntAlphabet(int minInt, int maxInt); 
	explicit finiteIntAlphabet(const finiteIntAlphabet& other); 
	virtual finiteIntAlphabet& operator=(const finiteIntAlphabet &other);
	virtual ~finiteIntAlphabet();

	virtual alphabet* clone() const { return new finiteIntAlphabet(*this); }
	int unknown() const  {return -2;}
	int gap() const  {errorMsg::reportError("The method indel::gap() is used"); return -1;} // What is it for ? I don't need this !!!
	int min() const {return _min;}
	int max() const {return _max;}
	int size() const {return _size;}
	int getCompositeIdsNum() {return _relationsProbs.size() - size();};
	int stringSize() const {return 1000;} //The string size is undefined should check this 
	int relations(const int charInSeq, const int charToCheck) const;
	MDOUBLE relationsProbs(const int charInSeq, const int charToCheck) const;
	int fromChar(const string& str, const int pos) const;
	//int fromChar(const char s) const;
	string fromInt(const int in_id) const;
	vector<int> fromString(const string& str) const;
	bool isSpecific(const int id) const {return (id >= 0 && id < size());}
	bool isComposite(int id) const; 
	int getCompositeId(const string& compName) const; //return the composite id or -1 if the composite does not exist
	void getCompositeParts(int comId, Vint& compIds, Vdouble& compProbs) const; //return the ids and probabilities that constitute the composite id
	string setCompositeId(const string& taxaName, const string& compName); //creates a new composite ID. returns the compName as inserted in _id2Str

	//convert from the corrected id to the real count (e.g., id = 0 represents count of 1)
	int id2Count(int id) const;
	//convert from the corrected id to the real count (e.g., id = 0 represents count of 1)
	int count2Id(int count) const;


private:
	void copy(const finiteIntAlphabet& other);
	string getCompositeStr(vector<pair<int, MDOUBLE> > comId) const;
	string getCompositeStr(const string& name) const; //get the composite name as represented in _id2Str
private:
	int _max;
	int _min;
	int _size; //number of possible states
//	VVint _relations; //each entry represent the possible value of the id. 
					//each entry is a vector of size _size with 1 indicate possible identity. 
					//so that if id==11 represent 0 or 3 then _relations[11][0]=1 and _relations[11][3]=1 and all other places are zero.
					//for IDs between 0 and _size this is the identity matrix: _relations[x][x]=1, _relations[x][!x]=0
	VVdouble _relationsProbs;
	vector<string> _id2str; //maps the id to the composite count string
	map<string, int> _str2id;//maps the composite count string to id
};
#endif
