#include "chrCountFormat.h"
#include "someUtil.h"
#include "errorMsg.h"
#include "ConversionUtils.h"
#include "sequence.h"
#include <algorithm>
using namespace std;

sequenceContainer chrCountFormat::read(istream &infile, finiteIntAlphabet* alph) {
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}


sequenceContainer chrCountFormat::readUnAligned(istream &infile, finiteIntAlphabet* pAlph) {
	sequenceContainer mySeqData;

	vector<string> seqFileData;
	putFileIntoVectorStringArray(infile,seqFileData);
	if (seqFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in fasta format");
	}

	vector<string>::iterator it1;
	int localid=0;
	for (it1 = seqFileData.begin(); it1!= seqFileData.end(); ) {
		*it1 = takeCharOutOfString(" \t", *it1);
		if (it1->empty()) {++it1;continue; }// empty line continue

		string remark;
		string name;

		if ((*it1)[0] == '>') {
			string::const_iterator itstrtmp = (*it1).begin();
			itstrtmp++;
			while (itstrtmp != (*it1).end()) {
				name+= *itstrtmp;
				itstrtmp++;
			}
			name = trim(name);

			//for (string::iterator i = name.begin(); i!=(name.end()-2);++i) {
			//	*i=*(i+1); // removing the ">". should be done more elegant...
			//}
			++it1;
		} else {
			LOG(0,<<"problem in line: "<<*it1<<endl);
			errorMsg::reportError("Error reading fasta file, error finding sequence name starting with >",1);
		}
		while (it1->empty()) it1++; // empty line continue
		
		string str;
		//while (it1!= seqFileData.end()) {
		//	if ((*it1)[0] == '>') break;
			str+=*it1;
			++it1;
		//}
		// remove spaces form str;
		str = takeCharOutOfString(" \t", str);
		vector<string> strVec;
		splitString(str, strVec, "_"); 
		Vint idsVec;
		if (strVec.size() > 1)
		{ //composite id ("17_34") or with probabilities ("17=0.6_34=0.4")
			str = pAlph->setCompositeId(name, str);
			//bool bProbVec = false;
			//if (str.find(":") != string::npos)
			//	bProbVec = true;
   //         for (int i = 0; i < strVec.size(); ++i)
			//{
			//	string s = strVec[i];
			//	if (bProbVec)
			//	{
			//		string count, probStr;
			//		splitString2(strVec[i], ":", count, probStr);
			//		int correctedId = atoi(count.c_str()) -1;
			//		MDOUBLE prob = atof(probStr.c_str());



			//	}
			//	else 
			//	{
   //                 int correctedId = atoi(strVec[i].c_str()) -1;
   //                 idsVec.push_back(correctedId);
			//		if (pAlph->getCompositeId(idsVec) == -1)
			//		pAlph->setCompositeId(idsVec);
			//	}
			//}
		}
		mySeqData.add(sequence(str,name,remark,localid,pAlph));
		localid++;
	}

	return mySeqData;
}


void chrCountFormat::write(ostream &out, const sequenceContainer& sd) {
	for (sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();it5!=sd.constTaxaEnd();++it5) {
		out<<">"<<(it5)->name()<<endl;
		out<<it5->toString()<<endl;
	}
}

//int chrCountFormat::getMaxIntInFile(const string &inFilName)
pair<int,int> chrCountFormat::getMinMaxCountInFile(const string &inFilName)
{
	ifstream inFile(inFilName.c_str());
	if (inFile.fail())
		errorMsg::reportError("cannot open file: " + inFilName);
	vector<string> seqFileData;
	putFileIntoVectorStringArray(inFile,seqFileData);
	inFile.close();
	if (seqFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in chrCountFormat format");
	}
	Vint allCounts;
	vector<string>::iterator it1;
	int localid=0;
	for (it1 = seqFileData.begin(); it1!= seqFileData.end(); ) {
		*it1 = takeCharOutOfString(" \t", *it1);
		if (it1->empty())
		{// empty line continue
			++it1;
			continue; 
		}
		if ((*it1)[0] == '>') 
		{
			++it1;
		} else {
			LOG(0,<<"problem in line: "<<*it1<<endl);
			errorMsg::reportError("Error reading fasta file, error finding sequence name starting with >",1);
		}
		while (it1->empty()) 
			it1++; // empty line continue
		
		string str;
		//while (it1!= seqFileData.end()) {
		//	if ((*it1)[0] == '>') break;
			str+=*it1;
			++it1;
		//}
		// remove spaces form str;
		str = takeCharOutOfString(" \t", str);
		vector<string> strVec;
		splitString(str, strVec, "_"); 
		if (strVec.size() > 1)
		{ //composite id ("17_34")
            for (int i = 0; i < strVec.size(); ++i)
			{
				int count = atoi(strVec[i].c_str());
				if (count > 0)
                    allCounts.push_back(count);
			}
		}
		else
		{
			int count = atoi(str.c_str());
			if (count > 0)
                allCounts.push_back(count);
		}
	}
	sort(allCounts.begin(), allCounts.end());
	pair<int,int> res(allCounts[0],allCounts[allCounts.size()-1]);
	return res;
}
