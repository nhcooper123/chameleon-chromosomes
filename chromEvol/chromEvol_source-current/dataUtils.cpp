#include "dataUtils.h"
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
#include "fastaFormat.h"
#include <fstream>
#include <string>
#include <map>
#include <sstream>
using namespace std;


void convertToPhylipFormat(const string& inMsaFileName,  const string& outMsaFileName, const string& outConvertTableFileName)
{
	nucleotide alph;
	ifstream in(inMsaFileName.c_str());
	sequenceContainer inSc = recognizeFormat::read(in, &alph);
	in.close();

	ofstream outConvertF(outConvertTableFileName.c_str());
	for (int seq = 0; seq < inSc.numberOfSeqs(); ++seq)
	{
		string oldName = inSc[seq].name();
		string newName = "taxa_" + int2string(seq);
		inSc[seq].setName(newName);
		outConvertF<<oldName<<" "<<newName<<endl;
	}
	outConvertF.close();
	ofstream outMsaF(outMsaFileName.c_str());
	phylipFormat::write(outMsaF, inSc);
	outMsaF.close();
}


//gets a nexus tree and convert it to newick.
void convertNexusTreeToNewickFormat(const string& nexusTreeFileName,  const string& outTreeFileName)
{
	ifstream in(nexusTreeFileName.c_str());
	vector<string> lines;
	putFileIntoVectorStringArray(in, lines);
	in.close();
	map<string, string> tmpToRealName;
	
	vector<string>::const_iterator it1 = lines.begin();
	while ( ( (*it1).find("Translate")  == -1) && ((*it1).find("TRANSLATE")  == -1) &&(it1 != lines.end()))
	{ //check for the word Translate
		++it1;
	}
	if (it1 == lines.end())
		errorMsg::reportError("cannot find word Translate");
	bool bEndTable = false;
	for (++it1; bEndTable != true; ++it1)
	{
		if (it1->empty())
			continue;
		if ((*it1).find(";")  != -1)
			bEndTable = true;
		string realName, translate;
		//splitString2(*it1, " \t", translate, realName);
		vector<string> strVec;
		splitString(*it1, strVec, "\t "); 
		if (strVec.size() != 2)
		{
			if (bEndTable) //in case there was a line with only ";"
				break;
			errorMsg::reportError("unable to read line");
		}
		translate = strVec[0];
		realName = strVec[1];
		translate = takeCharOutOfString(" \t", translate);
		realName = takeCharOutOfString(",; ", realName);

		
		//int intTans = atoi(translate.c_str());
		//translate = int2string(intTans);
		tmpToRealName[translate] = realName;
	}
	
	//find first tree
	while ( ( (*it1).find("TREE")  == -1) && (it1 != lines.end()))
	{ 
		++it1;
	}
	//parse trees
	int n = 0;
	for (it1; ((*it1).find("END")  == -1); ++it1, ++n)
	{
		if (it1->empty())
			continue;
		vector<string> strVec;
		splitString(*it1, strVec, "="); 
		if (strVec.size() != 2)
			errorMsg::reportError("unable to read TREE line");
		string prefixStr = strVec[0];
		string treeStr = strVec[1];
		treeStr = takeCharOutOfString(" \t", treeStr);

		stringstream treeStream(treeStr);
		tree inTree(treeStream);
		//convert to real Names
		vector <tree::nodeP> leaves;
		inTree.getAllLeaves(leaves, inTree.getRoot());
		map<string, string>::const_iterator itr;
		for (int i=0; i< leaves.size(); ++i) {
			string tmpName = leaves[i]->name();
			if ((itr = tmpToRealName.find(tmpName)) ==  tmpToRealName.end())
			{
				errorMsg::reportError("cannot find tmp name");
			}
			string realName = itr->second;
			leaves[i]->setName(realName);
		}
		
		//in case of multiple trees:
		//string outTreeName = outDir + "//" + int2string(n) + ".tree";
		//inTree.output(outTreeName);
		inTree.output(outTreeFileName);
		break;
	}
}

//given a nexus tree file with many trees. 
//seperate the trees and output each one in a seperate file (newick format) inside the output directory
//have to make sure that the last line of the translation table is ";"
void seperateNexusTrees(const string& nexusTreesFileName, const string& outDir)
{
	ifstream in(nexusTreesFileName.c_str());
	vector<string> lines;
	putFileIntoVectorStringArray(in, lines);
	in.close();
	vector<string>::const_iterator it1 = lines.begin();
	
	while ( ( (*it1).find("BEGIN TREES")  == -1) && (it1 != lines.end()))
	{ //check for the word BEGIN TREES
		++it1;
	}

	while ( ( (*it1).find("TRANSLATE")  == -1) && (it1 != lines.end()))
	{ //check for the word Translate
		++it1;
	}
	if (it1 == lines.end())
		errorMsg::reportError("cannot find word TRANSLATE (upper case) or BEGIN TREES");
	//translate table: fill the map tmpToRealName
	//translation continue until there is a line with only ";"
	map<string, string> tmpToRealName;
	for (++it1; ((*it1).find(";")  == -1); ++it1)
	{
		if (it1->empty())
			continue;
		string realName, translate;
		//splitString2(*it1, " \t", translate, realName);
		vector<string> strVec;
		splitString(*it1, strVec, "\t "); 
		if (strVec.size() != 2)
			errorMsg::reportError("unable to read line");
		string s1 = strVec[0];
		string s2 = strVec[1];
		translate = strVec[0];
		realName = strVec[1];
		translate = takeCharOutOfString(" \t", translate);
		realName = takeCharOutOfString(",; ", realName);
		tmpToRealName[translate] = realName;
	}

	//find first tree
	while ( ( (*it1).find("TREE")  == -1) && (it1 != lines.end()))
	{ 
		++it1;
	}
	//parse trees
	int n = 0;
	for (it1; ((*it1).find("END")  == -1); ++it1, ++n)
	{
		if (it1->empty())
			continue;
		string prefixStr, treeStr;
		vector<string> strVec;
		splitString(*it1, strVec, "="); 
		if (strVec.size() != 2)
			errorMsg::reportError("unable to read TREE line");
		string s1 = strVec[0];
		string s2 = strVec[1];
		prefixStr = strVec[0];
		treeStr = strVec[1];
		treeStr = takeCharOutOfString(" \t", treeStr);

		stringstream treeStream(treeStr);
		tree inTree(treeStream);
		//convert to real Names
		vector <tree::nodeP> leaves;
		inTree.getAllLeaves(leaves, inTree.getRoot());
		map<string, string>::const_iterator itr;
		for (int i=0; i< leaves.size(); ++i){
			string tmpName = leaves[i]->name();
		if ((itr = tmpToRealName.find(tmpName)) ==  tmpToRealName.end())
		{
            errorMsg::reportError("cannot find tmp name");
		}
		string realName = itr->second;
			leaves[i]->setName(realName);
		}
		
		string outTreeName = outDir + "//" + int2string(n) + ".tree";
		inTree.output(outTreeName);
	}
}

void translatePhylipTree(const string& inTreeFileName,  const string& inConvertTableFileName, const string& outTreeFileName)
{
	ifstream in(inConvertTableFileName.c_str());
	vector<string> lines;
	putFileIntoVectorStringArray(in, lines);
	map<string, string> tmpToRealName;
	for (int i=0; i < lines.size();i++)
	{
		string realName, translate;
		splitString2(lines[i], " ", realName, translate);
		translate = takeCharOutOfString(" \t", translate);
		realName = takeCharOutOfString(" \t", realName);
		tmpToRealName[translate] = realName;
	}

	in.close();

	tree inTree(inTreeFileName);
	vector <tree::nodeP> leaves;
	inTree.getAllLeaves(leaves, inTree.getRoot());
	map<string, string>::const_iterator itr;
	for (int i=0; i< leaves.size(); ++i){
		string tmpName = leaves[i]->name();
		if ((itr = tmpToRealName.find(tmpName)) ==  tmpToRealName.end())
		{
			errorMsg::reportError("cannot find tmp name");
		}
		string realName = itr->second;
		leaves[i]->setName(realName);
	}
	inTree.output(outTreeFileName);

}


void rootTree(tree& inTree, const string& rootAt)
{
	if (!(rootAt == "")){
		tree::nodeP myroot = inTree.findNodeByName(chrNumberOptions::_rootAt); //returns NULL if not found
		if (myroot){
			inTree.rootAt(myroot);
			LOG(chrNumberOptions::_logValue, <<"tree rooted at "<<myroot->name()<<" id, "<<myroot->id()<<endl);
			LOG(chrNumberOptions::_logValue, <<"sons of root are "<<inTree.getRoot()->getSon(0)->name()<<" , "<<inTree.getRoot()->getSon(1)->name()<<" , "<<inTree.getRoot()->getSon(2)->name()<<endl);
			return;
		}
	}
	
}


void removeTaxa(tree& inTree,  const string& missingListFile){
	
	ifstream m_in(missingListFile.c_str());
	vector<string> names;
	putFileIntoVectorStringArray(m_in,names);

	for (int i=0; i<names.size();i++){
		tree::nodeP thisNode = inTree.findNodeByName(names[i]);
		if (thisNode) {
			if ((thisNode->father()->isRoot()) && (inTree.getRoot()->getNumberOfSons() == 2))
			{
				//in case the tree was rooted and the removed leaf was one of the root' sons:
				//we have to remove the root and reroot the tree at the second root son
				tree::nodeP pRoot = inTree.getRoot();
				tree::nodeP otherSonOfRoot;
				if (thisNode == pRoot->getSon(0))
					otherSonOfRoot = pRoot->getSon(1);
				else
					otherSonOfRoot = pRoot->getSon(0);

				tree newTree;
				newTree.createRootNode();
				newTree.getRoot()->setName(otherSonOfRoot->name());
				newTree.recursiveBuildTree(newTree.getRoot(),otherSonOfRoot->getSon(0));
				newTree.recursiveBuildTree(newTree.getRoot(),otherSonOfRoot->getSon(1));
				inTree = newTree;
			}
			else
                inTree.removeLeaf(thisNode);
		}
		else
		{
			cerr<<"Error in removeTaxa. Cannot find: "<<names[i]<<" in tree"<<endl;
			//errorMsg::reportError(err);
		}
	}
}

//the counts are given in the leaf name as name_count.
//create the countsFile in a fasta format and put taxa with no counts in missingListFile
void getCountsFromTreeNames(tree& inTree, const string countsFileName, const string& missingListFileName)
{
	ofstream countsF(countsFileName.c_str());
	ofstream missingF(missingListFileName.c_str());

	vector <tree::nodeP> leaves;
	inTree.getAllLeaves(leaves, inTree.getRoot());
	for (int i=0; i< leaves.size(); ++i){
		bool found = false;
		string name = leaves[i]->name();
		int pos = name.find_last_of("_");
		if (pos != -1)
		{
			string suffix =  name.substr(pos+1);
			int count;
			string::const_iterator itBegin = suffix.begin();
			string::const_iterator itEnd = suffix.end();
			if (fromStringIterToInt(itBegin, itEnd, count) ==true)
			{
				found = true;
				countsF<<">"<<name<<endl;
				countsF<<count<<endl;
			}
		}
		if (found == false)
		{
			missingF<<name<<endl;
		}
	}
	countsF.close();
	missingF.close();
	
}


void buildScFromTree(tree &tr,sequenceContainer &origSc, map<string, int> &name2idSc, sequenceContainer &newSc)
{
	vector<tree::nodeP> nodesVecTree1;
	tr.getAllLeaves(nodesVecTree1, tr.getRoot());
	cout<<"building new sequence container, root is "<<tr.getRoot()->name()<<endl;
	for (int i = 0; i<nodesVecTree1.size(); ++i){
		int idInOrigSc = -1;
		string name = nodesVecTree1[i]->name();
		map<string,int>::iterator iter = name2idSc.find(name);
		if( iter == name2idSc.end() ) {
			string errorMs = "error in buildScFromTree, cannot find ";
			errorMs+=nodesVecTree1[i]->name();
			errorMs+=" in sequence container";
			errorMsg::reportError(errorMs);
		}
		else {
			idInOrigSc = iter->second;
		}
		int content = origSc[idInOrigSc][0];//first position
		sequence ss (int2string(content),name,"",i,origSc.getAlphabet());
		newSc.add(ss);
	}
}

sequenceContainer concatenateMsa(vector<string>& msaNames, alphabet* pAlph)
{
	if (msaNames.size() < 2)
		errorMsg::reportError("less than 2 MSAs were given to concatenate"); 
	ifstream inFile(msaNames[0].c_str());
	sequenceContainer sc = recognizeFormat::read(inFile, pAlph);
	inFile.close();
	for (int i = 1; i < msaNames.size(); ++i)
	{
		ifstream inF(msaNames[i].c_str());
		sequenceContainer tmpSc = recognizeFormat::read(inF, pAlph);
		inF.close();
		sc.concatenate(tmpSc);
	}
	return sc;

}

//void sampleSeqs(string& inTreeFileName, string& inMSAFileName, MDOUBLE inSamplePercent, string& outTreeFileName, string& outMSAFileName)
//{
//	tree inTree(inTreeFileName);
//	threeStateAlphabet alph;
//	ifstream inMSAFile(inMSAFileName.c_str());
//	sequenceContainer inSc = recognizeFormat::read(inMSAFile, &alph);
//	inMSAFile.close();
//
//	int seqNumToSample = static_cast<int>(inSamplePercent * inSc.numberOfSeqs());
//
//	///////sample sequences from MSA
//	//create a new sc with the sampled sequences and a vector of names of leaves to be removed from the tree
//	if (seqNumToSample > inSc.numberOfSeqs())
//		errorMsg::reportError("the number of requested seqeuences is larger than the number of sequences in the MSA");
//	cerr<<"eee";
//	sequenceContainer newSc;
//	Vint vec2Add(inSc.numberOfSeqs(),0);
//	Vstring seqToRemove; 
//	int n = 0;
//	while (n < seqNumToSample)
//	{
//		int seqPlaceToAdd = talRandom::giveIntRandomNumberBetweenZeroAndEntry(inSc.numberOfSeqs());
//		if (vec2Add[seqPlaceToAdd] == 0){
//			vec2Add[seqPlaceToAdd] = 1;
//			n++;
//		}
//        
//	}
//	for (int i = 0; i<vec2Add.size();i++){
//		if (vec2Add[i] == 1)
//			newSc.add(inSc[i]);
//		else
//			seqToRemove.push_back(inSc[i].name());
//	}
//
//	//filter tree
//	for (int i=0; i<seqToRemove.size();i++){
//		tree::nodeP thisNode = inTree.findNodeByName(seqToRemove[i]);
//		if (thisNode) {
//			inTree.removeLeaf(thisNode);
//		}
//		else
//		{
//			errorMsg::reportError("cannot find: " + seqToRemove[i]);
//		}
//	}
//	inTree.output(outTreeFileName);
//	ofstream outMSA(outMSAFileName.c_str());
//	fastaFormat::write(outMSA, newSc);
//	outMSA.close();
//}

