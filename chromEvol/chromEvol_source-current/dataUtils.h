#ifndef __DATA_UTILS___H
#define __DATA_UTILS___H
#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include <map>
using namespace std;

void convertToPhylipFormat(const string& inMsaFileName,  const string& outMsaFileName, const string& outConvertTableFileName);
void translatePhylipTree(const string& inTreeFileName,  const string& inConvertTableFileName, const string& outTreeFileName);
void rootTree(tree& inTree, const string& rootAt);
void removeTaxa(tree& inTree,  const string& missingListFile);
void convertNexusTreeToNewickFormat(const string& nexusTreeFileName,  const string& outTreeFileName);
void seperateNexusTrees(const string& nexusTreesFileName, const string& outDir);
void getCountsFromTreeNames(tree& inTree, const string countsFileName, const string& missingListFileName);
void buildScFromTree(tree &tr,sequenceContainer &origSc, map<string, int> &name2idSc, sequenceContainer &newSc);
sequenceContainer concatenateMsa(vector<string>& msaNames, alphabet* pAlph);

#endif
