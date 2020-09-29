#ifndef _EVAL_PP_DIST_MNG_H
#define _EVAL_PP_DIST_MNG_H

#include "definitions.h"
#include "tree.h"


class evalPPDistMng {

public:
	explicit evalPPDistMng();
	virtual ~evalPPDistMng();
    void run();
	

private:
	void init(); 
	void readExpectations();
	void calcExpectedAndObservedCounts();
	void calcExpectedAndObservedCountsLeaves(); //chech external edges versus internal
	void calculateEdgeHeights();
	void calcChiStat();
	void getTree();
	MDOUBLE adjustExpectationsAndTreeLengthExcludingLongExternals(const tree &tr, MDOUBLE cutoff); //calculate sum of branch lengths beside those external branches that are longer than the given cutoff
	
private:
	tree _tree;
	VVdouble _expChanges; // for each edge, the expected number of changes for [0-3]: gain, loss, dupl, halfDupl
	Vdouble _expChangesSet1;
	Vdouble _expChangesSet2;
	Vdouble _obsChangesSet1;
	Vdouble _obsChangesSet2;
	Vdouble _cutoffs; // a vector of size 4: lowerHeightSet1,upperHeightSet1,lowerHeightSet2,upperHeightSet2
	vector<pair<MDOUBLE, MDOUBLE> > _edgeHeights; //holds the lower and upper heights for each edge - specified by the lower node id
};
#endif
