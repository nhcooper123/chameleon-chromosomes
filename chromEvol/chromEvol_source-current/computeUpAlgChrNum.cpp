#include "definitions.h"
#include "computeUpAlgChrNum.h"
#include "treeIt.h"
#include "seqContainerTreeMap.h"
#include "logFile.h"
#include "finiteIntAlphabet.h"
#include <iostream>
#include <cassert>
using namespace std;

void computeUpAlgChrNum::fillComputeUp(const tree& et,
								 const sequenceContainer& sc,
								 const int pos,
								 const computePijHom& pi,
								 suffStatGlobalHomPos& ssc) {

	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),pi.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter;
		//if (mynode->isLeaf()) {
		//	for(letter=0; letter<pi.alphabetSize();letter++) {
		//		const int seqID = sctm.seqIdOfNodeI(mynode->id());
		//		doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
		//		ssc.set(mynode->id(),letter,val);
		//	}
		//}
        if (mynode->isLeaf()) {
			MDOUBLE sum = 0.0;
			for(letter=0; letter<pi.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				//doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				doubleRep val = static_cast<const finiteIntAlphabet*>(sc.getAlphabet())->relationsProbs(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
				sum += convert(val);
			}
			/*		for(letter=0; letter<pi.alphabetSize();letter++) {
				doubleRep val = ssc.get(mynode->id(),letter);
				ssc.set(mynode->id(),letter,val/sum);
			}*/
		}

		else {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				doubleRep total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons();++i){				
					doubleRep prob=0.0;
					for(int letInSon=0; letInSon<pi.alphabetSize();letInSon++) {
						prob += ssc.get(mynode->getSon(i)->id(), letInSon)*
							pi.getPij(mynode->getSon(i)->id(),letter,letInSon);
					}
				total_prob*=prob;
				}
				ssc.set(mynode->id(),letter,total_prob);
			}
		}
	}
}
