#ifndef ___COMPUTE_UP_ALG_CHRNUM
#define ___COMPUTE_UP_ALG_ALG_CHRNUM

#include "definitions.h"
#include "tree.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"


class computeUpAlgChrNum {
public: 
	void fillComputeUp(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc);

};
#endif


