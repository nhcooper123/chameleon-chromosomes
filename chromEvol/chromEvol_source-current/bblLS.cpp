#include "bblLS.h"
#include "numRec.h"


bblLS::bblLS()
{}


MDOUBLE bblLS::optimizeBranches(tree& et, const sequenceContainer &intronSc, const sequenceContainer &exonSc,
						exonIntronModel& model, int maxIter, MDOUBLE epsilon, MDOUBLE curL)
{
	if (curL == NULL)
        _treeLikelihood = model.compLogLikelihood(et, intronSc, exonSc);
	else
		_treeLikelihood  = curL;
	cerr<<"============================="<<endl;
	cerr<<"ll before bbl = "<<_treeLikelihood<<endl;
	vector<tree::nodeP> nodesV;
	et.getAllNodes(nodesV,et.getRoot());
	MDOUBLE prevIterL = VERYSMALL;
	for (int iter = 0; iter < maxIter; ++iter) 
	{
        if (_treeLikelihood < prevIterL + epsilon)
			return _treeLikelihood; //likelihood converged
		prevIterL = _treeLikelihood;
		LOG(3,<<"bbl iteration: "<<iter<<" begin"<<endl;);
		MDOUBLE paramFound;
		MDOUBLE oldBl;
		MDOUBLE newL;
		for (int i=0; i<nodesV.size(); i++)
		{
			if (nodesV[i]->isRoot()) 
				continue;
			oldBl = nodesV[i]->dis2father();
			newL = -brent(0.0, oldBl, MAX_BRANCH_LENGTH, evalBranch(nodesV[i],&et, intronSc, exonSc, model), epsilon, &paramFound); 
			if (newL >= _treeLikelihood) 
			{
				_treeLikelihood = newL;
				nodesV[i]->setDisToFather(paramFound);
			} 
			else //likelihood went down!
			{
				nodesV[i]->setDisToFather(oldBl); //return to previous BL
				LOG(3,<<"like went down when branch was optimized: "<<endl;);
				LOG(3,<<"BL Found... "<<paramFound<<"...LL="<<newL<<"..."<<endl;);
				LOG(3,<<"BL old... "<<oldBl<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			}			
			LOG(3,<<"BL old... "<<oldBl<<" BL done... "<<nodesV[i]->dis2father()<<"...LL="<<_treeLikelihood<<"..."<<endl;);
		}
	}
	return _treeLikelihood;
}

MDOUBLE bblLS::optimizeBranches(tree& et, const sequenceContainer &sc, codonPositionModel& model, int maxIter, MDOUBLE epsilon, MDOUBLE curL)
{
	if (curL == NULL)
        _treeLikelihood = model.compLogLikelihood(et,sc);
	else
		_treeLikelihood  = curL; //alon curL means current likelihood
	cerr<<"============================="<<endl;
	cerr<<"ll before bbl = "<<_treeLikelihood<<endl;
	vector<tree::nodeP> nodesV;
	et.getAllNodes(nodesV,et.getRoot());
	MDOUBLE prevIterL = VERYSMALL; //alon VERYSMALL is about -inf
	for (int iter = 0; iter < maxIter; ++iter) 
	{
        if (_treeLikelihood < prevIterL + epsilon)
			return _treeLikelihood; //likelihood converged
		prevIterL = _treeLikelihood;
		LOG(3,<<"bbl iteration: "<<iter<<" begin"<<endl;);
		MDOUBLE paramFound;
		MDOUBLE oldBl;
		MDOUBLE newL;
		for (int i=0; i<nodesV.size(); i++)
		{
			if (nodesV[i]->isRoot()) 
				continue;
			oldBl = nodesV[i]->dis2father();
			newL = -brent(0.0, oldBl, MAX_BRANCH_LENGTH, evalCodonPositionBranch(nodesV[i],&et, sc, model), epsilon, &paramFound); 
			if (newL >= _treeLikelihood) 
			{
				_treeLikelihood = newL;
				nodesV[i]->setDisToFather(paramFound);
			} 
			else //likelihood went down!
			{
				nodesV[i]->setDisToFather(oldBl); //return to previous BL
				LOG(3,<<"like went down when branch was optimized: "<<endl;);
				LOG(3,<<"BL Found... "<<paramFound<<"...LL="<<newL<<"..."<<endl;);
				LOG(3,<<"BL old... "<<oldBl<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			}			
			LOG(3,<<"BL old... "<<oldBl<<" BL done... "<<nodesV[i]->dis2father()<<"...LL="<<_treeLikelihood<<"..."<<endl;);
		}
	}
	return _treeLikelihood;
}






MDOUBLE bblLS::alon_optimizeBranches(tree& et, const sequenceContainer &sc,codonPositionModelLight& model, int maxIter, MDOUBLE epsilon, MDOUBLE curL)
{
	epsilon = epsilon;

	if (curL == NULL)
        _treeLikelihood = model.compLogLikelihood(et,sc); //set current likelihood
	else
		_treeLikelihood  = curL;

	MDOUBLE prevIterL = VERYSMALL; //set liklihood of previos iteration to -INF so we not stop before first stage
	vector<MDOUBLE> gradient;

	MDOUBLE max_branch_length = 10;
	MDOUBLE min_branch_length = 0;


		tree leftPoint = et;
		vector<tree::nodeP> nodesLeft;
		leftPoint.getAllNodes(nodesLeft,leftPoint.getRoot());
		tree rightPoint = et;
		vector<tree::nodeP> nodesRight;
	
	bool hit_the_edge = false;
	for(int iter = 0; iter < maxIter; ++iter)
	{
		if (_treeLikelihood < prevIterL + epsilon && !(hit_the_edge))
		{
			cout<<endl<<"done gradient round _treeLikelihood is:  "<<_treeLikelihood<<endl;
			et = leftPoint;
			return _treeLikelihood;
		}
			
		else
		{
			prevIterL = _treeLikelihood;
			hit_the_edge = true; // so the next iteration will have to change it back
		}

		nodesLeft.clear();
		leftPoint.getAllNodes(nodesLeft,leftPoint.getRoot());
		nodesRight.clear();
		rightPoint.getAllNodes(nodesRight,rightPoint.getRoot());
		
		gradient = alon_aid_find_gradient(leftPoint, sc, model, nodesLeft, epsilon/1000, _treeLikelihood, min_branch_length, max_branch_length); //todo alon chang to get the gradient in the right size!
		

		cout<<"gradient:   "<<endl;
		for(int i=0; i<gradient.size(); i++)
		{
			cout<<endl<<gradient[i]<<endl;			
		}


		cout<<"point before:   "<<endl;
		for(int i=0; i<nodesLeft.size(); i++)
		{
			cout<<endl<<nodesLeft[i]->dis2father();
			
		}

		for(int i=0; i<nodesRight.size(); i++)
		{
			nodesRight[i]->setDisToFather(nodesLeft[i]->dis2father()+gradient[i]);	//make sure it is not negative!
		}



		cout<<endl<<"point after:   "<<endl;
		for(int i=0; i<nodesRight.size(); i++)
		{
			cout<<endl<<nodesRight[i]->dis2father();
			
		}
		
		cout<<endl<<"xxxxxxxxxxxxxxxxxxxx"<<endl;
		cout<<endl<<_treeLikelihood<<endl;
		MDOUBLE rightPointLikelihood = 1000;

		//do line minimization
		MDOUBLE leftPointLikelihood = _treeLikelihood;
		string valid = "left";

		while(alon_aid_distance_of_two_points(leftPoint,rightPoint)>epsilon) //alon todo not value diff but distance diff + add var that remembers if we hit a wall last time
		{
		alon_aid_gradient_walk_step(leftPoint,rightPoint,sc,model,leftPointLikelihood,rightPointLikelihood,valid);
		_treeLikelihood = leftPointLikelihood;
		if (valid=="left") // we don't go with the gradient all the way.
			hit_the_edge = false;

		cout<<endl<<leftPointLikelihood;
		cout<<endl<<rightPointLikelihood;
		cout<<endl<<valid;
		}

		leftPoint = rightPoint; //in case right point is on the edge and left point is not.

	}
	et = leftPoint;
	cout<<endl<<"done1";
	return _treeLikelihood;

}

//the gradient is not only the direction but also size. Here, we want to find the gradient vector that will hit the edge of the parameter space. 
vector<MDOUBLE> bblLS::alon_aid_find_gradient(tree& et, const sequenceContainer &sc,codonPositionModelLight& model, vector<tree::nodeP> nodesV, MDOUBLE epsilon, MDOUBLE &curL, MDOUBLE min_value, MDOUBLE max_value)
{
	nodesV.clear();
	et.getAllNodes(nodesV,et.getRoot());
	curL = model.compLogLikelihood(et,sc);
	MDOUBLE tempL;
	MDOUBLE oldDisToFather;
	vector<MDOUBLE> gradient;
	gradient.push_back(0); //for the root
	for(int i=1; i<nodesV.size(); i++)  //we start from 1 because nodesV[0] is the root and has no braches above it
	{
		oldDisToFather = nodesV[i]->dis2father();
		nodesV[i]->setDisToFather(oldDisToFather+epsilon);
		tempL = model.compLogLikelihood(et,sc);
		if ((tempL>curL && (oldDisToFather+epsilon)>max_value) || (tempL<curL && (oldDisToFather-epsilon)<min_value)) //so we won't go off the edge
		{
			gradient.push_back(0);
		}
		else
		{
		gradient.push_back((tempL-curL)/epsilon); 
		}
		nodesV[i]->setDisToFather(oldDisToFather);
	}


	//resize gradient
	for(int i=0; i<gradient.size(); i++)
	{
		cout<<endl<<"nodesV[i] "<<nodesV[i]->dis2father();
	}

	vector<MDOUBLE> limit;
	MDOUBLE size = 1000/epsilon; 
	MDOUBLE temp_size = 1000/epsilon;
	//here we want to find the shortest addition of the gradient from the current point that will bring us to the edge of the parameter space. 
	//size is the multiplier we want to find (we will multiply the gradient by this number)
	for(int i=0; i<gradient.size(); i++)
	{
		//if (gradient[i]==0) do nothing
//		if (gradient[i]>0)
		//	temp_size = (max_value-nodesV[i]->dis2father())/gradient[i]; //the division by gradient[i] will rescale by 1.0
		//if (gradient[i]<0)
		//	temp_size = (min_value-nodesV[i]->dis2father())/gradient[i];
		MDOUBLE step_wanted;
		if (gradient[i]>0)
			step_wanted = max_value-nodesV[i]->dis2father(); //the division by gradient[i] will rescale by 1.0
		if (gradient[i]<0)
			step_wanted = min_value-nodesV[i]->dis2father();
		temp_size = step_wanted/gradient[i];


		if(temp_size<size && temp_size>0)
			size = temp_size;
		cout<<"temp_size:  "<<temp_size<<"size:  "<<size;

	}

	for(int i=0; i<gradient.size(); i++)
	{
		cout<<endl<<"gradient[i] "<<gradient[i];
		gradient[i] = gradient[i]*size;
		cout<<endl<<"gradient[i] "<<gradient[i];
	}

	return gradient;
}

MDOUBLE bblLS::alon_aid_distance_of_two_points(tree leftPoint, tree rightPoint){
	
	vector<tree::nodeP> nodesLeft;
	vector<tree::nodeP> nodesRight;
	leftPoint.getAllNodes(nodesLeft,leftPoint.getRoot());
	rightPoint.getAllNodes(nodesRight,rightPoint.getRoot());

	int dimension = nodesLeft.size();
	MDOUBLE dist=0;
		for(int i=0;i<dimension;i++){
			dist+=pow((nodesLeft[i]->dis2father()-nodesRight[i]->dis2father()),2);
		}
		dist=sqrt(dist);
		return dist;
}

void bblLS::alon_aid_gradient_walk_step(tree& leftPoint, tree& rightPoint, const sequenceContainer &sc,codonPositionModelLight& model, MDOUBLE &left_likelihood, MDOUBLE &right_likelihood, string &valid){

	tree middlePoint = leftPoint;
	vector<tree::nodeP> nodesLeft;
	vector<tree::nodeP> nodesRight;
	vector<tree::nodeP> nodesMiddle;
	leftPoint.getAllNodes(nodesLeft,leftPoint.getRoot());
	rightPoint.getAllNodes(nodesRight,rightPoint.getRoot());
	middlePoint.getAllNodes(nodesMiddle,middlePoint.getRoot());
	if(valid != "left")
	{
		cout<<endl<<"computing left";
		left_likelihood = model.compLogLikelihood(leftPoint,sc);
	}
	if(valid != "right")
	{
		cout<<endl<<"computing right";
		right_likelihood = model.compLogLikelihood(rightPoint,sc);
	}

	for(int i=0; i<nodesLeft.size(); i++)
	{
			nodesMiddle[i]->setDisToFather((nodesRight[i]->dis2father()+nodesLeft[i]->dis2father())/2);	
	}

	if(left_likelihood > right_likelihood)
	{
		rightPoint=middlePoint;
		valid = "left";
	}
	else
	{
		leftPoint=middlePoint;
		valid = "right";
	}
}




MDOUBLE bblLS::optimizeBranches(tree& et, const sequenceContainer &sc, codonPositionModelLight& model, int maxIter, MDOUBLE epsilon, MDOUBLE curL)
{
	if (curL == NULL)
        _treeLikelihood = model.compLogLikelihood(et,sc);
	else
		_treeLikelihood  = curL;
	cerr<<"============================="<<endl;
	cout<<"ll before bbl "<<time(0)<<endl; //alon
	cerr<<"ll before bbl = "<<_treeLikelihood<<endl;
	vector<tree::nodeP> nodesV;
	et.getAllNodes(nodesV,et.getRoot());
	MDOUBLE prevIterL = VERYSMALL;
	
	for (int iter = 0; iter < maxIter; ++iter) 
	{
		cout<<"bblLS::optimizeBranches iteration: "<<iter<<"of "<<maxIter<<", time: "<<time(0)<<endl; //alon
        if (_treeLikelihood < prevIterL + epsilon)
			return _treeLikelihood; //likelihood converged
		prevIterL = _treeLikelihood;
		LOG(3,<<"bbl iteration: "<<iter<<" begin"<<endl;);
		MDOUBLE paramFound;
		MDOUBLE oldBl;
		MDOUBLE newL;
		for (int i=0; i<nodesV.size(); i++) //for all nodes in the tree
		{
			//alon
			cout<<endl<<"this is an INNER branch length vector:   "<<endl;
			for(int i=0; i<nodesV.size(); i++)
			{
				cout<<endl<<nodesV[i]->dis2father();
			
			}
			//alon
			cout<<endl<<"bblLS::optimizeBranches inner-loop iteration: "<<i<<" of "<<nodesV.size()<<", time: "<<time(0)<<endl; //alon
			if (nodesV[i]->isRoot()) // if node is root we can't compute the distance from its father
				continue;
			oldBl = nodesV[i]->dis2father();
			newL = -brent(0.0, oldBl, MAX_BRANCH_LENGTH, evalCodonPositionBranchLight(nodesV[i],&et, sc, model), epsilon, &paramFound); //finds a new likelihood after optimizing the length of the brance that starts with node i
			if (newL >= _treeLikelihood) 
			{
				_treeLikelihood = newL;
				nodesV[i]->setDisToFather(paramFound);
			} 
			else //likelihood went down!
			{
				nodesV[i]->setDisToFather(oldBl); //return to previous BL
				LOG(3,<<"like went down when branch was optimized: "<<endl;);
				LOG(3,<<"BL Found... "<<paramFound<<"...LL="<<newL<<"..."<<endl;);
				LOG(3,<<"BL old... "<<oldBl<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			}			
			LOG(3,<<"BL old... "<<oldBl<<" BL done... "<<nodesV[i]->dis2father()<<"...LL="<<_treeLikelihood<<"..."<<endl;);
		}

		//alon
		cout<<endl<<"this is the TEMP branch length vector:   "<<endl;
		for(int i=0; i<nodesV.size(); i++)
		{
			cout<<endl<<nodesV[i]->dis2father();
			
		}
		//alon
	}

	//alon
	cout<<endl<<"this is the FINAL branch length vector:   "<<endl;
	for(int i=0; i<nodesV.size(); i++)
		{
			cout<<endl<<nodesV[i]->dis2father();
			
		}
	int a;
	cin>>a;
	//alon
	return _treeLikelihood;
}

MDOUBLE evalBranch::operator()(MDOUBLE x)
{
	_pNode->setDisToFather(x);
	MDOUBLE LL = _model.compLogLikelihood(*_pTree,_scIntron,_scExon);
	return -LL;
}

MDOUBLE evalCodonPositionBranch::operator()(MDOUBLE x)
{
	_pNode->setDisToFather(x);
	MDOUBLE LL = _model.compLogLikelihood(*_pTree,_sc);
	return -LL;
}

MDOUBLE evalCodonPositionBranchLight::operator()(MDOUBLE x)
{
	_pNode->setDisToFather(x);
	MDOUBLE LL = _model.compLogLikelihood(*_pTree,_sc);
	return -LL;
}


