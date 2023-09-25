// ######################################################################################

// # Copyright 2023, Nicolau Andrés-Thió

// # Permission is hereby granted, free of charge, to any person obtaining a copy
// # of this software and associated documentation files (the "Software"), to deal
// # in the Software without restriction, including without limitation the rights
// # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// # copies of the Software, and to permit persons to whom the Software is
// # furnished to do so, subject to the following conditions:

// # The above copyright notice and this permission notice shall be included in all
// # copies or substantial portions of the Software.

// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// # SOFTWARE.

// ######################################################################################

#ifndef AUX_SOLVERS_CPP
#define AUX_SOLVERS_CPP

#include "aux_solvers.hpp"



AuxSolver::AuxSolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, 10000, min, randomSeed, printInfo){}


AuxSolver::AuxSolver(Function* function, int maxEval, bool min, int randomSeed, bool printInfo) :
	function_(function),
	min_(min),
	randomSeed_(randomSeed),
	printInfo_(printInfo),
	maxEval_(maxEval){

	if(randomSeed_ != 0){
		mt19937 gen(randomSeed);
		randomGenerator_ = gen;
	}else{
		random_device rd;
		mt19937 gen(rd());
		randomGenerator_ = gen;
	}

	sampleGenerator_ = new SampleGenerator(function_, randomSeed_);
	prePrint_ = "";
}

AuxSolver::~AuxSolver(){
	delete sampleGenerator_;
}


void AuxSolver::reseedRandom(int newSeed){
	randomSeed_ = newSeed;
	if(randomSeed_ != 0){
		mt19937 gen(randomSeed_);
		randomGenerator_ = gen;
	}else{
		random_device rd;
		mt19937 gen(rd());
		randomGenerator_ = gen;
	}
	delete sampleGenerator_;
	sampleGenerator_ = new SampleGenerator(function_, randomSeed_);
}

void AuxSolver::updateProblem(Function* function, bool min){
	function_ = function;
	min_ = min;
	// Adding logic for the starting generator seed to change as well
	// Only do so if it is not 0
	if(randomSeed_){randomSeed_++;}
	// Also need to update generator
	sampleGenerator_->updateProblem(function);
}

VectorXd AuxSolver::optimise(){
	printf("Should not be calling this AuxSolver base class to optimise!! Stopping now...\n");
	exit(0);
}

int AuxSolver::betterPoint(VectorXd point1, double val1, VectorXd point2, double val2){
	if(abs(val1 - val2) < TOL){return 0;}
	else if((min_ && val1 < val2) | (!min_ && val1 > val2)){return -1;}
	else{return 1;}
}

ARSsolver::ARSsolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, 10000, min, randomSeed, printInfo),
	numSearch_(10){
}

ARSsolver::ARSsolver(Function* function, int numSearch, int maxEval, bool min, int randomSeed, bool printInfo):
	AuxSolver(function, maxEval, min, randomSeed, printInfo),
	numSearch_(numSearch){
	}

ARSsolver::~ARSsolver(){
}

VectorXd ARSsolver::optimise(){
	int numSphere = numSearch_;
	vector<VectorXd> centers = sampleGenerator_->randomLHS(numSphere);
	return optimise(centers);
}

VectorXd ARSsolver::optimise(vector<VectorXd> centers){
	int numSphere = numSearch_;
	int maxEval = maxEval_;
	// If allocating some evaluations to local search, decrease the maximum number of iterations
	int d = function_->d_;
	double maxRadius = 1;
	vector<double> radii(numSphere, maxRadius);
	vector<double> vals;
	// Store best point so far
	double bestVal = function_->evaluate(centers[0]);
	VectorXd bestPoint = centers[0];
	vals.push_back(bestVal);
	for(int i = 1; i < numSphere; i++){
		vals.push_back(function_->evaluate(centers[i]));
		if(betterPoint(bestPoint, bestVal, centers[i], vals[i]) == 1){
			bestPoint = centers[i];
			bestVal = vals[i];
		}
	}
	int evals = 0;
	if(printInfo_){
		printf("\r%sRunning ARS, evaluations %d/%d, best point has value %.20f at point (", prePrint_.c_str(), evals, maxEval, bestVal);
		for(int i = 0; i < d; i++){
			printf("%.4f", bestPoint[i]);
			if(i < d-1){printf(", ");}
			else{printf(")");}
		}
	}
	// For each iteration, deal with each sphere
	while(evals < maxEval){
		for(int i = 0; i < numSphere; i++){
			// Did not read paper properly the first time! Correct implementation is calling the function below.
			VectorXd testPoint = findNewPoint(centers[i], radii[i]);
			double testPointValue = function_->evaluate(testPoint);

			if(betterPoint(centers[i], vals[i], testPoint, testPointValue) == 1){
				centers[i] = testPoint;
				vals[i] = testPointValue;
				if(betterPoint(bestPoint, bestVal, testPoint, testPointValue) == 1){
					bestPoint = testPoint;
					// f_min = vals[i];
					bestVal = testPointValue;
				}
				radii[i] = maxRadius;
			
			}else{
				radii[i] = radii[i] / 2;
				if(radii[i] < TOL){radii[i] = maxRadius;}
				
			}
			evals++;
			if(evals >= maxEval){break;}
		}
		if(printInfo_){
			printf("\r%sRunning ARS, evaluations %d/%d, best point has value %.20f at point (", prePrint_.c_str(), evals, maxEval, bestVal);
			for(int i = 0; i < d; i++){
				printf("%.4f", bestPoint[i]);
				if(i < d-1){printf(", ");}
				else{printf(")                         ");}
			}
		}
	}
	if(printInfo_){printf("\n");}
	return bestPoint;
}


VectorXd ARSsolver::dBallRandomSample(int d){
	// First define randoms
	normal_distribution<double> normalDis(0, 1);
    // Define direction
    VectorXd direction = VectorXd(d);
    for(int i = 0; i < d; i++){
    	direction(i) = normalDis(randomGenerator_);
    }
    double norm = direction.norm();
    // Define radius
    double radius = pow((double) randomGenerator_() / randomGenerator_.max(), 1.0/d);
    // double radius = 1;
    // Combine into a point
    for(int i = 0; i < d; i++){
    	double final_value = radius * direction(i)/norm;
    	direction(i) = final_value;	
    }
    return direction;
}

VectorXd ARSsolver::findNewPoint(VectorXd point, double radius){
	VectorXd newPoint(function_->d_);
	// Cycle through the entries, randomly choose a point
	for(int i = 0; i < function_->d_; i++){
		double minVal = max(function_->lowerBound_[i], point(i) - radius * (function_->upperBound_[i] - function_->lowerBound_[i]));
		double maxVal = min(function_->upperBound_[i], point(i) + radius * (function_->upperBound_[i] - function_->lowerBound_[i]));
		double range = maxVal - minVal;
		newPoint(i) = minVal + range * ((double) randomGenerator_() / randomGenerator_.max());
	}
	return newPoint;
}





#endif