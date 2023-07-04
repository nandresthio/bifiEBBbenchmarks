#ifndef SAMPLE_GENERATOR_CPP
#define SAMPLE_GENERATOR_CPP


#include "sample_generator.hpp"


SampleGenerator::SampleGenerator(Function* function, int randomSeed, bool printInfo):
	function_(function),
	printInfo_(printInfo){

	if(randomSeed != 0){
		mt19937 gen(randomSeed);
		randomGenerator_ = gen;
	}else{
		random_device rd;
		mt19937 gen(rd());
		randomGenerator_ = gen;
	}

	prePrint_ = "";
}

SampleGenerator::~SampleGenerator(){}

vector<VectorXd> SampleGenerator::randomSample(int sampleSize){
	vector<VectorXd> sample;
	for(int i = 0; i < sampleSize; i++){
		// Get the random number
		VectorXd newSample = VectorXd(function_->d_);
		for(int j = 0; j < function_->d_; j++){
			newSample(j) = function_->lowerBound_[j] + (function_->upperBound_[j] - function_->lowerBound_[j]) * ((double) randomGenerator_() / randomGenerator_.max());
		}
		sample.push_back(newSample);
	}
	return sample;

}

vector<VectorXd> SampleGenerator::randomLHS(int sampleSize){

	vector<VectorXd> sample;
	sample.reserve(sampleSize);
	
	// Need to define an ordering for each dimension
	int d = function_->d_;
	vector< vector<int> > orderings;
	orderings.reserve(d);
	for(int i = 0; i < d; i++){
		vector<int> ordering;
		ordering.reserve(sampleSize);
		for(int j = 0; j < sampleSize; j++){
			ordering.push_back(j);
		}
   		shuffleIntVector(ordering, randomGenerator_);
   		orderings.push_back(ordering);
	}
	// For each sample, get a random value within allowed range
	for(int i = 0; i < sampleSize; i++){
		VectorXd newSample = VectorXd(d);
		// For each dimension, get the range that is "allowed" and get a random number
		for(int j = 0; j < d; j++){
			int position = orderings[j][i];
			double curr_lower_bound = function_->lowerBound_[j] + position*(function_->upperBound_[j] - function_->lowerBound_[j])/sampleSize;
			double curr_upper_bound = function_->lowerBound_[j] + (position+1)*(function_->upperBound_[j] - function_->lowerBound_[j])/sampleSize;
			newSample(j) = (curr_lower_bound + (curr_upper_bound - curr_lower_bound) * ((double) randomGenerator_() / randomGenerator_.max()));
		}
		sample.push_back(newSample);
		if(printInfo_){printf("\r%sRandomLHS: Working on sample %d/%d", prePrint_.c_str(), (int)sample.size(), sampleSize);}
	}
	if(printInfo_){printf("\r%s                                            \r", prePrint_.c_str());}
	if(printInfo_ && prePrint_ != ""){printf("\r%s RandomLHS sample of size %d completed", prePrint_.c_str(), sampleSize);}
	return sample;
}



tuple<double, int, int> SampleGenerator::morrisMitchellCriterion(vector<VectorXd> &pointSet){
	// Simply find the min distance between the points
	double minDist = DBL_MAX;
	if((int)pointSet.size() <= 1){
		printf("This should not happen, asking for Morris Mitchell point using set of 1 point or less! Stopping now...");
		exit(0);
	}
	int firstIndex = -1;
	int secondIndex = -1;
	for(int i = 0; i < (int)pointSet.size(); i++){
		for(int j = i+1; j < (int)pointSet.size(); j++){
			double distance = dist(pointSet[i], pointSet[j]);
			if(distance < minDist){
				minDist = distance;
				firstIndex = i;
				secondIndex = j;
			}
		}
	}
	return make_tuple(minDist, firstIndex, secondIndex);
}


void SampleGenerator::morrisMitchellLocalOptimal(vector<VectorXd> &pointSet, int maxImprovements){
	int n = (int)pointSet.size();
	if(n == 0){
		printf("Error, asking to make sample locally optimal with Morris-Mitchel criterion, but given sample is empty! Stopping now...");
		exit(0);
	}
	int d = (int)pointSet[0].size();
	if(n == 1){
		return;
	}
	bool improved = true;
	tuple<double, int, int> info = morrisMitchellCriterion(pointSet);
	double distance = get<0>(info);
	int firstIndex = get<1>(info);
	int secondIndex = get<2>(info);
	double bestScore = distance;
	int firstSwapIndex;
	int secondSwapIndex;
	int dimensionSwap;
	int newFirstIndex;
	int newSecondIndex;
	int improvementNumber = 0;

	uniform_int_distribution<int>  distr(0, n);

	while(improved && (maxImprovements == 0 || improvementNumber < maxImprovements)){
		improved = false;
		for(int i = 0; i <= 1; i++){
			int index;
			if(i == 0){index = firstIndex;}
			else{index = secondIndex;}
			// Randomise "starting point" of points to check for swapping
			int start = distr(randomGenerator_);
			for(int j = 0; j < n; j++){
				int changingIndex = (start + j) % n;
				if(changingIndex == firstIndex || changingIndex == secondIndex){continue;}
				for(int k = 0; k < d; k++){
					double temp = pointSet[index](k);
					pointSet[index](k) = pointSet[changingIndex](k);
					pointSet[changingIndex](k) = temp;
					tuple<double, int, int> newInfo = morrisMitchellCriterion(pointSet);

					if(get<0>(newInfo) > bestScore){
						improved = true;
						bestScore = get<0>(newInfo);
						firstSwapIndex = index;
						secondSwapIndex = changingIndex;
						dimensionSwap = k;
						newFirstIndex = get<1>(newInfo);
						newSecondIndex = get<2>(newInfo);
					}

					temp = pointSet[index](k);
					pointSet[index](k) = pointSet[changingIndex](k);
					pointSet[changingIndex](k) = temp;

				}
				// If found a better point, stop search and swap
				if(improved){break;}
			}
			if(improved){break;}
		}
		if(improved){
			double temp = pointSet[firstSwapIndex][dimensionSwap];
			pointSet[firstSwapIndex][dimensionSwap] = pointSet[secondSwapIndex][dimensionSwap];
			pointSet[secondSwapIndex][dimensionSwap] = temp;
			firstIndex = newFirstIndex;
			secondIndex = newSecondIndex;
			improvementNumber++;
		}
		if(printInfo_){printf("\r%sMorris-Mitchell - Optimising sample, currently gone from %.4f to %.4f, improved %d times", prePrint_.c_str(), distance, bestScore, improvementNumber);}
	}
	if(printInfo_){printf("\r%s                                                                                                                      \r", prePrint_.c_str());}
	if(printInfo_ && prePrint_ != ""){printf("\r%sMorris-Mitchell optimisation completed, gone from %.4f to %.4f, improved %d times", prePrint_.c_str(), distance, bestScore, improvementNumber);}
}




vector<VectorXd> SampleGenerator::morrisMitchellSubset(vector<VectorXd> &pointSet, int subsetSize){
	int n = (int)pointSet.size();
	if(n < subsetSize){
		printf("Asked for a subset larger than the size of the original set! Initial set of size %d, wanted a set of size %d. Exiting now...\n", n, subsetSize);
		exit(0);
	}
	
	vector<int> ordering;
	ordering.reserve(pointSet.size());
	for(int i = 0; i < n; i++){
		ordering.push_back(i);
	}
	shuffleIntVector(ordering.begin(), ordering.end(), randomGenerator_);
	// Initial subset is read off from initial indexes in ordering
	vector<VectorXd> pointSubset;
	pointSubset.reserve(subsetSize);
	// Array of bools to keep track of what is in the subset, populate subset
	vector<bool> inSubset(n, false);
	for(int i = 0; i < subsetSize; i++){
		pointSubset.push_back(pointSet[ordering[i]]);
		inSubset[ordering[i]] = true;
	}
	// Calculate current score
	tuple<double, int, int> info = morrisMitchellCriterion(pointSubset);
	double distance = get<0>(info);
	int firstIndex = get<1>(info);
	int secondIndex = get<2>(info);
	double bestScore = distance;
	int inSwapIndex;
	int outSwapIndex;
	int newFirstIndex;
	int newSecondIndex;
	int improvementNumber = 0;

	// Iterate through entries and see if a single swap can improve the score
	bool improved = true;
	while(improved){
		improved = false;
		// Create temporary copy
		vector<VectorXd> pointSubsetCopy = pointSubset;
		// Want to check what happens when swapping out one of the two points which are closest to each other
		for(int i = 0; i <= 1; i++){
			int index;
			if(i == 0){index = firstIndex;}
			else{index = secondIndex;}
			// Iterate through the points not in the set
			for(int j = 0; j < n; j++){
				if(inSubset[j]){continue;}
				VectorXd temp = pointSubsetCopy[index];
				pointSubsetCopy[index] = pointSet[j];
				tuple<double, int, int> newInfo = morrisMitchellCriterion(pointSubsetCopy);
				if(get<0>(newInfo) > bestScore){
					improved = true;
					outSwapIndex = index;
					inSwapIndex = j;
					bestScore = get<0>(newInfo);
					newFirstIndex = get<1>(newInfo);
					newSecondIndex = get<2>(newInfo);
				}
				pointSubsetCopy[index] = temp;
			}

		}
		if(improved){
			pointSubset[outSwapIndex] = pointSet[inSwapIndex];
			firstIndex = newFirstIndex;
			secondIndex = newSecondIndex;
			inSubset[inSwapIndex] = true;
			inSubset[ordering[outSwapIndex]] = false;
			ordering[outSwapIndex] = inSwapIndex;
			improvementNumber++;

		}
		if(printInfo_){printf("\r%sMorris-Mitchell - Optimising sample subset, currently gone from %.4f to %.4f, improved %d times", prePrint_.c_str(), distance, bestScore, improvementNumber);}
	}
	if(printInfo_){printf("\r%s                                                                                                                      \r", prePrint_.c_str());}
	if(printInfo_){printf("\r%s Morris-Mitchell optimisation completed, gone from %.4f to %.4f, improved %d times", prePrint_.c_str(), distance, bestScore, improvementNumber);}
	return pointSubset;
}


void SampleGenerator::updateProblem(Function* function){
	function_ = function;
}










#endif