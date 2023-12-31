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

#ifndef AUX_FUNCTIONS_CPP
#define AUX_FUNCTIONS_CPP

#include "aux_functions.hpp"


void shuffleIntVector(vector<int> &inputVector, mt19937 &randomGenerator){
	for(int i = (int)inputVector.size() - 1; i > 0; i--){
		uniform_int_distribution<int> uniform(0,i);
		int temp = inputVector[i];
		int newIndex = uniform(randomGenerator);
		inputVector[i] = inputVector[newIndex];
		inputVector[newIndex] = temp;
	}
}

void shuffleDoubleVector(vector<double> &inputVector, mt19937 &randomGenerator){
	vector<double> copy = inputVector;
	vector<int> indices;
	indices.reserve((int)inputVector.size());
	for(int i = 0; i < (int)inputVector.size(); i++){
		indices.push_back(i);
	}
	shuffleIntVector(indices, randomGenerator);
	for(int i = 0; i < (int)inputVector.size(); i++){
		inputVector[i] = copy[i];
	}
}

void scalePoint(VectorXd &point, Function* function){
	for(int i = 0; i < (int)point.size(); i++){
		point(i) = (point(i) - function->lowerBound_[i]) / (function->upperBound_[i] - function->lowerBound_[i]);
	}
}

void scalePoints(vector<VectorXd> &points, Function* function){
	for(int i = 0; i < (int)points.size(); i++){
		scalePoint(points[i], function);
	}
}

void unscalePoint(VectorXd &point, Function* function){
	for(int i = 0; i < (int)point.size(); i++){
		point(i) = point(i) * (function->upperBound_[i] - function->lowerBound_[i]) + function->lowerBound_[i];
	}
}

void unscalePoints(vector<VectorXd> &points, Function* function){
	for(int i = 0; i < (int)points.size(); i++){
		unscalePoint(points[i], function);
	}
}

double normalPDF(double value){
   return (1/sqrt(2*M_PI)) * exp(-1/2 * value * value);
}


double normalCDF(double value){
   return 0.5 * erfc(-value * M_SQRT1_2);
}

bool checkIfRepeated(vector<VectorXd> &points, VectorXd point){
	for(int i = 0; i < (int)points.size(); i++){
		if((point - points[i]).norm() < TOL){return true;}
	}
	return false;
}

double dist(VectorXd &point1, VectorXd &point2){
	if(point1.size() != point2.size()){
		printf("Comparing distance between two points of different dimensions!\n");
		exit(0);
	}
	VectorXd temp = point1 - point2;
	return temp.norm();
}


double weightedCorrelationCoefficient(vector<double> &dataSet1, vector<double> &dataSet2, vector<double> &weights, bool square){
	// First check the lengths are all correct
	if(dataSet1.size() != dataSet2.size() || dataSet1.size() != weights.size()){
		printf("Incorrect usage of weighted correlation, inputed vectors are of different length! Exiting now...\n");
		exit(0);
	}
	int n = (int)dataSet1.size();
	if(n == 1){return 1;}
	// Apply formula!
	double sumWeights = 0.0;
	double nominatorMean1 = 0.0;
	double nominatorMean2 = 0.0;
	for(int i = 0; i < n; i++){
		sumWeights += weights[i];
		nominatorMean1 += weights[i] * dataSet1[i];
		nominatorMean2 += weights[i] * dataSet2[i];
	}
	double mean1 = nominatorMean1 / sumWeights;
	double mean2 = nominatorMean2 / sumWeights;

	double nominatorVar1 = 0.0;
	double nominatorVar2 = 0.0;
	for(int i = 0; i < n; i++){
		nominatorVar1 += weights[i] * pow(dataSet1[i] - mean1, 2);
		nominatorVar2 += weights[i] * pow(dataSet2[i] - mean2, 2);
	}
	double var1 = nominatorVar1 / sumWeights;
	double var2 = nominatorVar2 / sumWeights;
	
	// If both are flat lines, perfect correlation
	if(abs(var1) < TOL && abs(var2) < TOL){return 1.0;}
	// If only one of them is a flat line, not correlated at all
	if(abs(var1) < TOL || abs(var2) < TOL){return 0.0;}

	double nominatorCorr = 0.0;
	for(int i = 0; i < n; i++){
		nominatorCorr += weights[i] * (dataSet1[i] - mean1) * (dataSet2[i] - mean2);
	}
	if(square){
		return pow(nominatorCorr / (sumWeights * sqrt(var1) * sqrt(var2)), 2);
	}else{
		return nominatorCorr / (sumWeights * sqrt(var1) * sqrt(var2));
	}
	
}

double relativeRootMeanSquaredError(vector<double> &dataSet1, vector<double> &dataSet2){
	if((int)dataSet1.size() != (int)dataSet2.size()){
		printf("Incorrect usage of RMSE, inputed vectors are of different length! Exiting now...\n");
		exit(0);
	}
	double min = DBL_MIN, max = -DBL_MAX;
	int n = (int)dataSet1.size();
	double total = 0.0;
	for(int i = 0; i < n; i++){
		total += pow(dataSet1[i] - dataSet2[i], 2);
		if(dataSet1[i] < min){min = dataSet1[i];}
		if(dataSet1[i] > max){max = dataSet1[i];}	
	}

	total = total / n;
	total = sqrt(total);
	if(abs(min - max) > TOL){total = total / (max - min);}
	return total;
}

vector<double> calculateFunctionFeatures(BiFidelityFunction* function, int sampleSize, int seed, double r, vector<double> pVals, vector<VectorXd> givenSample, vector<double> highFiVals, vector<double> lowFiVals){
	// Work out sample size
	int size;
	vector<VectorXd> sample;
	vector<double> highSample;
	vector<double> lowSample;

	if(givenSample.size() == 0){
		// size = sampleSizeMult * function->d_;
		size = sampleSize;
		// Will need a generator
		SampleGenerator* generator = new SampleGenerator(function, seed, false);
		// Get the sample
		sample = generator->randomLHS(size);
		delete generator;
		highSample = function->evaluateMany(sample);
		lowSample = function->evaluateManyLow(sample);
	
	}else{
		size = givenSample.size();
		sample = givenSample;

		if((int)highFiVals.size() != size){
			
			highSample = function->evaluateMany(sample);
		}else{
			highSample = highFiVals;
		}

		if((int)lowFiVals.size() != size){
			lowSample = function->evaluateManyLow(sample);
		}else{
			lowSample = lowFiVals;
		}	
	}


	// Calculate features
	// Define all weights as 1 to get normal correlation coefficient
	vector<double> weightsNormal(size, 1);
	double correlationCoefficient = weightedCorrelationCoefficient(highSample, lowSample, weightsNormal);
	// double sampleCorrelation = weightedCorrelationCoefficient(highSample, lowSample, weightsNormal, false);
	double relativeError = relativeRootMeanSquaredError(highSample, lowSample);
	// Calculate LCC values
	pair< vector<double>, vector<double> > localCorrelations = calculateLocalCorrelations(function, r, pVals, sample, highSample, lowSample);
	// Store features
	int n = (int)localCorrelations.first.size();
	vector<double> results(2 + 2*n, 0.0);
	results[0] = correlationCoefficient;
	// results[1] = sampleCorrelation;
	results[1] = relativeError;
	for(int i = 0; i < n; i++){
		results[2 + i] = localCorrelations.first[i];
		results[2 + n + i] = localCorrelations.second[i];
	}
	
	// int n = (int)localCorrelations.first.size();
	// vector<double> results(2 + n, 0.0);
	// results[0] = correlationCoefficient;
	// results[1] = relativeError;
	// for(int i = 0; i < n; i++){
	// 	results[2 + i] = localCorrelations.first[i];
	// }


	return results;
}

pair<vector<double>, vector<double> > calculateLocalCorrelations(BiFidelityFunction* function, double r, vector<double> pVals, vector<VectorXd> sample, vector<double> highSample, vector<double> lowSample){
	int sampleSize = sample.size();
	if(sample.size() != highSample.size() || sample.size() != lowSample.size()){
		printf("Weird, calling calculate local correlations with different sized vectors for high, low and sample. Stopping now...\n");
		exit(0);
	}

	// So seems like the first thing to do is to scale the locations so that different dimension ranges don't mess with distances
	scalePoints(sample, function);
	// Now the points lie in [0,1]^d, so the radius of the hyperball with center at (0.5, ... , 0.5) 
	// which contains the domain has a radius of 0.5 * sqrt(d)
	// The hyperball which has a volume of size r relative to the original hyperball has a radius 
	// of r^{1/d} * sqrt(d)
	double maxDistRelativeVolume =  pow(r, 1.0 / function->d_) * sqrt(function->d_);

	// Also will use max distance and return both
	// Since scaled, this will be sqrt d, times the ratio
	double maxDist = r * sqrt(function->d_);

	// Calculate distance for which neighbourhood applies
	// double maxDist = 0.0;
	// for(int i = 0; i < function->d_; i++){
	// 	maxDist += pow(function->upperBound_[i] - function->lowerBound_[i], 2);
	// }
	// maxDist = pow(r, 1.0 / function->d_) * sqrt(maxDist);
	// Cycle through each sample and calculate and store local correlation
	vector<double> localCorrValues(sampleSize, 0.0);
	vector<double> weights;
	vector<double> localHighSample;
	vector<double> localLowSample;
	weights.reserve(sampleSize);
	localHighSample.reserve(sampleSize);
	localLowSample.reserve(sampleSize);

	vector<double> localCorrValuesRelative(sampleSize, 0.0);
	vector<double> weightsRelative;
	vector<double> localHighSampleRelative;
	vector<double> localLowSampleRelative;
	weightsRelative.reserve(sampleSize);
	localHighSampleRelative.reserve(sampleSize);
	localLowSampleRelative.reserve(sampleSize);
	for(int i = 0; i < sampleSize; i++){
		VectorXd point = sample[i];
		// Cycle through all points, if closer than max dist, add to local vector
		for(int j = 0; j < sampleSize; j++){
			double localDist = dist(point, sample[j]);
			if(localDist <= maxDist){
				localHighSample.push_back(highSample[j]);
				localLowSample.push_back(lowSample[j]);
				weights.push_back(1.0 - localDist / maxDist);
			}
			if(localDist <= maxDistRelativeVolume){
				localHighSampleRelative.push_back(highSample[j]);
				localLowSampleRelative.push_back(lowSample[j]);
				weightsRelative.push_back(1.0 - localDist / maxDistRelativeVolume);
			}
			
			
		}
		// Calculate local correlation
		localCorrValues[i] = weightedCorrelationCoefficient(localLowSample, localHighSample, weights);
		localCorrValuesRelative[i] = weightedCorrelationCoefficient(localLowSampleRelative, localHighSampleRelative, weightsRelative);
		// Done! Clear vectors for next iteration
		weights.clear();
		localHighSample.clear();
		localLowSample.clear();
	}
	// Now process values to get features
	vector<double> localCorrs((int)pVals.size(), 0.0);
	vector<double> localCorrsRelative((int)pVals.size(), 0.0);
	
	// First LCC^r_p values
	// Add to count if larger than or equal to cut off
	for(int i = 0; i < sampleSize; i++){
		for(int j = 0; j < (int)pVals.size(); j++){
			if(localCorrValues[i] >= pVals[j]){localCorrs[j]++;}
			if(localCorrValuesRelative[i] >= pVals[j]){localCorrsRelative[j]++;}
		}
	}
	// Divide by total and done
	for(int j = 0; j < (int)pVals.size(); j++){localCorrs[j] = localCorrs[j] / sampleSize;}
	for(int j = 0; j < (int)pVals.size(); j++){localCorrsRelative[j] = localCorrsRelative[j] / sampleSize;}
		

	// Now calculate LCC^r_mean, LCC^r_sd and LCC^r_coeff
	double lccMean = 0;
	double lccMeanRelative = 0;
	for(int i = 0; i < sampleSize; i++){
		lccMean += localCorrValues[i];
		lccMeanRelative += localCorrValuesRelative[i];
	}
	lccMean = lccMean / sampleSize;
	lccMeanRelative = lccMeanRelative / sampleSize;

	double lccSD = 0;
	double lccSDRelative = 0;
	for(int i = 0; i < sampleSize; i++){
		lccSD += pow(localCorrValues[i] - lccMean, 2);
		lccSDRelative += pow(localCorrValuesRelative[i] - lccMeanRelative, 2);
		
	}
	lccSD = sqrt(lccSD / (sampleSize - 1));
	lccSDRelative = sqrt(lccSDRelative / (sampleSize - 1));


	double lccCoeff = lccSD;
	double lccCoeffRelative = lccSDRelative;
	
	if(lccMean > TOL){lccCoeff = lccCoeff / lccMean;}
	if(lccMeanRelative > TOL){lccCoeffRelative = lccCoeffRelative / lccMeanRelative;}
	
	localCorrs.push_back(lccMean);
	localCorrs.push_back(lccSD);
	localCorrs.push_back(lccCoeff);

	localCorrsRelative.push_back(lccMeanRelative);
	localCorrsRelative.push_back(lccSDRelative);
	localCorrsRelative.push_back(lccCoeffRelative);

	// Unscale points before I forget and I use them in the future
	unscalePoints(sample, function);
	
	return make_pair(localCorrs, localCorrsRelative);
}

void printPoint(VectorXd point){
	printf("(");
	for(int i = 0; i < point.size(); i++){
		printf("%.4f",point(i));
		if(i < point.size()-1){printf(", ");}
		else{printf(")");}
	}
}

#endif