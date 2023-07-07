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

#include <Rcpp.h>

#include "libraries.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"
#include "aux_solvers.hpp"
#include "input_output.hpp"
#include "aux_functions.hpp"



// This file contains interface functions to call the C++ functions implemented in this package by an R script
// using Rcpp.

// Function which transfors point with data type VectorXd to a point with data type vector<double> 
vector<double> vectorXdtoVectorDouble(VectorXd point){
	vector<double> newPoint;
	for(int i = 0; i < point.size(); i++){
		newPoint.push_back(point(i));
	}
	return newPoint;
}

// Function which transforms a vector of VectorXd points to a vector of vector<double> points
vector< vector<double> > multipleVectorXdtoVectorDouble(vector<VectorXd> points){
	vector< vector<double> > newPoints;
	for(int i = 0; i < (int)points.size(); i++){
		newPoints.push_back(vectorXdtoVectorDouble(points[i]));
	}
	return newPoints;
}

// Function which transforms a point with data type vector<double> to a point with data type VectorXd
VectorXd vectorDoubletoVectorXd(vector<double> point){
	VectorXd newPoint(point.size());
	for(int i = 0; i < (int)point.size(); i++){
		newPoint(i) = point[i];
	}
	return newPoint;
}

// Function which transforms a vector of vector<double> points to a vector of VectorXd points
vector< VectorXd > multipleVectorDoubletoVectorXd(vector< vector<double> > points){
	vector< VectorXd > newPoints;
	for(int i = 0; i < (int)points.size(); i++){
		newPoints.push_back(vectorDoubletoVectorXd(points[i]));
	}
	return newPoints;
}


// Function which either returns the f_h, f_l or f_h - f_l function value (specified by a level of 0, 1 or 2 respectively).
// The function sampled is specified by its function name. For fast initialisation of the function, the maximum and minimum function values
// can be passed so that they don't need to be found for disturbance based functions
// [[Rcpp::export]]
double sampleFunction(string functionName, vector<double> &vectorPoint, int level = 0, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	// Super clunky, but initialise function, then sample it, then return the value
	VectorXd point = Map<VectorXd, Unaligned>(vectorPoint.data(), vectorPoint.size());
	for(int i = 0; i < point.size(); i++){
	}
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	double val;
	if(level == 0){
		val = function->evaluate(point);
	}
	else if(level == 1){
		val = function->evaluateLow(point);
	}
	else if(level == 2){
		val = function->evaluate(point) - function->evaluateLow(point);
	}else{
		printf("Problem with level chosen!\n");
		return 0;
	}
	delete function;
	return val;
}

// Function which returns the dimension of the instance specified by functionName
// [[Rcpp::export]]
int functionDimension(string functionName, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	int dimension = function->d_;
	delete function;
	return dimension;
}

// Function which returns the lower bound of the domain of the instance specified by functionName
// [[Rcpp::export]]
vector<double> functionLowerBound(string functionName, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	vector<double> bound = function->lowerBound_;
	delete function;
	return bound;
}

// Function which returns the upper bound of the domain of the instance specified by functionName
// [[Rcpp::export]]
vector<double> functionUpperBound(string functionName, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	vector<double> bound = function->upperBound_;
	delete function;
	return bound;
}

// Function which returns the minimum and maximum function values the instance specified by functionName
// [[Rcpp::export]]
vector<double> functionMinMax(string functionName, int seed, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	ARSsolver* auxSolver = new ARSsolver(function, 10, 5000, false, seed, false);
	VectorXd best = auxSolver->optimise();
	double fMax = function->evaluate(best);
	auxSolver->updateProblem(function, true);
	best = auxSolver->optimise();
	double fMin = function->evaluate(best);
	delete auxSolver;
	delete function;
	vector<double> range;
	range.push_back(fMin);
	range.push_back(fMax);
	return range;
}

// Function which returns feature values of the instance specified by function name, specifically the features CC, RRMSE, and LCC features
// [[Rcpp::export]]
vector<double> functionBasicFeatures(string functionName, int seed, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	// printf("In basic features, got %d %.5f %.5f\n", knowOptVals, knownFmin, knownFmax);
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	return calculateFunctionFeatures(function, 5000, seed);
}

// Same as above, but the locations of the sample are given as part of the input
// [[Rcpp::export]]
vector<double> functionBasicFeaturesWithSample(string functionName, vector< vector<double> > sample, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	// Need to change vector<double> to VectorXd
	vector< VectorXd > convertedSample = multipleVectorDoubletoVectorXd(sample);
	return calculateFunctionFeatures(function, 0, 0, 0.2, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975}, convertedSample);
}

// Same as above, but the locations of the sample as well as the function values are given as part of the input
// [[Rcpp::export]]
vector<double> functionBasicFeaturesWithSampleAndVals(string functionName, vector< vector<double> > sample, vector<double> fHigh, vector<double> fLow, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	// Need to change vector<double> to VectorXd
	vector< VectorXd > convertedSample = multipleVectorDoubletoVectorXd(sample);
	return calculateFunctionFeatures(function, 0, 0, 0.2, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975}, convertedSample, fHigh, fLow);
}

// Returns a sample of size lowFiBudget and a subset sample of size highFiBudget fit to the sample space of the specified functionName
// [[Rcpp::export]]
vector< vector< vector<double> > > functionSample(string functionName, int seed, int highFiBudget, int lowFiBudget, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0){
	BiFidelityFunction* function = processFunctionName(functionName, knowOptVals, knownFmin, knownFmax);
	// Here want to access the sample from the appropriate file
	string samplePlanFilename = "data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
	vector< VectorXd > sampledPointsLow = readPointsFromFile(samplePlanFilename, lowFiBudget, function->d_);
	
	samplePlanFilename = "data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
	vector< VectorXd > sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
	if((int)sampledPoints.size() == 0){
		SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, false);
		sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
		delete sampleGenerator;
	}


	unscalePoints(sampledPointsLow, function);
	unscalePoints(sampledPoints, function);

	// Now need to change VectorXd to vector<double>
	vector< vector<double> > sampledPointsLowVector = multipleVectorXdtoVectorDouble(sampledPointsLow);
	vector< vector<double> > sampledPointsVector = multipleVectorXdtoVectorDouble(sampledPoints);
	vector< vector< vector<double> > > points;
	points.push_back(sampledPointsLowVector);
	points.push_back(sampledPointsVector);
	return points;
}



