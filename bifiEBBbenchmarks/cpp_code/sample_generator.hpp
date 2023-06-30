#ifndef SAMPLE_GENERATOR_HPP
#define SAMPLE_GENERATOR_HPP

#include "libraries.hpp"
#include "aux_functions.hpp"
#include "functions.hpp"

class SampleGenerator{
	public:

	SampleGenerator(Function* function, int randomSeed, bool printInfo = false);

	~SampleGenerator();

	vector<VectorXd> randomSample(int sampleSize);

	vector<VectorXd> randomLHS(int sampleSize);

	tuple<double, int, int> morrisMitchellCriterion(vector<VectorXd> &pointSet);

	void morrisMitchellLocalOptimal(vector<VectorXd> &pointSet, int maxImprovements = 0);

	vector<VectorXd> morrisMitchellSubset(vector<VectorXd> &pointSet, int subsetSize);

	void updateProblem(Function* function);


	mt19937 randomGenerator_;
	Function* function_;
	bool printInfo_;

	string prePrint_;

};

#endif