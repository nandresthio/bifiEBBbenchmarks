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

#ifndef SAMPLE_GENERATOR_HPP
#define SAMPLE_GENERATOR_HPP

#include "libraries.hpp"
#include "aux_functions.hpp"
#include "functions.hpp"


// Class used to generate samples given a function instance (which contains dimension and bound information)
// If a seed other than 0 is specified, the behaviour is random; otherwise it is reproducible across operating systems
class SampleGenerator{
	public:

	SampleGenerator(Function* function, int randomSeed, bool printInfo = false);

	~SampleGenerator();

	// Generates a random sample of specified size; a check is performance to assert that no
	// to points are on top of one another
	vector<VectorXd> randomSample(int sampleSize);

	// Generates a Latin Hypercube Sampling plan of size n = sampleSize. This is done by 
	// breaking up each dimension into n cells, and placing n points so that no to points lie in the 
	// same cell of the same dimension. The location within the cell is uniformly randomly chosen.
	vector<VectorXd> randomLHS(int sampleSize);

	// Finds the morris-mitchell criterion of a set of points, namely it returns the minimum distance
	// between any pair of points, and the indices of the two closest points
	tuple<double, int, int> morrisMitchellCriterion(vector<VectorXd> &pointSet);

	// Function which makes a set of points locally optimal according to the morris-mitchell criterion.
	// Starting with a set of points, it swaps two coordinates of two points if doing so increases the 
	// minimum distance between any pair of points. This is done until such a swap no longer exists,
	// at which point the set is considered locally optimal.
	void morrisMitchellLocalOptimal(vector<VectorXd> &pointSet, int maxImprovements = 0);

	// Function which finds a locally optimal subset of a set of points. This is done by starting with a subset
	// of the points, and swapping a point inside the set with one outside the set as long as this 
	// increases the minimum distance between any pair of points in the subset. Once this can no longer be done,
	// the subset is said to be locally optimal.
	vector<VectorXd> morrisMitchellSubset(vector<VectorXd> &pointSet, int subsetSize);

	// Helper method which updates the function of the sample generator
	void updateProblem(Function* function);


	mt19937 randomGenerator_;
	Function* function_;
	bool printInfo_;

	string prePrint_;

};

#endif