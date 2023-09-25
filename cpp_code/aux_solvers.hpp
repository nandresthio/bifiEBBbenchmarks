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

#ifndef AUX_SOLVERS_HPP
#define AUX_SOLVERS_HPP

#include "libraries.hpp"
#include "aux_functions.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"

// Parent class of an black-box solver used for auxiliary optimisation such as training a model or finding the next sample point
// Takes in the function to be optimised, whether the problem is a minimisation or maximisation, the random seed for 
// reproducibiltiy purposes, and whether to print information as it runs
class AuxSolver{
	public: 

	AuxSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	AuxSolver(Function* function, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~AuxSolver();

	// Update the problem which the solver will be run on; both the function and whether it will be a max or min problem
	void updateProblem(Function* function, bool min = true);

	// Reseed the random generator.
	void reseedRandom(int newSeed);

	// Check whether a potential point is better or worse; implemented so that derived class can override it if relevant
	// Here, if two points have a near enough objective function value the function returns 0, it returns -1 if point 1
	// is better, and it returns 1 if point 2 is better.
	virtual int betterPoint(VectorXd point1, double val1, VectorXd point2, double val2);

	// Main function call to be overriden by derived classes
	virtual VectorXd optimise();

	Function* function_;						// The function which ARS will optimise.
	SampleGenerator* sampleGenerator_;			// Random generator which creates initial points inside sample space
	bool min_;									// Whether this is an optimisation or a maximisation problem
	int randomSeed_;							// Random seed used for reproducibility
	mt19937 randomGenerator_;					// Random number generator
	bool printInfo_;							// Whether to print the information as the optimisation is running
	string prePrint_;							// Used to add a message to the print statement when printInfo_ = true

	int maxEval_;								// Number of evaluations the search is allowed
	
};




// Implementation of Accelerated Random Search (ARS), a black box solver used for the auxiliary optimisation problems
// when constructing surrogate models. For details, consult
// Appel MJ, Labarre R and Radulovic D (2004): "On accelerated random search"
class ARSsolver : public AuxSolver{
	public:

	ARSsolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	ARSsolver(Function* function, int numSphere, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~ARSsolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Method to be called if want to specify the starting centers of the balls. 
	// Called by the main optimise function after the calls are randomly initialised.
	VectorXd optimise(vector<VectorXd> centers);

	// Returns a point within a hypersphere of dimension d and radius 1.
	VectorXd dBallRandomSample(int d);

	// Returns a point so that the sup norm of the new point and point is less than or equal to radius,
	// and the new point lies inside the sample space of the function.
	VectorXd findNewPoint(VectorXd point, double radius);

	int numSearch_;							// Number of spheres employed by the search
	
};



#endif