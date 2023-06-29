#ifndef AUX_SOLVERS_HPP
#define AUX_SOLVERS_HPP

#include "libraries.hpp"
#include "aux_functions.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"


class AuxSolver{
	public: 

	AuxSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~AuxSolver();

	// Update the problem which ARS will be run on.
	void updateProblem(Function* function, bool min = true);

	// Reseed the random generator.
	void reseedRandom(int newSeed);

	int betterPoint(VectorXd point1, double val1, VectorXd point2, double val2);

	virtual VectorXd optimise();

	// VectorXd findGoodFirstPoint();

	Function* function_;						// The function which ARS will optimise.
	SampleGenerator* sampleGenerator_;			// Random generator which creates initial points inside sample space
	bool min_;									// Whether this is an optimisation or a maximisation problem
	int randomSeed_;							// Random seed used for reproducibility
	mt19937 randomGenerator_;					// Random number generator
	bool printInfo_;							// Whether to print the information as the optimisation is running
	string prePrint_;							// Used to add a message to the print statement when printInfo_ = true
	
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

	VectorXd optimise(vector<VectorXd> centers);

	// Returns a point within a hypersphere of dimension d and radius 1.
	VectorXd dBallRandomSample(int d);

	// Returns a point so that the sup norm of the new point and point is less than or equal to radius,
	// and the new point lies inside the sample space of the function.
	VectorXd findNewPoint(VectorXd point, double radius);

	int numSearch_;							// Number of spheres employed by the search
	int maxEval_;							// Number of evaluations the search is allowed

	

};



#endif