#ifndef INPUT_OUTPUT_HPP
#define INPUT_OUTPUT_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "functions.hpp"

// Functions used to interpret the information passed as input when calling the project.

// Function which takes a function name and outputs the correct bifidelity source function. For more information on the available functions,
// consult the functions.hpp file.
BiFidelityFunction* processFunctionName(string name, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0);

// Need further processing when specifying a COCO function, this function deals with that.
BiFidelityFunction* processSpecialFunctionName(string name, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0);

// Function which takes in an experiment specification as a string and processes the information. It then calls executeExperiment to run the specifications.
// r and pVals are specifications used when calculating the LCC feature. 
// void processExperiment(string outputFilename, string instructionLine, double r, vector<double> pVals);
void processExperiment(string outputFilename, string instructionLine, int problemType);

vector<VectorXd> readPointsFromFile(string filename, int pointsNum, int dimension);


void writePointsToFile(string filename, vector<VectorXd> points, int pointsNum, int dimension);


// Executes the experiment as speficied by "processExperiment". 
// double executeExperiment(string filename, string functionName, string technique, int highFiBudget, int lowFiBudget, int seed, double r, vector<double> pVals,
// 							bool printSolverInfo = true, bool printAllInfo = true, int testSampleSize = 1000, int auxSolverIterations = 1000);
vector<double> executeExperiment(string filename, string functionName, string technique, int problemType, int budget, double lowFiBudgetORcostRatio, int seed,
									bool printInfo = true, bool printSolverInfo = false, int testSampleSize = 5000, int auxSolverIterations = 1000);


// Generates a high and low fidelity sample. A sample of size "lowFiBudget" (or "highFiBudget" if 'lowFiBudget' = 0) is chosen by finding a random LHS sample plan
// and doing local swaps to find a morris-mitchell locally optimal sample.
// If 'lowFiBudget' > 0, a subset of size 'highFiBudget' is chosen to also be morris-mitchell locally optimal.
// Note if the sample plan already exists (i.e. a file with the right name exists) this work is skipped.
void generateAndSaveSample(Function* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo);


// Either reads in or generates initial sample. If the sample cannot be read in from the file, a random LHS sample is generated instead.
pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(Function* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo);


void plotAndAnalysieFunction(string functionName, int highFiBudget, int lowFiBudget, int seed, bool printInfo, bool scaling = true);




#endif