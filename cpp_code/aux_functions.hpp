#ifndef AUX_FUNCTIONS_HPP
#define AUX_FUNCTIONS_HPP

#include "libraries.hpp"
#include "sample_generator.hpp"
#include "functions.hpp"

// Auxiliary functions to be used by different methods troughout the project.

void shuffleIntVector(vector<int> &vector, mt19937 randomGenerator);

void shuffleDoubleVector(vector<double> &vector, mt19937 randomGenerator);




void scalePoint(VectorXd &point, Function* function);

void scalePoints(vector<VectorXd> &points, Function* function);

void unscalePoint(VectorXd &point, Function* function);

void unscalePoints(vector<VectorXd> &points, Function* function);



double normalPDF(double value);

double normalCDF(double value);

bool checkIfRepeated(vector<VectorXd> &points, VectorXd point);

// Returns the Euclidean distance between two points
double dist(VectorXd &point1, VectorXd &point2);

// Weighted correlation of two datasets given a set of weights. For more information, consult
// Andres-Thio N, Munoz MA, Smith-Miles K (2022): "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances"
double weightedCorrelationCoefficient(vector<double> &dataSet1, vector<double> &dataSet2, vector<double> &weights, bool square = true);

// Relative Root Mean Squared Error (RRMSE) used to calculate accuracy of trained surrogate model,
// as well as as a feature of a bi-fidelity source. For more information, consult
// Andres-Thio N, Munoz MA, Smith-Miles K (2022): "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances"
double relativeRootMeanSquaredError(vector<double> &dataSet1, vector<double> &dataSet2);

// Calculates the Correlation Coefficient (CC), Relative Root Mean Squared Error (RRMSE) and Local Correlation Coefficients (LCC) of a bi-fidelity source.
// Does so by taking a sample of size sampleSizeMult * d, where d is the dimension of the bi-fidelity source. For more information, consult
// Andres-Thio N, Munoz MA, Smith-Miles K (2022): "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances"
vector<double> calculateFunctionFeatures(BiFidelityFunction* function, int sampleSize, int seed, double r = 0.2, vector<double> pVals = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975}, vector<VectorXd> givenSample = {}, vector<double> highFiVals = {}, vector<double> lowFiVals = {});

// For a bi-fidelity source, this function calculates LCC^r_p, LCC^r_sd and LCC^r_coeff for given r and a set of p values.
// Does so by taking a sample of size sampleSizeMult * d, where d is the dimension of the bi-fidelity source, and calculating 
// the weighted correlation at each of those points. For more information, consult
// Andres-Thio N, Munoz MA, Smith-Miles K (2022): "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances"
pair<vector<double>, vector<double> > calculateLocalCorrelations(BiFidelityFunction* function, double r, vector<double> pVals, vector<VectorXd> sample, vector<double> highSample, vector<double> lowSample);

// Useful function for debugging, takes care of printing point using the right formatting.
void printPoint(VectorXd point);


#endif