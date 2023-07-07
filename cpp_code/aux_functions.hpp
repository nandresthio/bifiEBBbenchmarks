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

#ifndef AUX_FUNCTIONS_HPP
#define AUX_FUNCTIONS_HPP

#include "libraries.hpp"
#include "sample_generator.hpp"
#include "functions.hpp"



// Implementation of the shuffle function which allows for replication across operating systems.
// This function takes a vector of integers and randomly reorders them.
void shuffleIntVector(vector<int> &vector, mt19937 &randomGenerator);

// Same as above, but this function takes in a vector of doubles.
void shuffleDoubleVector(vector<double> &vector, mt19937 &randomGenerator);

// Function which takes a point in the hypercube domain of a black-box function, and maps it
// to the unit hypercube [0,1]^d
void scalePoint(VectorXd &point, Function* function);

// Function which scales a vector of points by iterative calling the scalePoint function
void scalePoints(vector<VectorXd> &points, Function* function);

// Function which takes a point in the unit hypercube [0,1]^d and scales it to the hypercube
// domain of the supplied function.
void unscalePoint(VectorXd &point, Function* function);

// Function which unscales a vector of points by iteratively calling the unscale function
void unscalePoints(vector<VectorXd> &points, Function* function);

// Returns the probability density function of a standard normal distribution and point "value"
double normalPDF(double value);


// Returns the cumulative density function of a standard normal distribution and point "value"
double normalCDF(double value);

// A function which checks whether "point" is contained within the vector "points"
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
// Note that the actual radius taken is r^{1/d}
pair<vector<double>, vector<double> > calculateLocalCorrelations(BiFidelityFunction* function, double r, vector<double> pVals, vector<VectorXd> sample, vector<double> highSample, vector<double> lowSample);

// Useful function for debugging, takes care of printing point using the right formatting.
void printPoint(VectorXd point);


#endif