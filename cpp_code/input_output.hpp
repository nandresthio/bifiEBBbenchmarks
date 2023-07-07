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

#ifndef INPUT_OUTPUT_HPP
#define INPUT_OUTPUT_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "functions.hpp"


// Function provided as an interface so that only an instance name needs to be provided,
// and an instantiation of the class is returned. 
// The names of functions which have been implemented can be found in the folder data/availableFunctions.
BiFidelityFunction* processFunctionName(string name, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0);

// If calling "proccessFunctionName" for a disturbance function, this function is called for further processing of the function name
BiFidelityFunction* processSpecialFunctionName(string name, bool knowOptVals = false, double knownFmin = 0, double knownFmax = 0);

// Function which takes as input a file name, a number of points and the dimension of those points
// and returns a vector containing the points. Note no error checking is performed.
vector<VectorXd> readPointsFromFile(string filename, int pointsNum, int dimension);

// Function which takes as input a set of points and the dimension of the points, and outputs them to a file
// specified by filename
void writePointsToFile(string filename, vector<VectorXd> points, int pointsNum, int dimension);

// Generates a high and low fidelity sample. A sample of size "lowFiBudget" (or "highFiBudget" if 'lowFiBudget' = 0)
// is chosen by finding a random LHS sample plan and doing local swaps to find a morris-mitchell locally optimal sample,
// that is a sample for which swapping the coordinate of two points does not increase the minimum distance between 
// any pair of points. If 'lowFiBudget' > 0, a subset of size 'highFiBudget' is chosen to also be morris-mitchell locally optimal.
// Note if the sample plan already exists (i.e. a file with the right name exists) this work is skipped.
// The sample is saved within the hypercube [0,1]^d, where d is the dimension of the specified function
void generateAndSaveSample(Function* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo);

// Either reads in or generates a sample of size lowFiBudget and a subset of size highFiBudget,
// or if lowFiBudget = 0 it generates or reads in a sample of size highFiBudget.
// The seed specifies what the random generator should be seeded with for reproducibility purposes, 
// as it might have already been generated and saved to a file.
// Note that if reading from a file, the assumption is that the sample lies within [0,1]^d, and the 
// sample is linearly transformed to lie inside the sample space of the specified function.
pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(Function* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo);


#endif