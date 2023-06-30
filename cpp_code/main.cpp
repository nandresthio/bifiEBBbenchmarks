#ifndef MAIN_CPP
#define MAIN_CPP

#include "libraries.hpp"
#include "input_output.hpp"

int main(int argc, char *argv[]){

	// The code in this file shows how the repo can be used
	
	// This line creates a bifidelity function class which can be queried for both high and low fidelity function values
	BiFidelityFunction* function = processFunctionName("ToalBranin0.50");
	
	// The dimension and sampling domain can be accessed as follows:
	printf("Working with a function of dimension %d, and sampling domain ", function->d_);
	for(int i = 0; i < function->d_; i++){
		printf("[%.2f,%.2f]", function->lowerBound_[i], function->upperBound_[i]);
		if(i < function->d_ - 1){printf("x");}
	}
	printf("\n");


	// The following line creates and saves a Morris-Mitchel optimal Latin Hypercube Sampling (LHS) plan of the dimension of the 
	// given function (in this case, dimension d = 2) of size 10. Note the sampling plan is scaled and saved in the hypercube [0,1]^d
	// instead of the actual sampling domain of the function.
	// Note that as this sample already exists in the data folder, when this runs the existing sample is checked for optimality.
	printf("\nCreate and save a sample of size 10:\n");
	int sampleSize = 10;
	int subsetSampleSize = 0;
	int seed = 1;
	bool printInfo = true;
	generateAndSaveSample(function, sampleSize, subsetSampleSize, seed, printInfo);

	// The following both finds a Morris-Mitchel optimal LHS plan of size 10, as well as a Morris-Mitchel sample subset of size 5.
	// Note again that as these already exist, only optimality is checked.
	printf("\nCreate and save a sample of size 10, and a subset sample of size 5:\n");
	sampleSize = 10;
	subsetSampleSize = 5;
	generateAndSaveSample(function, sampleSize, subsetSampleSize, seed, printInfo);

	// The following reads in or generates a high and low fidelity sample of the specified sizes, where the high fidelity sample is a subset
	// of the low fidelity sample. If the low fidelity sample size is set to 0, the high fidelity sample is taken as a 
	// Morris-Mitchel LHS sample.
	// If the sampling plans exist in the data/samplePlans folder, it is read in; otherwise it is generated.
	// Note that if the sample is read in, it is scaled to the sampling domain of the specified function.
	printf("\nRead in a high and low fidelity sample of sizes 5 and 10 respectively:\n");
	int highFiBudget = 5;
	int lowFiBudget = 10;
	pair<vector<VectorXd>, vector<VectorXd> > points = readInOrGenerateInitialSample(function, highFiBudget, lowFiBudget, seed, printInfo);
	


	// The following lines print the high fidelity function value at the high fidelity function points,
	// and the low fidelity function values at the low fidelity function points
	printf("\nShow high and low objective function values at the respective sample points:\n");
	vector<VectorXd> highFiSamples = points.first;
	vector<VectorXd> lowFiSamples = points.second;
	for(int i = 0; i < (int)lowFiSamples.size(); i++){
		printPoint(lowFiSamples[i]);
		printf(": Low fidelity value of %.2f\n", function->evaluateLow(lowFiSamples[i]));
	}

	printf("\n");

	for(int i = 0; i < (int)highFiSamples.size(); i++){
		printPoint(highFiSamples[i]);
		printf(": High fidelity value of %.2f\n", function->evaluate(highFiSamples[i]));
	}
	

	return 0;
}



#endif