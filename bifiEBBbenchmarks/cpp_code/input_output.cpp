#ifndef INPUT_OUTPUT_CPP
#define INPUT_OUTPUT_CPP

#include "input_output.hpp"

BiFidelityFunction* processFunctionName(string name, bool knowOptVals, double knownFmin, double knownFmax){
	if(name.compare("LiuPedagogical") == 0){
		return new LiuPedagogicalFunction();
	
	}else if(name.compare("ShiGramacyLee") == 0){
		return new ShiGramacyLeeFunction();

	}else if(name.compare("ShiCurrinSin") == 0){
		return new ShiCurrinSinFunction();

	}else if(name.compare("ShiHolsclaw") == 0){
		return new ShiHolsclawFunction();

	}else if(name.compare("ShiSantner") == 0){
		return new ShiSantnerFunction();

	}else if(name.compare("LiuBranin") == 0){
		return new LiuBraninFunction();

	}else if(name.compare("ShiBranin") == 0){
		return new ShiBraninFunction();

	}else if(name.compare("ShiNumberSix") == 0){
		return new ShiNumberSixFunction();

	}else if(name.compare("ShiNumberSeven") == 0){
		return new ShiNumberSevenFunction();

	}else if(name.compare("ShiBeale") == 0){
		return new ShiBealeFunction();

	}else if(name.compare("ShiStyblinskiTang") == 0){
		return new ShiStyblinskiTangFunction();

	}else if(name.compare("ShiCurrinExp") == 0){
		return new ShiCurrinExpFunction();

	}else if(name.compare("ShiLim") == 0){
		return new ShiLimFunction();

	}else if(name.compare("ShiGramacy") == 0){
		return new ShiGramacyFunction();

	}else if(name.compare("DongBohachevsky") == 0){
		return new DongBohachevskyFunction();

	}else if(name.compare("DongBooth") == 0){
		return new DongBoothFunction();

	}else if(name.compare("DongBranin") == 0){
		return new DongBraninFunction();

	}else if(name.compare("DongHimmelblau") == 0){
		return new DongHimmelblauFunction();

	}else if(name.compare("DongSixHumpCamelback") == 0){
		return new DongSixHumpCamelbackFunction();

	}else if(name.compare("XiongCurrinExp") == 0){
		return new XiongCurrinExpFunction();

	}else if(name.compare("MarchWillcoxRosenbrock1") == 0){
		return new MarchWillcoxRosenbrockFunction(1);

	}else if(name.compare("MarchWillcoxRosenbrock2") == 0){
		return new MarchWillcoxRosenbrockFunction(2);

	}else if(name.compare("MarchWillcoxRosenbrock3") == 0){
		return new MarchWillcoxRosenbrockFunction(3);

	}else if(name.compare("MarchWillcoxRosenbrock4") == 0){
		return new MarchWillcoxRosenbrockFunction(4);

	}else if(name.compare("MarchWillcoxRosenbrock5") == 0){
		return new MarchWillcoxRosenbrockFunction(5);

	}else if(name.compare("RajnarayanHartmannH3") == 0){
		return new RajnarayanHartmannH3Function();

	}else if(name.compare("ShiHartmannH3") == 0){
		return new ShiHartmannH3Function();

	}else if(name.compare("ShiDettePepelyshevExp") == 0){
		return new ShiDettePepelyshevExpFunction();

	}else if(name.compare("RajnarayanWoods") == 0){
		return new RajnarayanWoodsFunction();

	}else if(name.compare("ShiHartmannH4") == 0){
		return new ShiHartmannH4Function();

	}else if(name.compare("ShiPark") == 0){
		return new ShiParkFunction();

	}else if(name.compare("XiongParkFirst") == 0){
		return new XiongParkFirstFunction();

	}else if(name.compare("XiongParkSecond") == 0){
		return new XiongParkSecondFunction();

	}else if(name.compare("LiuStyblinskiTang") == 0){
		return new LiuStyblinskiTangFunction();

	}else if(name.compare("RajnarayanHartmannH6") == 0){
		return new RajnarayanHartmannH6Function();

	}else if(name.compare("ShiHartmannH6") == 0){
		return new ShiHartmannH6Function();

	}else if(name.compare("ShiRosenbrock") == 0){
		return new ShiRosenbrockFunction();

	}else if(name.compare("ParkHartmannH6") == 0){
		return new ParkHartmannH6Function();

	}else if(name.compare("ShiDettePepelyshev") == 0){
		return new ShiDettePepelyshevFunction();

	}else if(name.compare("LiuAckley10") == 0){
		return new LiuAckley10Function();

	}else if(name.compare("LiuEllipsoid") == 0){
		return new LiuEllipsoidFunction();

	}else if(name.compare("LiuDixonPrice") == 0){
		return new LiuDixonPriceFunction();

	}else if(name.compare("LiuAckley20") == 0){
		return new LiuAckley20Function();

	}else if(name.compare(0, 10, "ToalBranin") == 0){
		double a = stof(name.substr(10, 4));
		return new ToalBraninFunction(a);

	}else if(name.compare(0, 20, "NicoSongToalForretal") == 0){
		double a = stof(name.substr(20, 4));
		int dim = stof(name.substr(25));
		return new NicoSongToalForretalFunction(dim, a);

	}else if(name.compare(0, 16, "SongToalForretal") == 0){
		double a = stof(name.substr(16, 4));
		return new SongToalForretalFunction(a);

	}else if(name.compare(0, 12, "ToalPaciorek") == 0){
		double a = stof(name.substr(12, 4));
		return new ToalPaciorekFunction(a);

	}else if(name.compare(0, 14, "SongToalBranin") == 0){
		double a = stof(name.substr(14, 4));
		return new SongToalBraninFunction(a);

	}else if(name.compare(0, 14, "ToalHartmannH3") == 0){
		double a = stof(name.substr(14, 4));
		return new ToalHartmannH3Function(a);

	}else if(name.compare(0, 16, "SongToalColville") == 0){
		double a = stof(name.substr(16, 4));
		return new SongToalColvilleFunction(a);

	}else if(name.compare(0, 8, "ToalTrid") == 0){
		double a = stof(name.substr(8, 4));
		return new ToalTridFunction(a);

	}else if(name.compare(0, 5, "SOLAR") == 0){
		double fid = stof(name.substr(5, 4));
		int fileNum = 0;
		if(name.size() > 9){fileNum = stoi(name.substr(9));}
		return new SOLARFunction(fid, fileNum);

	}else if(name.compare(0, 13, "WangRastrigin") == 0){
		stringstream ss(name);
		string line;
		// skip name
		getline(ss, line, '-');
		// get error line
		getline(ss, line, '-');
		int dim = stoi(line.substr(1));
		getline(ss, line, '-');
		int error = stoi(line.substr(5));
		getline(ss, line, '-');
		double phi = stof(line.substr(3));
		return new WangRastriginFunction(dim, error, phi);
	
	}else if(name.compare(0, 12, "COCOfunction") == 0 || name.compare(0, 24, "DisturbanceBasedFunction") == 0){
		return processSpecialFunctionName(name, knowOptVals, knownFmin, knownFmax);
	}
	printf("Problem: Could not match function name!! Stopping here...\n");							
	return NULL;
}


BiFidelityFunction* processSpecialFunctionName(string name, bool knowOptVals, double knownFmin, double knownFmax){
	// Here should have a lot of information to deal with this in different ways, but for now just deal with the case I am interested in
	// Really, should first specify global noise, followed by parameters pertaining to global noise
	// Then local noise, followed by parameters pertaining to local noise
	// So let's do that
	stringstream ss(name);
	string line;
	// Get function number
	getline(ss, line, '-');
	int func;
	bool coco;
	COCOBiFunction* cocoFunction;
	disturbanceBasedBiFunction* disturbanceFunction;
	int dim;
	int seed;
	char disturbanceType;
	int disturbanceNum;
	double disturbanceHeight;
	double disturbanceRadius;
	int centres;
	int basicDisturbanceFrequency;
	double basicDisturbanceAmplitude;

	if(line.compare(0, 12, "COCOfunction") == 0){
		func = stoi(line.substr(12));
		coco = true;
	}else{
		func = stoi(line.substr(24));
		coco = false;
	}
	// Get dimension if this is a COCO function
	if(coco){
		getline(ss, line, '-');
		dim = stoi(line.substr(3));
	}
	// Get seed
	getline(ss, line, '-');
	seed = stoi(line.substr(4));
	// Get disturbance type
	getline(ss, line, '-');
	disturbanceType = line[4];
	disturbanceNum = stoi(line.substr(5,1));

	// Initialise and assign correct functions
	if(coco){
		cocoFunction = new COCOBiFunction(func, dim, seed);
		cocoFunction->disturbanceType_ = disturbanceType;
		cocoFunction->disturbanceNum_ = disturbanceNum;
	}else{
		disturbanceFunction = new disturbanceBasedBiFunction(func, seed);
		disturbanceFunction->disturbanceType_ = disturbanceType;
		disturbanceFunction->disturbanceNum_ = disturbanceNum;
	}

	
	if(disturbanceType == 'h'){
		// Get height
		getline(ss, line, '-');
		disturbanceHeight = stof(line.substr(6));
		// Get radius
		getline(ss, line, '-');
		disturbanceRadius = stof(line.substr(6));

		if(coco){
			cocoFunction->disturbanceHeight_ = disturbanceHeight;
			cocoFunction->disturbanceRadius_ = disturbanceRadius;
		}else{
			disturbanceFunction->disturbanceHeight_ = disturbanceHeight;
			disturbanceFunction->disturbanceRadius_ = disturbanceRadius;
		}

	}else if(disturbanceType == 's'){
		// Get number of sources and initialise points
		getline(ss, line, '-');
		centres = stoi(line.substr(7));
		// Get radius
		getline(ss, line, '-');
		disturbanceRadius = stof(line.substr(6));	
		if(coco){
			SampleGenerator* generator = new SampleGenerator(cocoFunction, seed, false);
			cocoFunction->disturbanceCentres_ = generator->randomLHS(centres);
			delete generator;
			cocoFunction->disturbanceRadius_ = disturbanceRadius;		
		}else{
			SampleGenerator* generator = new SampleGenerator(disturbanceFunction, seed, false);
			disturbanceFunction->disturbanceCentres_ = generator->randomLHS(centres);
			delete generator;
			disturbanceFunction->disturbanceRadius_ = disturbanceRadius;
		}
		
	}else{
		printf("Disturbance type %c not yet implemented, for now only have height (h) and source (s) based!\n", disturbanceType);
		return NULL;
	}
	// Get basic disturbance parameters, frequency and amplitude
	getline(ss, line, '-');
	basicDisturbanceFrequency = stoi(line.substr(4));
	getline(ss, line, '-');
	basicDisturbanceAmplitude = stof(line.substr(3));

	if(coco){
		cocoFunction->basicDisturbanceFrequency_ = basicDisturbanceFrequency;
		cocoFunction->basicDisturbanceAmplitude_ = basicDisturbanceAmplitude;

		if(knowOptVals){
			cocoFunction->fMax_ = knownFmax;
		}else{
			// Also need to know the maximum points, minimum is assigned already
			ARSsolver* auxSolver = new ARSsolver(cocoFunction, 10, 5000, false, seed, false);
			VectorXd best = auxSolver->optimise();
			cocoFunction->fMax_ = cocoFunction->evaluate(best);
			delete auxSolver;
		}
		return cocoFunction;
	
	}else{
		disturbanceFunction->basicDisturbanceFrequency_ = basicDisturbanceFrequency;
		disturbanceFunction->basicDisturbanceAmplitude_ = basicDisturbanceAmplitude;
		if(knowOptVals){
			disturbanceFunction->fMax_ = knownFmax;
			disturbanceFunction->fMin_ = knownFmin;
		}else{
			// Also need to know the maximum points, minimum is assigned already
			ARSsolver* auxSolver = new ARSsolver(disturbanceFunction, 10, 5000, false, seed, false);
			VectorXd best = auxSolver->optimise();
			disturbanceFunction->fMax_ = disturbanceFunction->evaluate(best);
			auxSolver->updateProblem(disturbanceFunction, true);
			best = auxSolver->optimise();
			disturbanceFunction->fMin_ = disturbanceFunction->evaluate(best);
			delete auxSolver;
		}
		return disturbanceFunction;
	}

}


vector<VectorXd> readPointsFromFile(string filename, int pointsNum, int dimension){
	ifstream samplePlanInput;
	vector<VectorXd> points;
	points.reserve(pointsNum);
	// First check if the file exists, if it does not return an empty vector
	samplePlanInput.open(filename);
	if(!samplePlanInput.is_open()){
		return {};
	}
	string line;
	for(int i = 0; i < pointsNum; i++){
		getline(samplePlanInput, line);
		VectorXd point(dimension);
		stringstream ss(line);
		string token;
		for(int j = 0; j < dimension; j++){
			getline(ss, token, ' ');
			point(j) = stof(token);
		}
		points.push_back(point);
	}    
    samplePlanInput.close();
    return points;
}


void writePointsToFile(string filename, vector<VectorXd> points, int pointsNum, int dimension){
	ofstream samplePlanOutput;
	// First check if the file exists, if it does not return an empty vector
	samplePlanOutput.open(filename);
	if(!samplePlanOutput.is_open()){printf("Error: Could not open output file %s to output experiment results! Stopping now...\n", filename.c_str()); exit(0);}
	for(int i = 0; i < pointsNum; i++){
		for(int j = 0; j < dimension; j++){
			samplePlanOutput << points[i](j) << " ";
		}
		samplePlanOutput << "\n";
	}
	samplePlanOutput.close();
}



void generateAndSaveSample(Function* function, int sampleSize, int subsetSampleSize, int seed, bool printInfo){
	clock_t tstart;
	tstart = clock();
	SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
	// Note will not really use this, just needed when initialising a model
	AuxSolver* auxSolver = new ARSsolver(function);
	vector<VectorXd> sampledPoints;
	vector<VectorXd> subsetSampledPoints;
	string samplePlanFilename;
	if(sampleSize < 2){
		printf("Please, specify a sample size of size 2 or larger!\n");
		return;
	}

	if(printInfo){
		printf("Generate sample: ");
		sampleGenerator->prePrint_ = "Generate sample: ";
	}
	samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(sampleSize) + "-s" + to_string(seed) + ".txt";
	sampledPoints = readPointsFromFile(samplePlanFilename, sampleSize, function->d_);
	if(sampledPoints.empty()){
		sampledPoints = sampleGenerator->randomLHS(sampleSize);
		scalePoints(sampledPoints, function);
	}
	// Gonna add something here to run a portion of the optimisation, then save it, then repeat until no improvement is found
	// First need to know the current "score" i.e. min distance between a pair of points.
	tuple<double, int, int> info = sampleGenerator->morrisMitchellCriterion(sampledPoints);
	double distance = get<0>(info);
	double nextDistance = distance;
	// No a do while loop
	while(true){
		sampleGenerator->morrisMitchellLocalOptimal(sampledPoints, 100);
		info = sampleGenerator->morrisMitchellCriterion(sampledPoints);
		nextDistance = get<0>(info);
		writePointsToFile(samplePlanFilename, sampledPoints, sampleSize, function->d_);
		if(abs(distance - nextDistance) > TOL){
			if(printInfo){printf("\nPausing sample optimisation, improved from %.4f to %.4f. Elapsed time since start: %.2f seconds. Saving to file and continuing...\n", distance, nextDistance, (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);}
			// Improved, want to rewrite file to save the current points, then continue
			distance = nextDistance;
		}else{
			// Otherwise optimisation is complete, get out of here!
			if(printInfo){printf("\nSample plan is optimised, final min point distance is %.4f.\n", nextDistance);}
			break;
		}
	}
	if(subsetSampleSize > sampleSize){printf("Can't find subset sample of size larger than actual sample!\n");}
	else if(subsetSampleSize > 0){
		if(printInfo){
			printf("Generate subset sample: ");
			sampleGenerator->prePrint_ = "Generate subset sample: ";
		}
		// Choose a subset
		samplePlanFilename = "../data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(sampleSize) + "-nH" + to_string(subsetSampleSize) + "-s" + to_string(seed) + ".txt";
		subsetSampledPoints = readPointsFromFile(samplePlanFilename, subsetSampleSize, function->d_);
		// Again, if empty need to find them and store them
		if(subsetSampledPoints.empty()){
			subsetSampledPoints = sampleGenerator->morrisMitchellSubset(sampledPoints, subsetSampleSize);
			writePointsToFile(samplePlanFilename, subsetSampledPoints, subsetSampleSize, function->d_);
			if(printInfo){printf("\n");}
		}else{
			printf("Read in.\n");
		}

	}
	
	if(printInfo){printf("Elapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);}
	
	delete sampleGenerator;
	delete auxSolver;
}




pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(Function* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo){
	clock_t tstart;
	tstart = clock();
	string samplePlanFilename;
	SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
	// Note will not really use this, just needed when initialising a model
	AuxSolver* auxSolver = new ARSsolver(function);

	vector<VectorXd> sampledPoints = {};
	vector<VectorXd> sampledPointsLow = {};
	pair<vector<VectorXd>, vector<VectorXd> > points;
	
	// First generate the sample, read it if the file exists
	if(lowFiBudget > 0){
		// First try to read both sets
		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPointsLow = readPointsFromFile(samplePlanFilename, lowFiBudget, function->d_);
		samplePlanFilename = "../data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		
		// If there is no low fidelity sample, need to choose both samples
		if(sampledPointsLow.empty()){
			if(printInfo){
				printf("Generating low fi data sample, choosing random LHS sample, NOT Morris-Mitchell locally optimal.\n");
				printf("Generate low fi data sample: ");
				sampleGenerator->prePrint_ = "Generate low fi data sample: ";
			}
			sampledPointsLow = sampleGenerator->randomLHS(lowFiBudget);
			if(printInfo){printf("\n");}
			scalePoints(sampledPointsLow, function);
			if(printInfo){
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
				printf("Generate high fi data sample: ");
				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
			}
			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
			if(printInfo){printf("\n");}
		}else if(sampledPoints.empty()){
			if(printInfo){
				printf("Read in sample locations for low fidelity sample only.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
				printf("Generate high fi data sample: ");
				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
			}
			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
			if(printInfo){printf("\n");}

		}else{
			if(printInfo){
				printf("Read in sample locations for low and high fidelity samples.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
			}
		}

		// Unscale both samples and return!
		unscalePoints(sampledPointsLow, function);
		unscalePoints(sampledPoints, function);
		
	
	}else{
		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		if(sampledPoints.empty()){
			if(printInfo){
				printf("Generating high fi data sample, choosing random LHS sample, NOT Morris-Mitchell locally optimal.\n");
				printf("Generate high fi data sample: ");
				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
			}
			sampledPoints = sampleGenerator->randomLHS(highFiBudget);
			if(printInfo){
				printf("Read in sample locations for low and high fidelity samples.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
			}
		}else{
			unscalePoints(sampledPoints, function);
			if(printInfo){
				printf("Read in sample locations for high fidelity samples.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
			}
		}
	}

	delete sampleGenerator;
	delete auxSolver;

	points = make_pair(sampledPoints, sampledPointsLow);

	return points;
}









#endif