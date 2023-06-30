#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP

#include "functions.hpp"
	
/*
============================================================================

                                PARENT FUNCTIONS

============================================================================
*/

Function::Function(int d, double lowBound, double upBound):
	d_(d){
		vector<double> lowerBound(d, lowBound);
		lowerBound_ = lowerBound;
		vector<double> upperBound(d, upBound);
		upperBound_ = upperBound;
}

Function::Function(int d, vector<double> lowerBound, vector<double> upperBound):
	d_(d),
	lowerBound_(lowerBound),
	upperBound_(upperBound){}

Function::~Function(){};

void Function::checkDimension(VectorXd &point){
	if(d_ != point.size()){
		printf("Called Function method on a point with wrong dimension! Stopping now...");
		exit(1);
	}
}

bool Function::pointWithinBounds(VectorXd &point){
	checkDimension(point);
	for(int i = 0; i < d_; i++){
		if(dimensionWithinBounds(point, i) != 0){
			return false;
		}
	}
	return true;
}

int Function::dimensionWithinBounds(VectorXd &point, int d){
	if(d >= point.size()){
		printf("Called dimension within bounds outside of point range! Stopping now...");
		exit(1);
	}
	if(point(d) >  upperBound_[d]){return 1;}
	else if(point(d) < lowerBound_[d]){return -1;}
	else{return 0;}
}

double Function::evaluate(VectorXd &point){
	return 0;
}

vector<double> Function::evaluateMany(vector<VectorXd> &points){
	vector<double> values;
	values.reserve(points.size());
	for(int i = 0; i < (int)points.size(); i++){
		values.push_back(evaluate(points[i]));
	}
	return values;
}

int Function::betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag){
	// if(abs(val1 - val2) < TOL){return 0;}
	// else if(val1 > val2){return 1;}
	// else{return -1;}
	return 0;
}

int Function::betterPoint(VectorXd &point1, VectorXd &point2, string flag){
	double val1 = evaluate(point1);
	double val2 = evaluate(point2);
	return betterPoint(point1, val1, point2, val2, flag);
}


BiFidelityFunction::BiFidelityFunction(int d, double lowBound, double upBound) :
	Function(d, lowBound, upBound){}

BiFidelityFunction::BiFidelityFunction(int d, vector<double> lowerBound, vector<double> upperBound) :
	Function(d, lowerBound, upperBound){}

BiFidelityFunction::~BiFidelityFunction(){}

double BiFidelityFunction::evaluate(VectorXd &point){
	return 0;
}

double BiFidelityFunction::evaluateLow(VectorXd &point){
	return 0;
}

vector<double> BiFidelityFunction::evaluateManyLow(vector<VectorXd> &points){
	vector<double> values;
	values.reserve(points.size());
	for(int i = 0; i < (int)points.size(); i++){
		values.push_back(evaluateLow(points[i]));
	}
	return values;
}





SOLARFunction::SOLARFunction(double fidelityLevel, int fileNum):
	// BiFidelityFunction(5, vector<double> {793, 2, 2, 0.01, 0.01}, vector<double> {995.0, 50.0, 30.0, 5.00, 5.00}),
	BiFidelityFunction(5, vector<double> {0, 0, 0, 0, 0}, vector<double> {1, 1, 1, 1, 1}),
	fidelityLevel_(fidelityLevel),
	fileNum_(fileNum){
		if(fidelityLevel_ <= 0 || fidelityLevel_ > 1){
			printf("Initialising SOLAR simulation function with fidelity outside of (0,1], this does not work! Stopping now...\n");
			exit(0);
		}
	}
SOLARFunction::~SOLARFunction(){}

VectorXd SOLARFunction::scalePoint(VectorXd point){
	point(0) = point(0)*(995 - 793) + 793;
	point(1) = point(1)*(50 - 2) + 2;
	point(2) = point(2)*(30 - 2) + 2;
	point(3) = point(3)*(5 - 0.01) + 0.01;
	point(4) = point(4)*(5 - 0.01) + 0.01;
	return point;
}

double SOLARFunction::callSimulation(VectorXd &inputPoint, double fidelity){
	VectorXd point = scalePoint(inputPoint);
	// Need to create a file, save the point to be evaluated, evaluate, then delete the file
	ofstream pointOutput;
	pointOutput.open("cpp_code\\solar\\point" + to_string(fileNum_) + ".txt");
	// File is now open, can write to it
	for(int i = 0; i < point.size(); i++){
		pointOutput << point(i) << " ";
	}
	pointOutput << "\n";
	pointOutput.close();
	string command = "cpp_code\\solar\\bin\\solar 10 cpp_code\\solar\\point" + to_string(fileNum_) + ".txt -fid=" + to_string(fidelity);
	// Got the code to run system and get output from
	// https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
	array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) {
        throw runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

	// Still need to extract the output!
	string filename = "cpp_code\\solar\\point" + to_string(fileNum_) + ".txt";
	remove(filename.c_str());

	// When the code executes, if a hidden constraint is violated, 1e+20 is returned. Replace this with upper bound on the function output
	double val = stof(result);

	if(val > 100000){
		val = 2500;
	}
	printf("%.2f\n", val);
	
	return val;
}

double SOLARFunction::evaluate(VectorXd &point){
	printf("High\n");
	return callSimulation(point, 1);
}

double SOLARFunction::evaluateLow(VectorXd &point){
	printf("Low\n");
	return callSimulation(point, fidelityLevel_);
}




NicoSongToalForretalFunction::NicoSongToalForretalFunction(int dim, double a):
	BiFidelityFunction(dim, 0, 1),
	a_(a){}
NicoSongToalForretalFunction::~NicoSongToalForretalFunction(){}
double NicoSongToalForretalFunction::evaluate(VectorXd &point){
	return pow(6 * point(0) - 2, 2) * sin(12 * point(0) - 4);
}
double NicoSongToalForretalFunction::evaluateLow(VectorXd &point){
	return (1 - a_ * a_ - 2 * a_) * evaluate(point) + 10 * (point(0) - 0.5) - 5;
}







ToalBraninFunction::ToalBraninFunction(double a):
	BiFidelityFunction(2, vector<double> {-5, 0}, vector<double> {10, 15}),
	a_(a){}
ToalBraninFunction::~ToalBraninFunction(){}
double ToalBraninFunction::evaluate(VectorXd &point){	
	return pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2) + 10*(1-1/(8*M_PI))*cos(point(0)) + 10;
}
double ToalBraninFunction::evaluateLow(VectorXd &point){
	return evaluate(point) - (a_ + 0.5)*pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2);
}

ToalPaciorekFunction::ToalPaciorekFunction(double a):
	BiFidelityFunction(2, vector<double> {0.3, 0.3}, vector<double> {1, 1}),
	a_(a){}
ToalPaciorekFunction::~ToalPaciorekFunction(){}
double ToalPaciorekFunction::evaluate(VectorXd &point){
	if(abs(point(0) * point(1)) < TOL){return 0.0;}
	return sin(1/(point(0) * point(1)));
}
double ToalPaciorekFunction::evaluateLow(VectorXd &point){
	if(abs(point(0) * point(1)) < TOL){return 0.0;}
	return evaluate(point) - 9 * a_ * a_ * cos(1/(point(0) * point(1)));
}

ToalHartmannH3Function::ToalHartmannH3Function(double a):
	BiFidelityFunction(3, vector<double> {0, 0, 0}, vector<double> {1, 1, 1}),
	a_(a){
		alpha_ = {1, 1.2, 3, 3.2};
		beta_ = {{3, 10, 30},
				{0.1, 10, 35},
				{3, 10, 30},
				{0.1, 10, 35}};
		p_ = {{0.3689, 0.1170, 0.2673},
			  {0.4699, 0.4387, 0.7470},
			  {0.1091, 0.8732, 0.5547},
			  {0.0381, 0.5743, 0.8828}};
}
ToalHartmannH3Function::~ToalHartmannH3Function(){}
double ToalHartmannH3Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 3; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}
double ToalHartmannH3Function::evaluateLow(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 3; j++){
			exponent += beta_[i][j] * pow(point(j) - (3.0/4.0) * p_[i][j] * (a_ + 1.0), 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}

ToalTridFunction::ToalTridFunction(double a):
	BiFidelityFunction(10, vector<double> {-100, -100, -100, -100, -100, -100, -100, -100, -100, -100}, 
							vector<double> {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}),
	a_(a){}
ToalTridFunction::~ToalTridFunction(){}
double ToalTridFunction::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 10; i++){
		total += pow(point(i) - 1, 2);
		if(i > 0){total -= point(i) * point(i-1);}
	}
	return total;
}
double ToalTridFunction::evaluateLow(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 10; i++){
		total += pow(point(i) - a_, 2);
		// This line is not the same as that presented by Toal, but it reproduces results!
		if(i > 0){total += (a_ - 0.65) * (i+1) * point(i) * point(i-1);}
		// This is the line as specified by Toal
		// if(i > 0){total -= (a_ - 0.65) * (i+1) * point(i) * point(i-1);}
	}
	return total;
}



SongToalForretalFunction::SongToalForretalFunction(double a):
	BiFidelityFunction(1, vector<double> {0}, 
							vector<double> {1}),
	a_(a){}
SongToalForretalFunction::~SongToalForretalFunction(){}
double SongToalForretalFunction::evaluate(VectorXd &point){
	return pow(6 * point(0) - 2, 2) * sin(12 * point(0) - 4);
}
double SongToalForretalFunction::evaluateLow(VectorXd &point){
	return (1 - a_ * a_ - 2 * a_) * evaluate(point) + 10 * (point(0) - 0.5) - 5;
}

SongToalBraninFunction::SongToalBraninFunction(double a):
	BiFidelityFunction(2, vector<double> {-5, 0}, vector<double> {10, 15}),
	a_(a){}
SongToalBraninFunction::~SongToalBraninFunction(){}
double SongToalBraninFunction::evaluate(VectorXd &point){
	return pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2) + 10*(1-1/(8*M_PI))*cos(point(0)) + 10;
}
double SongToalBraninFunction::evaluateLow(VectorXd &point){
	return evaluate(point) - (0.5 * a_ * a_ + a_ + 0.2) * pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2);
}

SongToalColvilleFunction::SongToalColvilleFunction(double a):
	BiFidelityFunction(4, vector<double> {-10, -10, -10, -10}, vector<double> {10, 10, 10, 10}),
	a_(a){}
SongToalColvilleFunction::~SongToalColvilleFunction(){}
double SongToalColvilleFunction::evaluate(VectorXd &point){
	return 100 * pow(pow(point(0), 2) - point(1), 2) + pow(point(0) - 1, 2) + pow(point(2) - 1, 2) + 90 * pow(pow(point(2), 2) - point(3), 2) + 10.1 * (pow(point(1) - 1, 2) + (pow(point(3) - 1, 2))) + 19.8 * (point(1) - 1) * (point(3) - 1);
}
double SongToalColvilleFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = a_ * a_ * point;
	return evaluate(tempPoint) - (a_ + 0.5) * (5 * point(0) * point(0) + 4 * point(1) * point(1) + 3 * point(2) * point(2) + point(3) * point(3));
}



MarchWillcoxRosenbrockFunction::MarchWillcoxRosenbrockFunction(int lowFiFunction):
	BiFidelityFunction(2, vector<double> {-5, -5}, vector<double> {10, 10}),
	lowFiFunction_(lowFiFunction){}
MarchWillcoxRosenbrockFunction::~MarchWillcoxRosenbrockFunction(){}
double MarchWillcoxRosenbrockFunction::evaluate(VectorXd &point){
	return pow(point(1) - pow(point(0), 2), 2) + pow(1 - point(0), 2);
}
double MarchWillcoxRosenbrockFunction::evaluateLow(VectorXd &point){
	if(lowFiFunction_ == 1){return 0.0;}
	else if(lowFiFunction_ == 2){return pow(point(0), 2) + pow(point(1), 2);}
	else if(lowFiFunction_ == 3){return pow(point(0), 4) + pow(point(1), 2);}
	else if(lowFiFunction_ == 4){return evaluate(point);}
	else if(lowFiFunction_ == 5){return -1*(pow(point(0), 2) + pow(point(1), 2));}
	else{
		printf("Calling MarchWillcoxRosenbrockFunction with wrong lowFiFunction option, should be int 1-5 but is %d! Stopping now...\n", lowFiFunction_);
		exit(0);
	}
	return 0.0;
}



RajnarayanHartmannH3Function::RajnarayanHartmannH3Function():
	BiFidelityFunction(3, vector<double> {0, 0, 0}, vector<double> {1, 1, 1}){
		alpha_ = {1, 1.2, 3, 3.2};
		beta_ = {{3, 10, 30},
				{0.1, 10, 35},
				{3, 10, 30},
				{0.1, 10, 35}};
		p_ = {{0.3689, 0.1170, 0.2673},
			  {0.4699, 0.4387, 0.7470},
			  {0.1091, 0.8732, 0.5547},
			  {0.0381, 0.5743, 0.8828}};
}
RajnarayanHartmannH3Function::~RajnarayanHartmannH3Function(){}
double RajnarayanHartmannH3Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 3; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}
double RajnarayanHartmannH3Function::evaluateLow(VectorXd &point){
	double extra = 0.0;
	for(int i = 0; i < 3; i++){
		extra += point(i) - 4;
	}
	return evaluate(point) + 0.5 * sin(10 * extra);
}

RajnarayanHartmannH6Function::RajnarayanHartmannH6Function():
	BiFidelityFunction(6, vector<double> {0, 0, 0, 0, 0, 0}, vector<double> {1, 1, 1, 1, 1, 1}){
		alpha_ = {1, 1.2, 3, 3.2};
		beta_ = {{10, 3, 17, 3.5, 1.7, 8},
				{0.05, 10, 17, 0.1, 8, 14},
				{3, 3.5, 1.7, 10, 17, 8},
				{17, 8, 0.05, 10, 0.1, 14}};
		p_ = {{0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
			  {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
			  {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
			  {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}};
}
RajnarayanHartmannH6Function::~RajnarayanHartmannH6Function(){}
double RajnarayanHartmannH6Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 6; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}
double RajnarayanHartmannH6Function::evaluateLow(VectorXd &point){
	double extra = 0.0;
	for(int i = 0; i < 6; i++){
		extra += point(i) - 4;
	}
	return evaluate(point) + 0.5 * sin(25 * extra);
}

RajnarayanWoodsFunction::RajnarayanWoodsFunction():
	BiFidelityFunction(4, vector<double> {-1, -1, -1, -1}, vector<double> {1, 1, 1, 1}){}
RajnarayanWoodsFunction::~RajnarayanWoodsFunction(){}
double RajnarayanWoodsFunction::evaluate(VectorXd &point){
	return 100 * pow(point(1) - point(0), 2) + 
			pow(1 - point(0), 2) + 
			90 * pow(point(3) - pow(point(2), 2), 2) + 
			pow(1 - point(2), 2) + 
			10.1 * (pow(1 - point(1), 2) + pow(1 - point(3), 2)) + 
			19.8 * (1 - point(1)) * (1 - point(3));
}
double RajnarayanWoodsFunction::evaluateLow(VectorXd &point){
	double extra = 0.0;
	for(int i = 0; i < 4; i++){
		extra += pow(point(i) - 0.5, 2);
	}
	return evaluate(point) + 10 * extra;
}



LiuEllipsoidFunction::LiuEllipsoidFunction():
	BiFidelityFunction(20, vector<double> {-30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30},
							vector<double> {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}){
		sH_ = {0.3, 0.4, 0.2, 0.6, 1, 0.9, 0.2, 0.8, 0.5, 0.7, 0.4, 0.3, 0.7, 1, 0.9, 0.6, 0.2, 0.8, 0.2, 0.5};
		sS_ = {1.8, 0.4, 2, 1.2, 1.4, 0.6, 1.6, 0.2, 0.8, 1, 1.3, 1.1, 2, 1.4, 0.5, 0.3, 1.6, 0.7, 0.3, 1.9};
}
LiuEllipsoidFunction::~LiuEllipsoidFunction(){}
double LiuEllipsoidFunction::evaluate(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += (i+1) * sH_[i] * pow(point(i) - sS_[i], 2.0);
	}
	return sum;
}
double LiuEllipsoidFunction::evaluateLow(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += (i+1) * point(i) * point(i);
	}
	return sum;
}

LiuDixonPriceFunction::LiuDixonPriceFunction():
	BiFidelityFunction(20, vector<double> {-30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30},
							vector<double> {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}){
		sS_ = {1.8, 0.5, 2, 1.2, 0.4, 0.2, 1.4, 0.3, 1.6, 0.6, 0.8, 1, 1.3, 1.9, 0.7, 1.6, 0.3, 1.1, 2, 1.4};
	}
LiuDixonPriceFunction::~LiuDixonPriceFunction(){}
double LiuDixonPriceFunction::evaluate(VectorXd &point){
	double sum = pow(point(0) - sS_[0] - 1, 2.0);
	for(int i = 1; i < d_; i++){
		sum += (i+1) * pow(2 * (point(i) - sS_[i]) * (point(i) - sS_[i]) - (point(i-1) - sS_[i-1]), 2);
	}
	return sum;
}
double LiuDixonPriceFunction::evaluateLow(VectorXd &point){
	double sum = pow(point(0) - 1, 2.0);
	for(int i = 1; i < d_; i++){
		sum += (i+1) * pow(2 * point(i) * point(i) - point(i-1), 2);
	}
	return sum;
}

LiuStyblinskiTangFunction::LiuStyblinskiTangFunction():
	BiFidelityFunction(5, vector<double> {-5, -5, -5, -5, -5},
							vector<double> {5, 5, 5, 5, 5}){
		sS_ = {0.28, 0.59, 0.47, 0.16, 0.32};
}
LiuStyblinskiTangFunction::~LiuStyblinskiTangFunction(){}
double LiuStyblinskiTangFunction::evaluate(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += pow(point(i) - sS_[i], 4) - 16 * pow(point(i) - sS_[i], 2) + 5 * (point(i) - sS_[i]);
	}
	return 0.5 * sum;
}
double LiuStyblinskiTangFunction::evaluateLow(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += pow(point(i), 4) - 16 * pow(point(i), 2) + 5 * point(i);
	}
	return 0.5 * sum;
}

LiuAckley10Function::LiuAckley10Function():
	BiFidelityFunction(10, vector<double> {-30, -30, -30, -30, -30, -30, -30, -30, -30, -30},
							vector<double> {30, 30, 30, 30, 30, 30, 30, 30, 30, 30}){
		sS_ = {1.3, 0.1, 1.4, 0.8, 1.7, 1, 1.5, 0.6, 2, 0.4};
		sF_ = 1.3;
}
LiuAckley10Function::~LiuAckley10Function(){}
double LiuAckley10Function::evaluate(VectorXd &point){
	double firstPow = 0.0;
	double secondPow = 0.0;
	for(int i = 0; i < d_; i++){
		firstPow += pow(point(i) - sS_[i], 2);
		secondPow += cos(2 * sF_ * M_PI * (point(i) - sS_[i]));
	}
	firstPow = -0.2 * sqrt(firstPow / d_);
	secondPow = secondPow / d_;
	return -20 * exp(firstPow) - exp(secondPow) + 20 + exp(1);
}
double LiuAckley10Function::evaluateLow(VectorXd &point){
	double firstPow = 0.0;
	double secondPow = 0.0;
	for(int i = 0; i < d_; i++){
		firstPow += point(i) * point(i);
		secondPow += cos(2 * M_PI * point(i));
	}
	firstPow = -0.2 * sqrt(firstPow / d_);
	secondPow = secondPow / d_;
	return -20 * exp(firstPow) - exp(secondPow) + 20 + exp(1);
}

LiuAckley20Function::LiuAckley20Function():
	BiFidelityFunction(20, vector<double> {-30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30},
							vector<double> {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}){
		sS_ = {1.2, 0.2, 1.4, 0.8, 1.8, 1, 1.6, 0.6, 2, 0.4, 1.3, 0.3, 1.5, 0.9, 1.9, 1.1, 1.7, 0.7, 2.1, 0.5};
		sF_ = 1.3;
}
LiuAckley20Function::~LiuAckley20Function(){}
double LiuAckley20Function::evaluate(VectorXd &point){
	double firstPow = 0.0;
	double secondPow = 0.0;
	for(int i = 0; i < d_; i++){
		firstPow += pow(point(i) - sS_[i], 2);
		secondPow += cos(2 * sF_ * M_PI * (point(i) - sS_[i]));
	}
	firstPow = -0.2 * sqrt(firstPow / d_);
	secondPow = secondPow / d_;
	return -20 * exp(firstPow) - exp(secondPow) + 20 + exp(1);
}
double LiuAckley20Function::evaluateLow(VectorXd &point){
	double firstPow = 0.0;
	double secondPow = 0.0;
	for(int i = 0; i < d_; i++){
		firstPow += point(i) * point(i);
		secondPow += cos(2 * M_PI * point(i));
	}
	firstPow = -0.2 * sqrt(firstPow / d_);
	secondPow = secondPow / d_;
	return -20 * exp(firstPow) - exp(secondPow) + 20 + exp(1);
}



WangRastriginFunction::WangRastriginFunction(int dimension, int error, double phi):
	BiFidelityFunction(dimension, -1, 1),
	error_(error),
	phi_(phi){
		if(error < 1 || error > 4){printf("Error when calling Wang Rastrigin function, undefined error! Should be between 1 and 4, got %d. Exiting now...\n", error); exit(0);}
		if(phi < 0 && phi > 10000 && abs(phi) > TOL && abs(10000 - phi) > TOL){printf("Error when calling Wang Rastrigin function, undefined phi! Should be between 0 and 10000, got %.2f. Exiting now...\n", phi); exit(0);}

		if(error == 1){
			thetaVal_ = 1 - 0.0001 * phi_;

		}else if(error == 2){
			thetaVal_ = exp(-0.00025 * phi);

		}else if(error == 3){
			if(phi < 1000){thetaVal_ = 1 - 0.0002 * phi;}
			else if(phi < 2000){thetaVal_ = 0.8;}
			else if(phi < 3000){thetaVal_ = 1.2 - 0.0002 * phi;}
			else if(phi < 4000){thetaVal_ = 0.6;}
			else if(phi < 5000){thetaVal_ = 1.4 - 0.0002 * phi;}
			else if(phi < 6000){thetaVal_ = 0.4;}
			else if(phi < 7000){thetaVal_ = 1.6 - 0.0002 * phi;}
			else if(phi < 8000){thetaVal_ = 0.2;}
			else if(phi < 9000){thetaVal_ = 1.8 - 0.0002 * phi;}
			else if(phi <= 10000){thetaVal_ = 0.0;}
			else{thetaVal_ = 0.0;}

		}else if(error == 4){
			thetaVal_ = 1 - 0.0001 * phi_;

		}

		aVal_ = thetaVal_;
		wVal_ = 10 * M_PI * thetaVal_;
		bVal_ = 0.5 * M_PI * thetaVal_;
	}
WangRastriginFunction::~WangRastriginFunction(){}
double WangRastriginFunction::evaluate(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += point(i) * point(i) + 1 - cos(10 * M_PI * point(i));
	}
	return sum;
}
double WangRastriginFunction::evaluateLow(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		if(error_ < 4){sum += aVal_ * cos(wVal_ * point(i) + bVal_ + M_PI);}
		// NOTE THIS SHOULD BE abs(x - xOpt) but for the given problem xOpt = 0
		else{sum += (1.0 - abs(point(i))) * aVal_ * cos(wVal_ * point(i) + bVal_ + M_PI);}
	}
	return evaluate(point) + sum;
}



LiuPedagogicalFunction::LiuPedagogicalFunction():
	BiFidelityFunction(1, vector<double> {0}, vector<double> {1}){}
LiuPedagogicalFunction::~LiuPedagogicalFunction(){}
double LiuPedagogicalFunction::evaluate(VectorXd &point){
	return 5 * pow(point(0), 2) * sin(12 * point(0));
}
double LiuPedagogicalFunction::evaluateLow(VectorXd &point){
	return evaluate(point) + (pow(point(0), 3) - 0.5) * sin(3 * point(0) - 0.5) + 4 * cos(2 * point(0));
}

LiuBraninFunction::LiuBraninFunction():
	BiFidelityFunction(2, vector<double> {-5, 0}, vector<double> {10, 15}){}
LiuBraninFunction::~LiuBraninFunction(){}
double LiuBraninFunction::evaluate(VectorXd &point){
	return pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2) + 10*(1-1/(8*M_PI))*cos(point(0)) + 10;
}
double LiuBraninFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) -= 2;
	tempPoint(1) -= 2;
	return 10 * sqrt(evaluate(tempPoint)) + 2 * (point(0) - 2.5) - 3 * (3 * point(1) - 7) - 1;
}



WuHimmelblauFunction::WuHimmelblauFunction():
	BiFidelityFunction(2, vector<double> {-3, -3}, vector<double> {3, 3}){}
WuHimmelblauFunction::~WuHimmelblauFunction(){}
double WuHimmelblauFunction::evaluate(VectorXd &point){
	return pow(pow(point(0), 2) + point(1) - 11, 2.0) + pow(point(0) + pow(point(1), 2.0) - 7, 2.0);
}
double WuHimmelblauFunction::evaluateLow(VectorXd &point){
	return pow(pow(0.5 * point(0), 2) + 0.8 * point(1) - 11, 2.0) + pow(point(0) + pow(point(1), 2.0) - 7, 2.0) + pow(point(1), 3) - pow(point(0) + 1, 2);
}



ShiGramacyLeeFunction::ShiGramacyLeeFunction():
	BiFidelityFunction(1, vector<double> {0}, vector<double> {1}){}
ShiGramacyLeeFunction::~ShiGramacyLeeFunction(){}
double ShiGramacyLeeFunction::evaluate(VectorXd &point){
	return sin(10 * M_PI * point(0)) / (2.0 * point(0)) + pow(point(0) - 1, 4.0);
}
double ShiGramacyLeeFunction::evaluateLow(VectorXd &point){
	return sin(10 * M_PI * point(0)) / point(0) + 2.0 * pow(point(0) - 1, 4.0);
}

ShiCurrinSinFunction::ShiCurrinSinFunction():
	BiFidelityFunction(1, vector<double> {0}, vector<double> {1}){}
ShiCurrinSinFunction::~ShiCurrinSinFunction(){}
double ShiCurrinSinFunction::evaluate(VectorXd &point){
	return sin(2 * M_PI * (point(0) - 0.1)) + pow(point(0), 2.0);
}
double ShiCurrinSinFunction::evaluateLow(VectorXd &point){
	return sin(2 * M_PI * (point(0) - 0.1));
}

ShiHolsclawFunction::ShiHolsclawFunction():
	BiFidelityFunction(1, vector<double> {0}, vector<double> {10}){}
ShiHolsclawFunction::~ShiHolsclawFunction(){}
double ShiHolsclawFunction::evaluate(VectorXd &point){
	return point(0) * sin(point(0)) / 10.0;
}
double ShiHolsclawFunction::evaluateLow(VectorXd &point){
	return (point(0) * sin(point(0)) + point(0)) / 10.0;
}

ShiSantnerFunction::ShiSantnerFunction():
	BiFidelityFunction(1, vector<double> {0}, vector<double> {1}){}
ShiSantnerFunction::~ShiSantnerFunction(){}
double ShiSantnerFunction::evaluate(VectorXd &point){
	return exp(-1.4 * point(0)) * cos(3.5 * M_PI * point(0));
}
double ShiSantnerFunction::evaluateLow(VectorXd &point){
	return exp(-1.4 * point(0)) * cos(3.5 * M_PI * point(0)) + 0.75 * pow(point(0), 2.0);
}

ShiBraninFunction::ShiBraninFunction():
	BiFidelityFunction(2, vector<double> {0, 0}, vector<double> {5, 5}){}
ShiBraninFunction::~ShiBraninFunction(){}
double ShiBraninFunction::evaluate(VectorXd &point){
	return pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2) + 10*(1-1/(8*M_PI))*cos(point(0)) + 10;
}
double ShiBraninFunction::evaluateLow(VectorXd &point){
	return (1.0 - 0.124 * M_PI) * cos(point(0));
}

ShiNumberSixFunction::ShiNumberSixFunction():
	BiFidelityFunction(2, vector<double> {-1, -1}, vector<double> {1, 1}){}
ShiNumberSixFunction::~ShiNumberSixFunction(){}
double ShiNumberSixFunction::evaluate(VectorXd &point){
	return 101 * pow(point(0), 2) + 101 * pow(pow(point(0), 2) + pow(point(1), 2), 2);
}
double ShiNumberSixFunction::evaluateLow(VectorXd &point){
	return pow(point(0), 2) + 100 * pow(pow(point(0), 2) + pow(point(1), 2), 4);
}

ShiNumberSevenFunction::ShiNumberSevenFunction():
	BiFidelityFunction(2, vector<double> {-5, -5}, vector<double> {10, 10}){}
ShiNumberSevenFunction::~ShiNumberSevenFunction(){}
double ShiNumberSevenFunction::evaluate(VectorXd &point){
	return pow(1.0 - 0.2 * point(1) + 0.05 * sin(4 * M_PI * point(1) - point(0)), 2.0) + pow(point(1) - 0.5 * sin(2 * M_PI * point(0)), 2.0);
}
double ShiNumberSevenFunction::evaluateLow(VectorXd &point){
	return pow(1.0 - 0.2 * point(1) + 0.05 * sin(4 * M_PI * point(1) - point(0)), 2.0) + 4.0 * pow(point(1) - 0.5 * sin(2 * M_PI * point(0)), 2.0);
}

ShiBealeFunction::ShiBealeFunction():
	BiFidelityFunction(2, vector<double> {-4.5, -4.5}, vector<double> {4.5, 4.5}){}
ShiBealeFunction::~ShiBealeFunction(){}
double ShiBealeFunction::evaluate(VectorXd &point){
	return pow(1.5 - point(0) + point(0) * point(1), 2.0) + pow(2.25 - point(0) + point(0) * pow(point(1), 2.0), 2.0) + pow(2.625 - point(0) + point(0) * pow(point(1), 3), 2.0);
}
double ShiBealeFunction::evaluateLow(VectorXd &point){
	return pow(1.5 - point(0) + point(0) * point(1), 2.0) + point(0) + point(1);
}

ShiStyblinskiTangFunction::ShiStyblinskiTangFunction():
	BiFidelityFunction(2, vector<double> {-4, -4}, vector<double> {4, 4}){}
ShiStyblinskiTangFunction::~ShiStyblinskiTangFunction(){}
double ShiStyblinskiTangFunction::evaluate(VectorXd &point){
	return pow(point(0), 4) - 16 * pow(point(0), 2) + 5 * point(0) + pow(point(1), 4) - 16 * pow(point(1), 2) + 5 * point(1);
}
double ShiStyblinskiTangFunction::evaluateLow(VectorXd &point){
	return pow(point(0), 4) - 16 * pow(point(0), 2) + pow(point(1), 4) - 16 * pow(point(1), 2);
}

ShiCurrinExpFunction::ShiCurrinExpFunction():
	BiFidelityFunction(2, vector<double> {0, 0}, vector<double> {0.5, 0.5}){}
ShiCurrinExpFunction::~ShiCurrinExpFunction(){}
double ShiCurrinExpFunction::evaluate(VectorXd &point){
	return (1.0 - exp(-1.0 / (2 * point(1)))) * (2300 * pow(point(0), 3) + 1900 * pow(point(0), 2) + 2092 * point(0) + 60) / (100 * pow(point(0), 3) + 500 * pow(point(0), 2) + 4 * point(0) + 20);
}
double ShiCurrinExpFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) = point(0) + 0.05;
	tempPoint(1) = point(1) + 0.05;
	double result = -2.0 * evaluate(tempPoint) / 5.0;

	tempPoint(0) = point(0) + 0.05;
	tempPoint(1) = max(0.0, point(1) - 0.05);
	result += 0.25 * evaluate(tempPoint);

	tempPoint(0) = point(0) - 0.05;
	tempPoint(1) = point(1) + 0.05;
	result += 0.25 * evaluate(tempPoint);

	tempPoint(0) = point(0) - 0.05;
	tempPoint(1) = max(0.0, point(1) - 0.05);
	result += 0.25 * evaluate(tempPoint);

	return result;
}

ShiLimFunction::ShiLimFunction():
	BiFidelityFunction(2, vector<double> {0, 0}, vector<double> {10, 10}){}
ShiLimFunction::~ShiLimFunction(){}
double ShiLimFunction::evaluate(VectorXd &point){
	return ((30.0 + 5.0 * point(0) * sin(5.0 * point(0))) * (4.0 + exp(-5 * point(1))) - 100.0) / 6.0;
}
double ShiLimFunction::evaluateLow(VectorXd &point){
	return ((30.0 + 5.0 * point(0) * sin(5.0 * point(0))) * (4.0 + 2.0 * exp(-5 * point(1)) / 5.0) - 100.0) / 6.0;
}

ShiGramacyFunction::ShiGramacyFunction():
	BiFidelityFunction(2, vector<double> {-2, -2}, vector<double> {2, 2}){}
ShiGramacyFunction::~ShiGramacyFunction(){}
double ShiGramacyFunction::evaluate(VectorXd &point){
	return point(0) * exp(-pow(point(0), 2) - pow(point(1), 2));
}
double ShiGramacyFunction::evaluateLow(VectorXd &point){
	return point(0) * exp(-pow(point(0), 2) - pow(point(1), 2)) + point(0) / 10.0;
}

ShiChengSanduFunction::ShiChengSanduFunction():
	BiFidelityFunction(2, vector<double> {-1, -1}, vector<double> {1, 1}){}
ShiChengSanduFunction::~ShiChengSanduFunction(){}
double ShiChengSanduFunction::evaluate(VectorXd &point){
	return cos(point(0) + point(1)) * exp(point(0) * point(1));
}
double ShiChengSanduFunction::evaluateLow(VectorXd &point){
	return cos(point(0) + point(1)) * exp(point(0) * point(1)) + cos(pow(point(0), 2) + pow(point(1), 2));
}

ShiHartmannH3Function::ShiHartmannH3Function():
	BiFidelityFunction(3, vector<double> {0, 0, 0}, vector<double> {1, 1, 1}){
		alpha_ = {1, 1.2, 3, 3.2};
		beta_ = {{3, 10, 30},
				{0.1, 10, 35},
				{3, 10, 30},
				{0.1, 10, 35}};
		p_ = {{0.3689, 0.1170, 0.2673},
			  {0.4699, 0.4387, 0.7470},
			  {0.1091, 0.8732, 0.5547},
			  {0.0381, 0.5743, 0.8828}};
}
ShiHartmannH3Function::~ShiHartmannH3Function(){}
double ShiHartmannH3Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 3; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}
double ShiHartmannH3Function::evaluateLow(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 3; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += exp(-1 * exponent);
	}
	return -1 * total;
}

ShiDettePepelyshevExpFunction::ShiDettePepelyshevExpFunction():
	BiFidelityFunction(3, vector<double> {0, 0, 0}, vector<double> {1, 1, 1}){}
ShiDettePepelyshevExpFunction::~ShiDettePepelyshevExpFunction(){}
double ShiDettePepelyshevExpFunction::evaluate(VectorXd &point){
	return 100 * (exp(-2.0 / pow(point(0), 1.75)) + exp(-2.0 / pow(point(1), 1.5)) + exp(-2.0 / pow(point(2), 1.25)));
}
double ShiDettePepelyshevExpFunction::evaluateLow(VectorXd &point){
	return 100 * (exp(-2.0 / pow(point(0), 1.75)) + exp(-2.0 / pow(point(1), 1.5)) + 0.2 * exp(-2.0 / pow(point(2), 1.25)));
}

ShiHartmannH4Function::ShiHartmannH4Function():
	BiFidelityFunction(4, vector<double> {0, 0, 0, 0}, vector<double> {1, 1, 1, 1}){
		alpha_ = {1, 1.2, 3, 3.2};
		beta_ = {{10, 3, 17, 3.5, 1.7, 8},
				{0.05, 10, 17, 0.1, 8, 14},
				{3, 3.5, 1.7, 10, 17, 8},
				{17, 8, 0.05, 10, 0.1, 14}};
		p_ = {{0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
			  {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
			  {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
			  {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}};
}
ShiHartmannH4Function::~ShiHartmannH4Function(){}
double ShiHartmannH4Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 4; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}
double ShiHartmannH4Function::evaluateLow(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 4; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += exp(-1 * exponent);
	}
	return -1 * total;
}

ShiParkFunction::ShiParkFunction():
	BiFidelityFunction(4, vector<double> {0, 0, 0, 0}, vector<double> {1, 1, 1, 1}){}
ShiParkFunction::~ShiParkFunction(){}
double ShiParkFunction::evaluate(VectorXd &point){
	return point(0) * (sqrt(1 + (point(0) + pow(point(2), 3.0)) * point(3) / pow(point(0), 2.0)) - 1.0)/ 2.0 + (point(0) + 3 * point(3)) * exp(1 + sin(point(2)));
}
double ShiParkFunction::evaluateLow(VectorXd &point){
	return 0.79 * (1.0 + sin(point(0)) / 10.0) * evaluate(point) - 2 * point(0) + pow(point(1), 2) + pow(point(2), 2) + 0.5;
}

ShiHartmannH6Function::ShiHartmannH6Function():
	BiFidelityFunction(6, vector<double> {0, 0, 0, 0, 0, 0}, vector<double> {1, 1, 1, 1, 1, 1}){
		alpha_ = {1, 1.2, 3, 3.2};
		beta_ = {{10, 3, 17, 3.5, 1.7, 8},
				{0.05, 10, 17, 0.1, 8, 14},
				{3, 3.5, 1.7, 10, 17, 8},
				{17, 8, 0.05, 10, 0.1, 14}};
		p_ = {{0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
			  {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
			  {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
			  {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}};
}
ShiHartmannH6Function::~ShiHartmannH6Function(){}
double ShiHartmannH6Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 6; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * total;
}
double ShiHartmannH6Function::evaluateLow(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 6; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += exp(-1 * exponent);
	}
	return -1 * total;
}

ShiRosenbrockFunction::ShiRosenbrockFunction():
	BiFidelityFunction(6, vector<double> {-5, -5, -5, -5, -5, -5}, vector<double> {10, 10, 10, 10, 10, 10}){}
ShiRosenbrockFunction::~ShiRosenbrockFunction(){}
double ShiRosenbrockFunction::evaluate(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		sum += 100 * pow(point(i + 1) - pow(point(i), 2), 2) + pow(point(i) - 1, 2);
	}
	return sum;
}
double ShiRosenbrockFunction::evaluateLow(VectorXd &point){
	double sum = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		sum += 100 * pow(point(i + 1) - point(i), 2) + 4 * pow(point(i) - 1, 2);
	}
	return sum;
}

ShiDettePepelyshevFunction::ShiDettePepelyshevFunction():
	BiFidelityFunction(8, vector<double> {0, 0, 0, 0, 0, 0, 0, 0}, vector<double> {1, 1, 1, 1, 1, 1, 1, 1}){}
ShiDettePepelyshevFunction::~ShiDettePepelyshevFunction(){}
double ShiDettePepelyshevFunction::evaluate(VectorXd &point){
	double sum = 4 * pow(point(0) - 2 + 8 * point(1) - 8 * pow(point(1), 2), 2) + pow(3 - 4 * point(1), 2) + 16 * sqrt(point(2) + 1) * pow(2 * point(2) - 1, 2);
	for(int i = 4; i <= 8; i++){
		double subSum = 0.0;
		for(int j = 3; j <= i; j++){
			subSum += point(j-1);
		}
		sum += i * log(1 + subSum);
	}
	return sum;
}
double ShiDettePepelyshevFunction::evaluateLow(VectorXd &point){
	double sum = 4 * pow(point(0) - 2 + 8 * point(1) - 8 * pow(point(1), 2), 2) + pow(3 - 4 * point(1), 2) + 16 * sqrt(point(2) + 1) * pow(2 * point(2) - 1, 2);
	for(int i = 4; i <= 8; i++){
		double subSum = 0.0;
		for(int j = 3; j <= i; j++){
			subSum += point(j-1);
		}
		sum += log(1 + subSum);
	}
	return sum;
}



DongBohachevskyFunction::DongBohachevskyFunction():
	BiFidelityFunction(2, vector<double> {-100, -100}, vector<double> {100, 100}){}
DongBohachevskyFunction::~DongBohachevskyFunction(){}
double DongBohachevskyFunction::evaluate(VectorXd &point){
	return pow(point(0), 2) + 2 * pow(point(1), 2) - 0.3 * cos(3 * M_PI * point(0)) - 0.4 * cos(4 * M_PI * point(1)) + 0.7;
}
double DongBohachevskyFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) = 0.7 * point(0);
	return evaluate(tempPoint) + point(0) * point(1) - 12;
}

DongBoothFunction::DongBoothFunction():
	BiFidelityFunction(2, vector<double> {-10, -10}, vector<double> {10, 10}){}
DongBoothFunction::~DongBoothFunction(){}
double DongBoothFunction::evaluate(VectorXd &point){
	return pow(point(0) + 2 * point(1) - 7, 2) + pow(2 * point(0) + point(1) - 5, 2);
}
double DongBoothFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) = 0.4 * point(0);
	return evaluate(tempPoint) + 1.7 * point(0) * point(1) - point(0) + 2 * point(1);
}

DongBraninFunction::DongBraninFunction():
	BiFidelityFunction(2, vector<double> {-5, 0}, vector<double> {10, 15}){}
DongBraninFunction::~DongBraninFunction(){}
double DongBraninFunction::evaluate(VectorXd &point){
	return pow(point(1) - 5.1*pow(point(0),2)/(4*pow(M_PI,2)) + 5*point(0)/M_PI - 6, 2) + 10*(1-1/(8*M_PI))*cos(point(0)) + 10 - 22.5 * point(1);
}
double DongBraninFunction::evaluateLow(VectorXd &point){
	return pow(0.7 * point(1) - 5.1*pow(0.7 * point(0),2)/(4*pow(M_PI,2)) + 5*0.7*point(0)/M_PI - 6, 2) + 10*(1-1/(8*M_PI))*cos(0.7*point(0)) + 10 - 15.75 * point(1) + 20 * pow(0.9 + point(0), 2) - 50;
}

DongHimmelblauFunction::DongHimmelblauFunction():
	BiFidelityFunction(2, vector<double> {-3, -3}, vector<double> {3, 3}){}
DongHimmelblauFunction::~DongHimmelblauFunction(){}
double DongHimmelblauFunction::evaluate(VectorXd &point){
	return pow(pow(point(0),2) + point(1) - 11, 2) + pow(pow(point(1), 2) + point(0) - 7, 2);
}
double DongHimmelblauFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) = 0.5 * point(0);
	tempPoint(1) = 0.8 * point(1);
	return evaluate(tempPoint) + pow(point(1), 3) - pow(point(0) + 1, 2);
}

DongSixHumpCamelbackFunction::DongSixHumpCamelbackFunction():
	BiFidelityFunction(2, vector<double> {-2, -2}, vector<double> {2, 2}){}
DongSixHumpCamelbackFunction::~DongSixHumpCamelbackFunction(){}
double DongSixHumpCamelbackFunction::evaluate(VectorXd &point){
	return 4 * pow(point(0), 2) - 2.1 * pow(point(0), 4) + pow(point(0), 6)/3 + point(0) * point(1) - 4 * pow(point(1), 2) + 4 * pow(point(1), 4);
}
double DongSixHumpCamelbackFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) = 0.7 * point(0);
	tempPoint(1) = 0.7 * point(1);
	return evaluate(tempPoint) + point(0) * point(1) - 15;
}



XiongCurrinExpFunction::XiongCurrinExpFunction():
	BiFidelityFunction(2, vector<double> {0, 0}, vector<double> {1, 1}){}
XiongCurrinExpFunction::~XiongCurrinExpFunction(){}
double XiongCurrinExpFunction::evaluate(VectorXd &point){
	return (1.0 - exp(-1.0 / (2 * point(1)))) * (2300 * pow(point(0), 3) + 1900 * pow(point(0), 2) + 2092 * point(0) + 60) / (100 * pow(point(0), 3) + 500 * pow(point(0), 2) + 4 * point(0) + 20);
}
double XiongCurrinExpFunction::evaluateLow(VectorXd &point){
	VectorXd tempPoint = point;
	tempPoint(0) = point(0) + 0.05;
	tempPoint(1) = point(1) + 0.05;
	double result = evaluate(tempPoint);

	tempPoint(0) = point(0) + 0.05;
	tempPoint(1) = max(0.0, point(1) - 0.05);
	result += evaluate(tempPoint);

	tempPoint(0) = point(0) - 0.05;
	tempPoint(1) = point(1) + 0.05;
	result += evaluate(tempPoint);

	tempPoint(0) = point(0) - 0.05;
	tempPoint(1) = max(0.0, point(1) - 0.05);
	result += evaluate(tempPoint);

	return result / 4.0;
}

XiongParkFirstFunction::XiongParkFirstFunction():
	BiFidelityFunction(4, vector<double> {0, 0, 0, 0}, vector<double> {1, 1, 1, 1}){}
XiongParkFirstFunction::~XiongParkFirstFunction(){}
double XiongParkFirstFunction::evaluate(VectorXd &point){
	return point(0) * (sqrt(1 + (point(0) + pow(point(2), 3.0)) * point(3) / pow(point(0), 2.0)) - 1.0)/ 2.0 + (point(0) + 3 * point(3)) * exp(1 + sin(point(2)));
}
double XiongParkFirstFunction::evaluateLow(VectorXd &point){
	return (1.0 + sin(point(0)) / 10.0) * evaluate(point) - 2 * point(0) + pow(point(1), 2) + pow(point(2), 2) + 0.5;
}

XiongParkSecondFunction::XiongParkSecondFunction():
	BiFidelityFunction(4, vector<double> {0, 0, 0, 0}, vector<double> {1, 1, 1, 1}){}
XiongParkSecondFunction::~XiongParkSecondFunction(){}
double XiongParkSecondFunction::evaluate(VectorXd &point){
	return 2 * exp(point(0) + point(1)) / 3 - point(3) * sin(point(2)) + point(2);
}
double XiongParkSecondFunction::evaluateLow(VectorXd &point){
	return 1.2 * evaluate(point) - 1;
}



ParkHartmannH6Function::ParkHartmannH6Function():
	BiFidelityFunction(6, vector<double> {0, 0, 0, 0, 0, 0}, vector<double> {1, 1, 1, 1, 1, 1}){
		alpha_ = {1, 1.2, 3, 3.2};
		alphaDash_ = {0.5, 0.5, 2.0, 4.0};
		beta_ = {{10, 3, 17, 3.5, 1.7, 8},
				{0.05, 10, 17, 0.1, 8, 14},
				{3, 3.5, 1.7, 10, 17, 8},
				{17, 8, 0.05, 10, 0.1, 14}};
		p_ = {{0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
			  {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
			  {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
			  {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}};
}
ParkHartmannH6Function::~ParkHartmannH6Function(){}
double ParkHartmannH6Function::evaluate(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 6; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		total += alpha_[i] * exp(-1 * exponent);
	}
	return -1 * (2.58 + total) / 1.94;
}
double ParkHartmannH6Function::evaluateLow(VectorXd &point){
	double total = 0.0;
	for(int i = 0; i < 4; i++){
		double exponent = 0.0;
		for(int j = 0; j < 6; j++){
			exponent += beta_[i][j] * pow(point(j) - p_[i][j], 2);
		}
		exponent = -1 * exponent;
		total += alphaDash_[i] * pow(exp(-4.0 / 9.0) + exp(-4.0 / 9.0) * (exponent + 4.0) / 9.0, 9);
	}
	return -1 * (2.58 + total) / 1.94;
}







disturbanceBasedBiFunction::disturbanceBasedBiFunction(int function, int seed):
	BiFidelityFunction(0, 0, 0),
	function_(function),
	seed_(seed){

	// First initialise a bififunction to be used as the high fidelity source
	if(function_ == 1){
		// Classic Co-Kriging example
		highFiFunction_ = new LiuPedagogicalFunction();
	
	}else if(function_ == 2){
		// ShiGramacyLee
		highFiFunction_ = new ShiGramacyLeeFunction();
	
	}else if(function_ == 3){
		// ShiCurrinSin
		highFiFunction_ = new ShiCurrinSinFunction();
	
	}else if(function_ == 4){
		// ShiHolsclaw
		highFiFunction_ = new ShiHolsclawFunction();
	
	}else if(function_ == 5){
		// ShiSantner
		highFiFunction_ = new ShiSantnerFunction();
	
	}else if(function_ == 6){
		// Branin
		highFiFunction_ = new LiuBraninFunction();
	
	}else if(function_ == 7){
		// ShiNumberSix
		highFiFunction_ = new ShiNumberSixFunction();
	
	}else if(function_ == 8){
		// ShiNumberSeven
		highFiFunction_ = new ShiNumberSevenFunction();
	
	}else if(function_ == 9){
		// ShiBeale
		highFiFunction_ = new ShiBealeFunction();
	
	}else if(function_ == 10){
		// ShiStyblinskiTang
		highFiFunction_ = new ShiStyblinskiTangFunction();
	
	}else if(function_ == 11){
		// ShiCurrinExp
		highFiFunction_ = new ShiCurrinExpFunction();
	
	}else if(function_ == 12){
		// ShiLim
		highFiFunction_ = new ShiLimFunction();
	
	}else if(function_ == 13){
		// ShiGramacy
		highFiFunction_ = new ShiGramacyFunction();
	
	}else if(function_ == 14){
		// DongBohachevsky
		highFiFunction_ = new DongBohachevskyFunction();
	
	}else if(function_ == 15){
		// DongBooth
		highFiFunction_ = new DongBoothFunction();
	
	}else if(function_ == 16){
		// DongHimmelblau
		highFiFunction_ = new DongHimmelblauFunction();
	
	}else if(function_ == 17){
		// DongSixHumpCamelback
		highFiFunction_ = new DongSixHumpCamelbackFunction();
	
	}else if(function_ == 18){
		// XiongCurrinExp
		highFiFunction_ = new XiongCurrinExpFunction();
	
	}else if(function_ == 19){
		// MarchWillcoxRosenbrock
		highFiFunction_ = new MarchWillcoxRosenbrockFunction(1);
	
	}else if(function_ == 20){
		// Paciorek
		highFiFunction_ = new ToalPaciorekFunction(0);
	
	}else if(function_ == 21){
		// HartmannH3
		highFiFunction_ = new RajnarayanHartmannH3Function();
	
	}else if(function_ == 22){
		// ShiDettePepelyshevExp
		highFiFunction_ = new ShiDettePepelyshevExpFunction();
	
	}else if(function_ == 23){
		// Woods
		highFiFunction_ = new RajnarayanWoodsFunction();
	
	}else if(function_ == 24){
		// HartmannH4
		highFiFunction_ = new ShiHartmannH4Function();
	
	}else if(function_ == 25){
		// Park
		highFiFunction_ = new ShiParkFunction();
	
	}else if(function_ == 26){
		// XiongParkSecond
		highFiFunction_ = new XiongParkSecondFunction();
	
	}else if(function_ == 27){
		// Colville
		highFiFunction_ = new SongToalColvilleFunction(0);
	
	}else if(function_ == 28){
		// LiuStyblinskiTang
		highFiFunction_ = new LiuStyblinskiTangFunction();
	
	}else if(function_ == 29){
		// Rastrigin Function (note parameters are ignored for the purposes of this function)
		highFiFunction_ = new WangRastriginFunction(5, 1, 1000);
	
	}else if(function_ == 30){
		// HartmannH6
		highFiFunction_ = new ParkHartmannH6Function();
	
	}else if(function_ == 31){
		// Rosenbrock
		highFiFunction_ = new ShiRosenbrockFunction();
	
	}else if(function_ == 32){
		// ShiDettePepelyshev
		highFiFunction_ = new ShiDettePepelyshevFunction();
	
	}else if(function_ == 33){
		// LiuAckley10
		highFiFunction_ = new LiuAckley10Function();
	
	}else if(function_ == 34){
		// Trid
		highFiFunction_ = new ToalTridFunction(0);
	
	}else if(function_ == 35){
		// Rastrigin Function (note parameters are ignored for the purposes of this function)
		highFiFunction_ = new WangRastriginFunction(10, 1, 1000);
	
	}else if(function_ == 36){
		// Ellipsoid
		highFiFunction_ = new LiuEllipsoidFunction();
	
	}else if(function_ == 37){
		// LiuDixonPrice
		highFiFunction_ = new LiuDixonPriceFunction();
	
	}else if(function_ == 38){
		// LiuAckley20
		highFiFunction_ = new LiuAckley20Function();
	
	}else{
		printf("Defined a disturbance based function but with an incorrect integer! Please choose from 1 to 40. Exiting now...\n");
		exit(0);
	}
	
	// Save the dimension and bounds based on the chosen high fi function
	d_ = highFiFunction_->d_;
	lowerBound_ = highFiFunction_->lowerBound_;
	upperBound_ = highFiFunction_->upperBound_;

	// Work out maximum distance of two points
	maxDist_ = 0.0;
	VectorXd furthestPoint(d_);
	for(int i = 0; i < d_; i++){
		furthestPoint(i) = upperBound_[i] - lowerBound_[i];
	}
	maxDist_ = furthestPoint.norm();
}

disturbanceBasedBiFunction::~disturbanceBasedBiFunction(){}

double disturbanceBasedBiFunction::evaluate(VectorXd &point){
	return highFiFunction_->evaluate(point);
}

double disturbanceBasedBiFunction::evaluateLow(VectorXd &point){
	// Check noise came out alright
	double value = evaluate(point);
	bool addNoise = false;
	double dist = 0.0;
	double mult = 0.0;

	if(isnan(value)){
		printf("Encountered nan value from high fi at point (");
		for(int i = 0; i < (int)point.size() - 1; i++){
			printf("%.2f,", point(i));
			
		}
		printf("%.2f)! Stopping now...\n", point((int)point.size() - 1));
		exit(0);
	}
	// Work out value of sine and cosine sine noise, will work out a global position
	// and then add it if required
	VectorXd minDist(d_);
	for(int i = 0; i < d_; i++){minDist(i) = lowerBound_[i];}
	double trav = (point - minDist).norm() / maxDist_;
	// Calculate whether noise needs to be added, and its impact (mult value)
	if(disturbanceType_ == 'h' && disturbanceNum_ == 1){
		// Height based
		double noiseCentre = fMin_ + disturbanceHeight_ * (fMax_ - fMin_);
		double absRadius = disturbanceRadius_ * (fMax_ - fMin_);
		dist = abs(value - noiseCentre);
		if(dist < absRadius){
			addNoise = true;
			mult = 1.0 - dist / absRadius;
		}

	}else if(disturbanceType_ == 'h' && disturbanceNum_ == 2){
		// Height based
		double noiseCentre = fMin_ + disturbanceHeight_ * (fMax_ - fMin_);
		double absRadius = disturbanceRadius_ * (fMax_ - fMin_);
		dist = abs(value - noiseCentre);
		if(dist > absRadius){
			addNoise = true;
			mult = 1.0 - ((fMax_ - fMin_) - dist) / ((fMax_ - fMin_) - absRadius);
		}

	}else if(disturbanceType_ == 's' && disturbanceNum_ == 1){
		// Noise centres, work out closest point
		dist = DBL_MAX;
		for(int i = 0; i < (int)disturbanceCentres_.size(); i++){
			if((disturbanceCentres_[i] - point).norm() < dist){dist = (disturbanceCentres_[i] - point).norm();}
		}
		double absRadius = disturbanceRadius_ * maxDist_;
		if(dist > absRadius){
			addNoise = true;
			mult = 1.0 - (maxDist_ - dist) / (maxDist_ - absRadius);
		}
	}else if(disturbanceType_ == 's' && disturbanceNum_ == 2){
		// Disturbance centres, work out closest point
		dist = DBL_MAX;
		for(int i = 0; i < (int)disturbanceCentres_.size(); i++){
			if((disturbanceCentres_[i] - point).norm() < dist){dist = (disturbanceCentres_[i] - point).norm();}
		}
		double absRadius = disturbanceRadius_ * maxDist_;
		if(dist < absRadius){
			addNoise = true;
			mult = 1.0 - dist / absRadius;
		}
	
	}else{
		printf("Asking for Disturbance Based low fidelity evaluation with undefined pair (%c,%d) of disturbance values! Stopping now...\n", disturbanceType_, disturbanceNum_);
		exit(0);
	}

	// Add disturbance if it applies
	if(addNoise){
		value += mult * addBasicDisturbance(trav);
	}
	if(isnan(value)){
		printf("Mult %.5f dist %.5f\n", mult, addBasicDisturbance(trav));
		printf("fMax %.5f fMin %.5f absRadius %.5f\n", fMax_, fMin_, disturbanceRadius_ * (fMax_ - fMin_));
		printf("(Updated) Encountered nan value at point (");
		for(int i = 0; i < (int)point.size() - 1; i++){
			printf("%.2f,", point(i));
			
		}
		printf("%.2f)! Stopping now...\n", point((int)point.size() - 1));
		exit(0);
	}
	return value;

}

double disturbanceBasedBiFunction::addBasicDisturbance(double value){
	return basicDisturbanceAmplitude_ * (fMax_ - fMin_) * cos(basicDisturbanceFrequency_ * 2 * M_PI * value) * sin(basicDisturbanceFrequency_ * 2 * M_PI * value * value);
	
}









COCOBiFunction::COCOBiFunction(int function, int dimension, int seed):
	BiFidelityFunction(dimension, -5, 5),
	function_(function),
	seed_(seed){

	initialiseConstants();
	chooseOptVals();
}

COCOBiFunction::~COCOBiFunction(){}

void COCOBiFunction::initialiseConstants(){
	// Get random generator
	random_device rd;
    mt19937 gen(rd());
    if(seed_){
        mt19937 gen2(seed_);
        gen = gen2;
    }
    randGen_ = gen;

    // Calculate rotations, save matrix multiplications for certain functions
	rotationQ_ = randomRotationMatrix(d_);
	rotationR_ = randomRotationMatrix(d_);
	if(function_ == 6 || function_ == 13){leftMultiplication_ = rotationQ_ * matrixTalpha(10.0) * rotationR_;}
	if(function_ == 7){leftMultiplication_ = matrixTalpha(10.0) * rotationR_;}
	if(function_ == 15){leftMultiplication_ = rotationR_ * matrixTalpha(10.0) * rotationQ_;}
	if(function_ == 16){leftMultiplication_ = rotationR_ * matrixTalpha(0.01) * rotationQ_;}
	if(function_ == 17){leftMultiplication_ = matrixTalpha(10.0) * rotationQ_;}
	if(function_ == 18){leftMultiplication_ = matrixTalpha(1000.0) * rotationQ_;}
	if(function_ == 23 || function_ == 24){leftMultiplication_ = rotationQ_ * matrixTalpha(100.0) * rotationR_;}

	// Initialise a random ordering of alphas for relevant functions
	if(function_ == 21 || function_ == 22){
		int num;
		if(function_ == 21){num = 101;}
		else{num = 21;}
		alphas_.reserve(num);
		localOptima_.reserve(num);
		matrixPerm_.reserve(num);
		// Save all alpha values and dummy vectors of optima
		for(int i = 0; i < num; i++){
			if(i == 0){
				if(function_ == 21){alphas_.push_back(1000);}
				else{alphas_.push_back(1000000);}
			}
			else{alphas_.push_back(pow(1000.0, 2.0 * (i - 1.0) / (num - 2.0)));}

			VectorXd localOptimum(d_);
			localOptima_.push_back(localOptimum);
		}
		// Shuffle all alpha values except first
   		shuffle(alphas_.begin() + 1, alphas_.end(), randGen_);
   		// Generate random optima
   		double smallRange;
   		double bigRange;
   		if(function_ == 21){smallRange = 4; bigRange = 5;}
		else{smallRange = 3.92; bigRange = 4.9;}

   		uniform_real_distribution<double> randomDist(-1 * bigRange, bigRange);
   		for(int i = 0; i < num; i++){
   			for(int j = 0; j < d_; j++){
   				if(i == 0){localOptima_[i](j) = smallRange * randomDist(randGen_) / bigRange;}
   				else{localOptima_[i](j) = randomDist(randGen_);}
   			}
   		}
   		// Generate random shuffle of matrix columns
   		for(int i = 0; i < num; i++){
   			vector<int> order;
   			order.reserve(d_);
   			for(int j = 0; j < d_; j++){
   				order.push_back(j);
   			}
   			shuffle(order.begin(), order.end(), randGen_);
   			matrixPerm_.push_back(order);
   		}
	}
}

void COCOBiFunction::chooseOptVals(){
	uniform_real_distribution<double> randomDist(-4, 4);
	VectorXd xOpt(d_);
	// Generate location of optimum
	for(int i = 0; i < d_; i++){
		if(function_ == 5){
			xOpt(i) = randomDist(randGen_);
			if(xOpt(i) > 0){xOpt(i) = 5.0;}
			else{xOpt(i) = -5.0;}
		
		}else if(function_ == 8){
			// Pick from [-3,3] instead of [-4,4]
			xOpt(i) = 3.0 * randomDist(randGen_) / 4.0;
		}else if(function_ == 9 || function_ == 19){
			// xOpt not used, but saving it so it can be analysed
			xOpt(i) = 1.0 / (2.0 * max(1.0, sqrt(d_) / 8.0));
		}else if(function_ == 20){
			xOpt(i) = randomDist(randGen_);
			if(xOpt(i) > 0){xOpt(i) = 0.5 * 4.2096874633;}
			else{xOpt(i) = -0.5 * 4.2096874633;}
		
		}else if(function_ == 24){
			xOpt(i) = randomDist(randGen_);
			if(xOpt(i) > 0){xOpt(i) = 2.5 / 2.0;}
			else{xOpt(i) = -2.5 / 2.0;}

		}else{
			xOpt(i) = randomDist(randGen_);
		}
	}
	if(function_ == 9 || function_ == 19){xOpt_ = rotationR_.inverse() * xOpt;}
	else if(function_ == 21 || function_ == 22){xOpt_ = localOptima_[0];}
	else{xOpt_ = xOpt;}
	// Generate optimum value
	cauchy_distribution<double> randomDist2(0, 100);
	fOpt_ = min(100.0, max(-100.0, randomDist2(randGen_)));
	// Calculate maximum distance between any pair of points (distance from endpoints of sample space)
	maxDist_ = 0.0;
	VectorXd furthestPoint(d_);
	for(int i = 0; i < d_; i++){
		furthestPoint(i) = upperBound_[i] - lowerBound_[i];
	}
	maxDist_ = furthestPoint.norm();
}



MatrixXd COCOBiFunction::randomRotationMatrix(int d){
    normal_distribution<double> randomDist(0, 1);
    // Create matrix, fill it in with random numbers
    MatrixXd rotation(d,d);
    for(int i = 0; i < d; i++){
    	for(int j = 0; j < d; j++){
    		rotation(i,j) = randomDist(randGen_);
    	}
    }
    // Apply Gram-Schmidt process, copy pasted and applied from COCO code
    for(int i = 0; i < d; i++){
        for(int j = 0; j < i; j++){
            double prod = 0;
            for(int k = 0; k < d; k++){
                prod += rotation(k,i) * rotation(k,j);
            }
            for(int k = 0; k < d; k++){
                rotation(k,i) -= prod * rotation(k,j);
            }
        }
        double prod = 0;
        for(int k = 0; k < d; k++){
            prod += rotation(k,i) * rotation(k,i);
        }
        for (int k = 0; k < d; k++){
            rotation(k,i) /= sqrt(prod);
        }
    }
    return rotation;
}

MatrixXd COCOBiFunction::matrixTalpha(double alpha){
	MatrixXd alphaMatrix(d_, d_);

	for(int i = 0; i < d_; i++){
		for(int j = i + 1; j < d_; j++){
			alphaMatrix(i,j) = 0;
			alphaMatrix(j,i) = 0;
		}
		double power;
		if(d_ == 1){power = 0.0;}
		else{power = i / (2.0*(d_ - 1.0));}
		alphaMatrix(i,i) = pow(alpha, power);
	}
	return alphaMatrix;
}

MatrixXd COCOBiFunction::matrixC(int index){
	MatrixXd cMatrix(d_, d_);
	double alpha = alphas_[index]; 

	for(int i = 0; i < d_; i++){
		for(int j = i + 1; j < d_; j++){
			cMatrix(i,j) = 0;
			cMatrix(j,i) = 0;
		}
		double power;
		if(d_ == 1){power = 0;}
		else{power = matrixPerm_[index][i] / (2.0*(d_ - 1.0));}
		cMatrix(i,i) = pow(alpha, power) / pow(alpha, 0.25);
	}
	return cMatrix;
}

double COCOBiFunction::fTosc(double val){
	double a = 0.1;
	if (val > 0){
        val = log(val)/a;
        val = pow(exp(val + 0.49*(sin(val) + sin(0.79*val))), a);
    
    }else if (val < 0){
        val = log(-val)/a;
        val = -pow(exp(val + 0.49*(sin(0.55 * val) + sin(0.31*val))), a);
    }
    return val;
}

void COCOBiFunction::fTosc(VectorXd &point){
    for(int i = 0; i < d_; i++){
    	point(i) = fTosc(point(i));
    }
}

void COCOBiFunction::fTasy(VectorXd &point, double beta){
	for(int i = 0; i < d_; i++){
		if(point(i) <= 0){continue;}
		double power;
		if(d_ == 1){power = 1;}
		else{power = 1.0 + beta * sqrt(point(i)) * i / (d_ - 1.0);}
		point(i) = pow(point(i), power);
	}
}

void COCOBiFunction::fTalpha(VectorXd &point, double alpha){
	for(int i = 0; i < d_; i++){
		double power;
		if(d_ == 1){power = 0.0;}
		else{power = i / (2.0*(d_ - 1.0));}
		point(i) = point(i) * pow(alpha, power);
	}
}

double COCOBiFunction::fPen(VectorXd &point){
	double result = 0.0;
	for(int i = 0; i < d_; i++){
		if(abs(point(i)) - 5.0 > 0.0){result += pow(abs(point(i)) - 5.0, 2.0);}
		
	}
	return result;
}

int COCOBiFunction::sign(double val){
	if(abs(val) < TOL){return 0;}
	if(val > 0){return 1;}
	return -1;
}



double COCOBiFunction::f1(VectorXd &point){
	VectorXd z = point - xOpt_;
	return pow(z.norm(), 2.0) + fOpt_;
}

double COCOBiFunction::f2(VectorXd &point){
	VectorXd z = point - xOpt_;
	fTosc(z);
	double result = 0;
	for(int i = 0; i < d_; i++){
		if(d_ == 1){result += z(i) * z(i);}
		else{result += pow(10.0, 6.0 * i / (d_ - 1.0)) * z(i) * z(i);}
	}
	result += fOpt_;
	return result;
}

double COCOBiFunction::f3(VectorXd &point){
	VectorXd z = point - xOpt_;
	fTosc(z);
	fTasy(z, 0.2);
	fTalpha(z, 10.0);
	double result = 0.0;
	for(int i = 0; i < d_; i++){
		result += cos(2.0 * M_PI * z(i));
	}
	result = 10.0 * (d_ - result) + pow(z.norm(), 2) + fOpt_;

	return result;
}

double COCOBiFunction::f4(VectorXd &point){
	VectorXd z = point - xOpt_;
	fTosc(z);
	for(int i = 0; i < d_; i++){
		double s;
		if(d_ == 1){s = 1;}
		else{s = pow(sqrt(10.0), i / (d_ - 1.0));}
		if(z(i) > 0 && (i % 2) == 0){
			s = s * 10.0;
		}
		z(i) = s * z(i);
	}
	double result = 0.0;
	double sum2 = 0.0;
	for(int i = 0; i < d_; i++){
		result += cos(2.0 * M_PI * z(i));
		sum2 += z(i) * z(i);
	}
	result = 10.0 * (d_ - result) + sum2 + 100.0 * fPen(point) + fOpt_;
	return result;
}

double COCOBiFunction::f5(VectorXd &point){
	VectorXd z = point;
	for(int i = 0; i < d_; i++){
		if(z(i) * xOpt_(i) > 25){z(i) = xOpt_(i);}
	}
	double result = 0.0;
	for(int i = 0; i < d_; i++){
		double s;
		if(d_ == 1){s = sign(xOpt_(i));}
		else{s = sign(xOpt_(i)) * pow(10.0, i / (d_ - 1.0));}
		result += 5.0 * abs(s) - s * z(i);
	}
	result += fOpt_;
	return result;
}

double COCOBiFunction::f6(VectorXd &point){
	VectorXd z = leftMultiplication_ * (point - xOpt_);
	double result = 0.0;
	for(int i = 0; i < d_; i++){
		double s = 1.0;
		if(z(i) * xOpt_(i) > 0){s = 100.0;}
		result += pow(s * z(i), 2.0);
	}
	result = fTosc(result);
	result = pow(result, 0.9) + fOpt_;
	return result;
}

double COCOBiFunction::f7(VectorXd &point){
	VectorXd zHat = leftMultiplication_ * (point - xOpt_);
	VectorXd zBar = zHat;
	for(int i = 0; i < d_; i++){
		if(abs(zBar(i)) > 0.5){zBar(i) = round(zBar(i));}
		else{zBar(i) = round(10.0 * zBar(i)) / 10.0;}
	}
	VectorXd z = rotationQ_ * zBar;
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		if(d_ == 1){sum += z(i) * z(i);}
		else{sum += pow(10.0, 2.0 * i / (d_ - 1.0)) * z(i) * z(i);}
	}
	return 0.1 * max(pow(10, -4) * zHat(0), sum) + fPen(point) + fOpt_;
}

double COCOBiFunction::f8(VectorXd &point){
	VectorXd z = max(1.0, sqrt(d_)/8.0) * (point - xOpt_);
	for(int i = 0; i < d_; i++){z(i) = z(i) + 1;}

	double result = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		result += 100 * pow(z(i) * z(i) - z(i+1), 2.0) + pow(z(i) - 1.0, 2.0);
	}
	result += fOpt_;
	return result;
}

double COCOBiFunction::f9(VectorXd &point){
	VectorXd z = max(1.0, sqrt(d_)/8.0) * rotationR_ * point;
	for(int i = 0; i < d_; i++){z(i) = z(i) + 0.5;}
	double result = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		result += 100 * pow(z(i) * z(i) - z(i+1), 2.0) + pow(z(i) - 1.0, 2.0);
	}
	result += fOpt_;
	return result;
}

double COCOBiFunction::f10(VectorXd &point){
	VectorXd z = rotationR_ * (point - xOpt_);
	fTosc(z);
	double result = 0.0;
	for(int i = 0; i < d_; i++){
		if(d_ == 1){result += z(i) * z(i);}
		else{result += pow(10.0, 6.0 * i / (d_ - 1.0)) * z(i) * z(i);}
	}
	result += fOpt_;
	return result;
}

double COCOBiFunction::f11(VectorXd &point){
	VectorXd z = rotationR_ * (point - xOpt_);
	fTosc(z);
	double result = pow(10.0, 6.0) * z(0) * z(0);
	for(int i = 1; i < d_; i++){
		result += z(i) * z(i);
	}
	result += fOpt_;
	return result;
}

double COCOBiFunction::f12(VectorXd &point){
	VectorXd zTemp = rotationR_ * (point - xOpt_);
	fTasy(zTemp, 0.5);
	VectorXd z = rotationR_ * zTemp;
	double result = z(0) * z(0);
	for(int i = 1; i < d_; i++){
		result += pow(10.0, 6.0) * z(i) * z(i);
	}
	// result = result / pow(10.0, 3.0);
	result += fOpt_;
	return result;
}

double COCOBiFunction::f13(VectorXd &point){
	VectorXd z = leftMultiplication_ * (point - xOpt_);
	double sum = 0.0;
	for(int i = 1; i < d_; i++){
		sum += z(i) * z(i);
	}
	double result = z(0) * z(0) + 100 * sqrt(sum) + fOpt_;
	return result;
}

double COCOBiFunction::f14(VectorXd &point){
	VectorXd z = rotationR_ * (point - xOpt_);
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		if(d_ == 1){sum += pow(abs(z(i)), 2.0);}
		else{sum += pow(abs(z(i)), 2.0 + 4.0 * i / (d_ - 1.0));}
	}
	double result = sqrt(sum) + fOpt_;
	return result;
}

double COCOBiFunction::f15(VectorXd &point){
	VectorXd zTemp = rotationR_ * (point - xOpt_);
	fTosc(zTemp);
	fTasy(zTemp, 0.2);
	VectorXd z = leftMultiplication_ * zTemp;
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += cos(2 * M_PI * z(i));
	}
	double result = 10.0 * (d_ - sum) + pow(z.norm(), 2) + fOpt_;
	return result;
}

double COCOBiFunction::f16(VectorXd &point){
	VectorXd zTemp = rotationR_ * (point - xOpt_);
	fTosc(zTemp);
	VectorXd z = leftMultiplication_ * zTemp;
	double f0 = 0.0;
	for(int k = 0; k <= 11; k++){
		f0 += pow(0.5, k) * cos(2.0 * M_PI * pow(3, k) * 0.5);
	}
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		for(int k = 0; k <= 11; k++){
			sum += pow(0.5, k) * cos(2.0 * M_PI * pow(3, k) * (z(i) + 0.5));
		}
	}
	double result = 10 * pow(sum / d_ - f0, 3.0) + 10.0 * fPen(point) / (double)d_ + fOpt_;
	return result;
}

double COCOBiFunction::f17(VectorXd &point){
	VectorXd zTemp = rotationR_ * (point - xOpt_);
	fTasy(zTemp, 0.5);
	VectorXd z = leftMultiplication_ * zTemp;
	double sum = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		double s = sqrt(z(i) * z(i) + z(i+1) * z(i+1));
		sum += sqrt(s) + sqrt(s) * pow(sin(50 * pow(s, 0.2)), 2.0);
	}
	double result;
	if(d_ == 1){result = 10 * fPen(point) + fOpt_;}
	else{result = pow(sum / (d_ - 1.0), 2.0) + 10 * fPen(point) + fOpt_;}
	return result;
}

double COCOBiFunction::f18(VectorXd &point){
	VectorXd zTemp = rotationR_ * (point - xOpt_);
	fTasy(zTemp, 0.5);
	VectorXd z = leftMultiplication_ * zTemp;
	double sum = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		double s = sqrt(z(i) * z(i) + z(i+1) * z(i+1));
		sum += sqrt(s) + sqrt(s) * pow(sin(50 * pow(s, 0.2)), 2.0);
	}
	double result;
	if(d_ == 1){result = 10 * fPen(point) + fOpt_;}
	else{result = pow(sum / (d_ - 1.0), 2.0) + 10 * fPen(point) + fOpt_;}
	return result;
}

double COCOBiFunction::f19(VectorXd &point){
	VectorXd z = max(1.0, sqrt(d_)/8.0) * rotationR_ * point;
	for(int i = 0; i < d_; i++){z(i) = z(i) + 0.5;}

	double sum = 0.0;
	for(int i = 0; i < d_ - 1; i++){
		double s = 100 * pow(z(i) * z(i) - z(i+1), 2.0) + pow(z(i) - 1.0, 2.0);
		sum += s / 4000.0 - cos(s);
	}
	double result;
	if(d_ == 1){result = fOpt_;}
	else{result = 10.0 * sum / (d_ - 1.0) + 10.0 + fOpt_;}
	return result;
}

double COCOBiFunction::f20(VectorXd &point){
	VectorXd xHat = 2 * point;
	for(int i = 0; i < d_; i++){
		if(xOpt_(i) < 0){xHat(i) = -1 * xHat(i);}
	}
	VectorXd zHat = xHat;
	for(int i = 1; i < d_; i++){
		zHat(i) = xHat(i) + 0.25 * (xHat(i-1) - 2.0 * abs(xOpt_(i-1)));
	}
	for(int i = 0; i < d_; i++){
		zHat(i) -= 2.0 * abs(xOpt_(i));
	}
	VectorXd z = matrixTalpha(10.0) * zHat;
	for(int i = 0; i < d_; i++){
		z(i) += 2.0 * abs(xOpt_(i));
	}
	z = 100.0 * z;
	// VectorXd z = 100.0 * (matrixTalpha(10.0) * (zHat - 2.0 * abs(xOpt_) + 2.0 * abs(xOpt_)));
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		sum += z(i) * sin(sqrt(abs(z(i))));
	}
	VectorXd zDiv = z / 100.0;
	double result = -1 * sum / (100.0 * d_) + 4.189828872724339 + 100.0 * fPen(zDiv) + fOpt_;
	return result;
}

double COCOBiFunction::f21(VectorXd &point){
	double max = -DBL_MAX;
	for(int i = 0; i < 101; i++){
		VectorXd diff = point - localOptima_[i];
		VectorXd matrixMult = diff.transpose() * rotationR_.transpose() * matrixC(i) * rotationR_ * diff;
		double w;
		if(i == 0){w = 10.0;}
		else{w = 1.1 + 8 * (i - 1.0) / 99.0;}
		double result = w * exp(-1.0 * matrixMult(0) / (2.0 * d_));
		if(result > max){max = result;}
	}
	double total = pow(fTosc(10 - max), 2.0) + fPen(point) + fOpt_;
	return total;
}

double COCOBiFunction::f22(VectorXd &point){
	double max = -DBL_MAX;

	for(int i = 0; i < 21; i++){
		VectorXd diff = point - localOptima_[i];
		VectorXd matrixMult = diff.transpose() * rotationR_.transpose() * matrixC(i) * rotationR_ * diff;
		double w;
		if(i == 0){w = 10.0;}
		else{w = 1.1 + 8 * (i - 1.0) / 19.0;}
		double result = w * exp(-1.0 * matrixMult(0) / (2.0 * d_));
		if(result > max){max = result;}
	}
	double total = pow(fTosc(10 - max), 2.0) + fPen(point) + fOpt_;
	return total;
}

double COCOBiFunction::f23(VectorXd &point){
	VectorXd z = leftMultiplication_ * (point - xOpt_);
	double prod = 1.0;
	for(int i = 0; i < d_; i++){
		double sum = 0.0;
		for(int j = 1; j <= 32; j++){
			sum += abs(pow(2.0, j) * z(i) - round(pow(2.0, j) * z(i)) ) / pow(2.0, j);
		}
		prod *= (1.0 + (i + 1.0) * sum);
	}
	double result = 10.0 * pow(prod, 10.0 / pow(d_, 1.2)) / (d_ * d_) - 10.0 / (d_ * d_) + fPen(point) + fOpt_;
	return result;
}

double COCOBiFunction::f24(VectorXd &point){
	// Define some constants
	double d = 1.0;
	double s = 1.0 - 1.0 / (2 * sqrt(d_ + 20.0) - 8.2);
	double mu0 = 2.5;
	double mu1 = -1.0 * sqrt((mu0 * mu0 - d)/s);

	VectorXd xHat = 2 * point;
	VectorXd zTemp = point;
	for(int i = 0; i < d_; i++){
		if(xOpt_(i) < 0){xHat(i) = -1.0 * xHat(i);}
		zTemp(i) = xHat(i) - mu0; 
	}
	VectorXd z = leftMultiplication_ * zTemp;
	double leftMin = 0.0;
	double rightMin = 0.0;
	double sum = 0.0;
	for(int i = 0; i < d_; i++){
		leftMin += pow(xHat(i) - mu0, 2.0);
		rightMin += pow(xHat(i) - mu1, 2.0);
		sum += cos(2.0 * M_PI * z(i));
	}
	rightMin = d * d_ + s * rightMin;
	double result = min(leftMin, rightMin) + 10.0 * (d_ - sum) + 10000 * fPen(point) + fOpt_;
	return result;
}



double COCOBiFunction::evaluate(VectorXd &point){
	if(function_ == 1){return f1(point);}
	else if(function_ == 2){return f2(point);}
	else if(function_ == 3){return f3(point);}
	else if(function_ == 4){return f4(point);}
	else if(function_ == 5){return f5(point);}
	else if(function_ == 6){return f6(point);}
	else if(function_ == 7){return f7(point);}
	else if(function_ == 8){return f8(point);}
	else if(function_ == 9){return f9(point);}
	else if(function_ == 10){return f10(point);}
	else if(function_ == 11){return f11(point);}
	else if(function_ == 12){return f12(point);}
	else if(function_ == 13){return f13(point);}
	else if(function_ == 14){return f14(point);}
	else if(function_ == 15){return f15(point);}
	else if(function_ == 16){return f16(point);}
	else if(function_ == 17){return f17(point);}
	else if(function_ == 18){return f18(point);}
	else if(function_ == 19){return f19(point);}
	else if(function_ == 20){return f20(point);}
	else if(function_ == 21){return f21(point);}
	else if(function_ == 22){return f22(point);}
	else if(function_ == 23){return f23(point);}
	else if(function_ == 24){return f24(point);}
	return 0.0;
}



double COCOBiFunction::evaluateLow(VectorXd &point){
	// Check noise came out alright
	double value = evaluate(point);
	bool addNoise = false;
	double dist = 0.0;
	double mult = 0.0;
	if(isnan(value)){
		printf("Encountered nan value from high fi at point (");
		for(int i = 0; i < (int)point.size() - 1; i++){
			printf("%.2f,", point(i));
			
		}
		printf("%.2f)! Stopping now...\n", point((int)point.size() - 1));
		exit(0);
	}
	// Work out value of sine and cosine sine noise, will work out a global position
	// and then add it if required
	VectorXd minDist(d_);
	for(int i = 0; i < d_; i++){minDist(i) = lowerBound_[i];}
	double trav = (point - minDist).norm() / maxDist_;
	// Calculate whether noise needs to be added, and its impact (mult value)
	if(disturbanceType_ == 'h' && disturbanceNum_ == 1){
		// Height based
		double noiseCentre = fOpt_ + disturbanceHeight_ * (fMax_ - fOpt_);
		double absRadius = disturbanceRadius_ * (fMax_ - fOpt_);
		dist = abs(value - noiseCentre);
		if(dist < absRadius){
			addNoise = true;
			mult = 1.0 - dist / absRadius;
		}

	}else if(disturbanceType_ == 'h' && disturbanceNum_ == 2){
		// Height based
		double noiseCentre = fOpt_ + disturbanceHeight_ * (fMax_ - fOpt_);
		double absRadius = disturbanceRadius_ * (fMax_ - fOpt_);
		dist = abs(value - noiseCentre);
		if(dist > absRadius){
			addNoise = true;
			mult = 1.0 - ((fMax_ - fOpt_) - dist) / ((fMax_ - fOpt_) - absRadius);
		}

	}else if(disturbanceType_ == 's' && disturbanceNum_ == 1){
		// Noise centres, work out closest point
		dist = DBL_MAX;
		for(int i = 0; i < (int)disturbanceCentres_.size(); i++){
			if((disturbanceCentres_[i] - point).norm() < dist){dist = (disturbanceCentres_[i] - point).norm();}
		}
		double absRadius = disturbanceRadius_ * maxDist_;
		if(dist > absRadius){
			addNoise = true;
			mult = 1.0 - (maxDist_ - dist) / (maxDist_ - absRadius);
		}
	}else if(disturbanceType_ == 's' && disturbanceNum_ == 2){
		// Disturbance centres, work out closest point
		dist = DBL_MAX;
		for(int i = 0; i < (int)disturbanceCentres_.size(); i++){
			if((disturbanceCentres_[i] - point).norm() < dist){dist = (disturbanceCentres_[i] - point).norm();}
		}
		double absRadius = disturbanceRadius_ * maxDist_;
		if(dist < absRadius){
			addNoise = true;
			mult = 1.0 - dist / absRadius;
		}
	
	}else{
		printf("Asking for COCO low fidelity evaluation with undefined pair (%c,%d) of disturbance values! Stopping now...\n", disturbanceType_, disturbanceNum_);
		exit(0);
	}

	// Add disturbance if it applies
	if(addNoise){
		value += mult * addBasicDisturbance(trav);
	}
	if(isnan(value)){
		printf("Encountered nan value at point (");
		for(int i = 0; i < (int)point.size() - 1; i++){
			printf("%.2f,", point(i));
			
		}
		printf("%.2f)! Stopping now...\n", point((int)point.size() - 1));
		exit(0);
	}
	return value;

}

double COCOBiFunction::addBasicDisturbance(double value){
	return basicDisturbanceAmplitude_ * (fMax_ - fOpt_) * cos(basicDisturbanceFrequency_ * 2 * M_PI * value) * sin(basicDisturbanceFrequency_ * 2 * M_PI * value * value);
	// return amplitude * cos(100 * 2 * M_PI * value) * sin(100 * 2 * M_PI * value * value);
}









#endif