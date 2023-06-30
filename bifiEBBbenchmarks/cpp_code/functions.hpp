#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "libraries.hpp"

// Definition of different bi-fidelity functions found in the literature.
// References to papers where the instances are defined are given at the beginning of each set of functions,
// but the whole list is first given here in order of appearance:

// 		- Toal DJ (2015) "Some considerations regarding the use of multi-fidelity kriging in the construction of surrogate models."
// 		- Song X, Lv L, Sun W, Zhang J (2019) "A radial basis function-based multi-fidelity surrogate model: exploring
// 				correlation between high-fidelity and low-fidelity models."
// 		- March A, Willcox K (2012) "Provably convergent multifidelity optimization algorithm not requiring high-fidelity derivatives"
// 		- Rajnarayan D, Haas A, Kroo I (2008) "A multifidelity gradient-free optimization method and application to aerodynamic design".
//		- Liu B, Koziel S, Zhang Q (2016) "A multi-fidelity surrogate-model-assisted evolutionary algorithm for
//				computationally expensive optimization problems"
//		- Wang H, Jin Y, Doherty J (2017) "A generic test suite for evolutionary multifidelity optimization."
//		- Liu H, Ong YS, Cai J, Wang Y (2018b) "Cope with diverse data structures in multi-fidelity modeling: a
//				gaussian process method."
//		- Wu Y, Hu J, Zhou Q, Wang S, Jin P (2020) "An active learning multi-fidelity metamodeling method based
// 				on the bootstrap estimator."
//		- Shi M, Lv L, Sun W, Song X (2020) "A multi-fidelity surrogate model based on support vector regression"
// 		- Dong H, Song B, Wang P, Huang S (2015) "Multi-fidelity information fusion based on prediction of kriging."
//		- Xiong S, Qian PZ, Wu CJ (2013) "Sequential design and analysis of high-accuracy and low-accuracy computer codes."
// 		- Park C, Haftka RT, Kim NH (2017) "Remarks on multi-fidelity surrogates".
//		- Hansen N, Auger A, Ros R, Mersmann O, Tuˇsar T, Brockhoff D (2020) "COCO: A platform for comparing
//				continuous optimizers in a black-box setting"
// 		- Andres-Thio N, Munoz MA, Smith-Miles K (2022): "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances"




// Parent Function  class. This class should not be instantiated, it is used as a template for future functions.
// Note sample space is assumed to be a hypercube.
class Function {
	public:

	// Constructor for function for which upper and lower bounds are the same for all dimensions.
	Function(int d, double lowBound, double upBound);

	// Constructor for function with different upper and lower bound for each dimensions; bounds are specified as vectors.
	Function(int d, vector<double> lowerBound, vector<double> upperBound);

	virtual ~Function();

	// Checks whether a point has the same dimension as the sample space of the function.
	void checkDimension(VectorXd &point);

	// Check whether a point is within the sample space.
	bool pointWithinBounds(VectorXd &point);

	// Check whether a point is within the sample space for a particular dimsneions.
	int dimensionWithinBounds(VectorXd &point, int d);

	// Main call to function, returns function value at point.
	virtual double evaluate(VectorXd &point);

	// Option to ask for objective valof multiple points at once.
	vector<double> evaluateMany(vector<VectorXd> &points);

	// Function which compares objective function at two points; it returns 0 if the value is equal, -1 if point1 is smaller, and 1 if point2 is smaller
	virtual int betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag = "");

	int betterPoint(VectorXd &point1, VectorXd &point2, string flag = "");

	int d_;							// Dimension of sample space of function i.e. number of input variables.
	vector<double> lowerBound_;		// Lower bound of sample space.
	vector<double> upperBound_;		// Upper bound of sample space.
};



// Parent bi-fidelity function class. This class should not be instantiated, it is used as a template.
class BiFidelityFunction : public Function {
	public:

	// Constructor for function for which upper and lower bounds are the same for all dimensions.
	BiFidelityFunction(int d, double lowBound, double upBound);
		
	// Constructor for function with different upper and lower bound for each dimensions; bounds are specified as vectors.
	BiFidelityFunction(int d, vector<double> lowerBound, vector<double> upperBound);

	virtual ~BiFidelityFunction();

	// Need to override this function but it should still do nothing.
	virtual double evaluate(VectorXd &point) override;

	// Returns value of low fidelity source at a point.
	virtual double evaluateLow(VectorXd &point);

	// Returns value of low fidelity source at many points.
	vector<double> evaluateManyLow(vector<VectorXd> &points);
};

// Function which queries the SOLAR simulation engine https://github.com/bbopt/solar
class SOLARFunction : public BiFidelityFunction {
	public:
	SOLARFunction(double fidelityLevel, int fileNum = 0);
	~SOLARFunction();
	VectorXd scalePoint(VectorXd point);
	double callSimulation(VectorXd &inputPoint, double fidelity);
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double fidelityLevel_;
	int fileNum_;
};


// Function identical to SongToalForretal, but of higher dimension.
// Dimensions other than d = 1 have no impact on the function value.
class NicoSongToalForretalFunction : public BiFidelityFunction {
	public:
	NicoSongToalForretalFunction(int dim, double a);
	~NicoSongToalForretalFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};



// The following functions are defined in 
// Toal DJ (2015) "Some considerations regarding the use of multi-fidelity kriging in the construction of surrogate models."
// They define a family of pairs of functions based on input parameter a.
class ToalBraninFunction : public BiFidelityFunction {
	public:
	ToalBraninFunction(double a);
	~ToalBraninFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};
class ToalPaciorekFunction : public BiFidelityFunction {
	public:
	ToalPaciorekFunction(double a);
	~ToalPaciorekFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};
class ToalHartmannH3Function : public BiFidelityFunction {
	public:
	ToalHartmannH3Function(double a);
	~ToalHartmannH3Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
	vector<double> alpha_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};
class ToalTridFunction : public BiFidelityFunction {
	public:
	ToalTridFunction(double a);
	~ToalTridFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};


// The following functions are defined in 
// Song X, Lv L, Sun W, Zhang J (2019) "A radial basis function-based multi-fidelity surrogate model: exploring
// correlation between high-fidelity and low-fidelity models."
// They define a family of pairs of functions based on input parameter a.
class SongToalForretalFunction : public BiFidelityFunction {
	public:
	SongToalForretalFunction(double a);
	~SongToalForretalFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};
class SongToalBraninFunction : public BiFidelityFunction {
	public:
	SongToalBraninFunction(double a);
	~SongToalBraninFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};
class SongToalColvilleFunction : public BiFidelityFunction {
	public:
	SongToalColvilleFunction(double a);
	~SongToalColvilleFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	double a_;
};

// The following functions are defined in 
// March A, Willcox K (2012) "Provably convergent multifidelity optimization algorithm not requiring high-fidelity derivatives"
// A single high fidelity function is paired with 5 different low fidelity functions.
class MarchWillcoxRosenbrockFunction : public BiFidelityFunction {
	public:
	MarchWillcoxRosenbrockFunction(int lowFiFunction);
	~MarchWillcoxRosenbrockFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	int lowFiFunction_;
};


// The following functions are defined in 
// Rajnarayan D, Haas A, Kroo I (2008) "A multifidelity gradient-free optimization method and application to aerodynamic design".
class RajnarayanHartmannH3Function : public BiFidelityFunction {
	public:
	RajnarayanHartmannH3Function();
	~RajnarayanHartmannH3Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> alpha_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};
class RajnarayanHartmannH6Function : public BiFidelityFunction {
	public:
	RajnarayanHartmannH6Function();
	~RajnarayanHartmannH6Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> alpha_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};
class RajnarayanWoodsFunction : public BiFidelityFunction {
	public:	
	RajnarayanWoodsFunction();
	~RajnarayanWoodsFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};


// The following functions are defined in 
// Liu B, Koziel S, Zhang Q (2016) "A multi-fidelity surrogate-model-assisted evolutionary algorithm for 
// computationally expensive optimization problems"
class LiuEllipsoidFunction : public BiFidelityFunction {
	public:	
	LiuEllipsoidFunction();
	~LiuEllipsoidFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> sH_;
	vector<double> sS_;
};
class LiuDixonPriceFunction : public BiFidelityFunction {
	public:	
	LiuDixonPriceFunction();
	~LiuDixonPriceFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> sS_;
};
class LiuStyblinskiTangFunction : public BiFidelityFunction {
	public:
	LiuStyblinskiTangFunction();
	~LiuStyblinskiTangFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> sS_;
};
class LiuAckley10Function : public BiFidelityFunction {
	public:
	LiuAckley10Function();
	~LiuAckley10Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> sS_;
	double sF_;
};
class LiuAckley20Function : public BiFidelityFunction {
	public:	
	LiuAckley20Function();
	~LiuAckley20Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> sS_;
	double sF_;
};



// The following functions are defined in 
// Wang H, Jin Y, Doherty J (2017) "A generic test suite for evolutionary multifidelity optimization."
// Using a single high fidelity function for which the dimension can be specified,
// a family of errors can be added to create a low fidelity function.
class WangRastriginFunction : public BiFidelityFunction {
	public:	
	WangRastriginFunction(int dimension, int error, double phi);
	~WangRastriginFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	int error_;
	double phi_;
	double aVal_;
	double wVal_;
	double bVal_;
	double thetaVal_;
};



// The following functions are defined in 
// Liu H, Ong YS, Cai J, Wang Y (2018b) "Cope with diverse data structures in multi-fidelity modeling: a
// gaussian process method."
class LiuPedagogicalFunction : public BiFidelityFunction {
	public:	
	LiuPedagogicalFunction();
	~LiuPedagogicalFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class LiuBraninFunction : public BiFidelityFunction {
	public:
	LiuBraninFunction();
	~LiuBraninFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};



// The following function is defined in
// Wu Y, Hu J, Zhou Q, Wang S, Jin P (2020) "An active learning multi-fidelity metamodeling method based
// on the bootstrap estimator."
class WuHimmelblauFunction : public BiFidelityFunction {
	public:
	WuHimmelblauFunction();
	~WuHimmelblauFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};



// The following functions are defined in
// Shi M, Lv L, Sun W, Song X (2020) "A multi-fidelity surrogate model based on support vector regression"
class ShiGramacyLeeFunction : public BiFidelityFunction {
	public:
	ShiGramacyLeeFunction();
	~ShiGramacyLeeFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiCurrinSinFunction : public BiFidelityFunction {
	public:
	ShiCurrinSinFunction();
	~ShiCurrinSinFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiHolsclawFunction : public BiFidelityFunction {
	public:
	ShiHolsclawFunction();
	~ShiHolsclawFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiSantnerFunction : public BiFidelityFunction {
	public:
	ShiSantnerFunction();
	~ShiSantnerFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiBraninFunction : public BiFidelityFunction {
	public:
	ShiBraninFunction();
	~ShiBraninFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiNumberSixFunction : public BiFidelityFunction {
	public:
	ShiNumberSixFunction();
	~ShiNumberSixFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiNumberSevenFunction : public BiFidelityFunction {
	public:
	ShiNumberSevenFunction();
	~ShiNumberSevenFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiBealeFunction : public BiFidelityFunction {
	public:
	ShiBealeFunction();
	~ShiBealeFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiStyblinskiTangFunction : public BiFidelityFunction {
	public:
	ShiStyblinskiTangFunction();
	~ShiStyblinskiTangFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiCurrinExpFunction : public BiFidelityFunction {
	public:
	ShiCurrinExpFunction();
	~ShiCurrinExpFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiLimFunction : public BiFidelityFunction {
	public:
	ShiLimFunction();
	~ShiLimFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiGramacyFunction : public BiFidelityFunction {
	public:
	ShiGramacyFunction();
	~ShiGramacyFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiChengSanduFunction : public BiFidelityFunction {
	public:
	ShiChengSanduFunction();
	~ShiChengSanduFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiHartmannH3Function : public BiFidelityFunction {
	public:
	ShiHartmannH3Function();
	~ShiHartmannH3Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> alpha_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};
class ShiDettePepelyshevExpFunction : public BiFidelityFunction {
	public:
	ShiDettePepelyshevExpFunction();
	~ShiDettePepelyshevExpFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiHartmannH4Function : public BiFidelityFunction {
	public:
	ShiHartmannH4Function();
	~ShiHartmannH4Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> alpha_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};
class ShiParkFunction : public BiFidelityFunction {
	public:
	ShiParkFunction();
	~ShiParkFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiHartmannH6Function : public BiFidelityFunction {
	public:
	ShiHartmannH6Function();
	~ShiHartmannH6Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> alpha_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};
class ShiRosenbrockFunction : public BiFidelityFunction {
	public:
	ShiRosenbrockFunction();
	~ShiRosenbrockFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class ShiDettePepelyshevFunction : public BiFidelityFunction {
	public:
	ShiDettePepelyshevFunction();
	~ShiDettePepelyshevFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};



// The following functions are defined in
// Dong H, Song B, Wang P, Huang S (2015) "Multi-fidelity information fusion based on prediction of kriging."
class DongBohachevskyFunction : public BiFidelityFunction {
	public:
	DongBohachevskyFunction();
	~DongBohachevskyFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class DongBoothFunction : public BiFidelityFunction {
	public:
	DongBoothFunction();
	~DongBoothFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class DongBraninFunction : public BiFidelityFunction {
	public:
	DongBraninFunction();
	~DongBraninFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class DongHimmelblauFunction : public BiFidelityFunction {
	public:
	DongHimmelblauFunction();
	~DongHimmelblauFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class DongSixHumpCamelbackFunction : public BiFidelityFunction {
	public:
	DongSixHumpCamelbackFunction();
	~DongSixHumpCamelbackFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};



// The following functions are defined in
// Xiong S, Qian PZ, Wu CJ (2013) "Sequential design and analysis of high-accuracy and low-accuracy computer codes."
class XiongCurrinExpFunction : public BiFidelityFunction {
	public:
	XiongCurrinExpFunction();
	~XiongCurrinExpFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class XiongParkFirstFunction : public BiFidelityFunction {
	public:
	XiongParkFirstFunction();
	~XiongParkFirstFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};
class XiongParkSecondFunction : public BiFidelityFunction {
	public:
	XiongParkSecondFunction();
	~XiongParkSecondFunction();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
};



// The following function is defined in
// Park C, Haftka RT, Kim NH (2017) "Remarks on multi-fidelity surrogates".
class ParkHartmannH6Function : public BiFidelityFunction {
	public:
	ParkHartmannH6Function();
	~ParkHartmannH6Function();
	virtual double evaluate(VectorXd &point) override;
	virtual double evaluateLow(VectorXd &point) override;
	vector<double> alpha_;
	vector<double> alphaDash_;
	vector< vector<double> > beta_;
	vector< vector<double> > p_;
};






class disturbanceBasedBiFunction : public BiFidelityFunction{
public:
	disturbanceBasedBiFunction(int function, int seed);

	~disturbanceBasedBiFunction();

	// Function which calls the appropriate definition based on the stored function_ value.
	virtual double evaluate(VectorXd &point) override;

	// Implementes the addition of a disturbance.
	virtual double evaluateLow(VectorXd &point) override;

	// Basic disturbance used is a trigonometric function.
	double addBasicDisturbance(double value);





	int function_;							// Function instanciated, between 1-24.
	int seed_;								// Random seed used for reproducibility.
	mt19937 randGen_;						// Generator of all random numbers of the class.

	BiFidelityFunction* highFiFunction_;
	
	double fMax_;							// Maxium value of the function. Not initialised, this variable exists so that it can be saved.
	double fMin_;
	double maxDist_;						// Maximum distance from any pair of points. Calculated as the distance between opposite endpoints of sample space hypercube.

	char disturbanceType_;					// Type of Disturbance, can be height based (h) or source based (s).
	int disturbanceNum_;					// Given disturbance type, the number (1 oe 2) specifies how it is applied.

	double disturbanceHeight_;				// Height at which the disturbance is centered, used for disturbance type h.
	double disturbanceRadius_;				// Radius of disturbance, used for both disturbance types.
	vector<VectorXd> disturbanceCentres_;	// Centres of the disturbances for disturbance type s.

	int basicDisturbanceFrequency_;			// Frequency of basic disturbance, affects behaviour of trigonometric functions.
	double basicDisturbanceAmplitude_;		// Amplitude of trigonometric functions in the basic disturbance.
	

};







// COCO functions defined in 
// Hansen N, Auger A, Ros R, Mersmann O, Tuˇsar T, Brockhoff D (2020) "COCO: A platform for comparing
//	continuous optimizers in a black-box setting"
//
// For implementation details, consult the document
// "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions"
//
// Low fidelity function creation is that preseneted in 
// Andres-Thio N, Munoz MA, Smith-Miles K (2022): "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances"
//
class COCOBiFunction : public BiFidelityFunction {
	public:

	COCOBiFunction(int function, int dimension, int seed);

	~COCOBiFunction();

	// Initialises matrices and constants which will be used by some of the functions.
	void initialiseConstants();

	// COCO functions are intialised with randomly chosen optimum location and value.
	// This function randomly chooses these values.
	void chooseOptVals();

	//  Random orthodonal (rotation) matrix creation.
	// As stated in "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions":
	// Orthogonal matrices are generated from standard normally distributed entries by Gram-Schmidt orthonormalization.
	// Columns and rows of an orthogonal matrix form an orthonormal basis
	MatrixXd randomRotationMatrix(int d);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	double fTosc(double val);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	void fTosc(VectorXd &point);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	void fTasy(VectorXd &point, double beta);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	void fTalpha(VectorXd &point, double alpha);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	MatrixXd matrixTalpha(double alpha);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	MatrixXd matrixC(int index);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	double fPen(VectorXd &point);

	// Consult "Real-Parameter Black-Box Optimization Benchmarking 2010: Presentation of the Noiseless Functions".
	int sign(double val);

	// Implementation of the 24 function of the COCO test suite.
	double f1(VectorXd &point);
	double f2(VectorXd &point);
	double f3(VectorXd &point);
	double f4(VectorXd &point);
	double f5(VectorXd &point);
	double f6(VectorXd &point);
	double f7(VectorXd &point);
	double f8(VectorXd &point);
	double f9(VectorXd &point);
	double f10(VectorXd &point);
	double f11(VectorXd &point);
	double f12(VectorXd &point);
	double f13(VectorXd &point);
	double f14(VectorXd &point);
	double f15(VectorXd &point);
	double f16(VectorXd &point);
	double f17(VectorXd &point);
	double f18(VectorXd &point);
	double f19(VectorXd &point);
	double f20(VectorXd &point);
	double f21(VectorXd &point);
	double f22(VectorXd &point);
	double f23(VectorXd &point);
	double f24(VectorXd &point);

	// Function which calls the appropriate definition based on the stored function_ value.
	virtual double evaluate(VectorXd &point) override;

	// Implementes the addition of a disturbance.
	virtual double evaluateLow(VectorXd &point) override;

	// Basic disturbance used is a trigonometric function.
	double addBasicDisturbance(double value);


	int function_;							// Function instanciated, between 1-24.
	int seed_;								// Random seed used for reproducibility.
	mt19937 randGen_;						// Generator of all random numbers of the class.
	MatrixXd rotationQ_;					// Stored orthogonal rotation matrix used by some functions.
	MatrixXd rotationR_;					// Another orthogonal rotation matrix for functions which use two.
	MatrixXd leftMultiplication_;			// Multiplication of matrices, saved for speed when evaluating certain functions.
	vector<double> alphas_;					// Vector used by some functions.
	vector< vector<int> > matrixPerm_;		// Permutation matrix used by some functions.
	vector<VectorXd> localOptima_;			// Set of local optima randomly generated, used by some functions.

	VectorXd xOpt_;							// Randomly generated position of the optimum.
	double fOpt_;							// Randomly generated optimum value, chosen to be between -100 and 100.
	double fMax_;							// Maxium value of the function. Not initialised, this variable exists so that it can be saved.
	double maxDist_;						// Maximum distance from any pair of points. Calculated as the distance between opposite endpoints of sample space hypercube.

	char disturbanceType_;					// Type of Disturbance, can be height based (h) or source based (s).
	int disturbanceNum_;					// Given disturbance type, the number (1 oe 2) specifies how it is applied.

	double disturbanceHeight_;				// Height at which the disturbance is centered, used for disturbance type h.
	double disturbanceRadius_;				// Radius of disturbance, used for both disturbance types.
	vector<VectorXd> disturbanceCentres_;	// Centres of the disturbances for disturbance type s.

	int basicDisturbanceFrequency_;			// Frequency of basic disturbance, affects behaviour of trigonometric functions.
	double basicDisturbanceAmplitude_;		// Amplitude of trigonometric functions in the basic disturbance.
	
};












#endif