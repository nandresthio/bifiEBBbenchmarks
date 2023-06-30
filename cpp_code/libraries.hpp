#ifndef LIBRARIES_HPP
#define LIBRARIES_HPP

#include <algorithm>
#include <cfloat>
#include <ctime>
#include <cmath>
// Uncomment this for everything except running rFeatureAnalysis.cpp from windows
#include <Eigen/Dense>
// #include "C:\Users\nandr\Documents\eigen-3.4.0\Eigen\Dense"

#include <fstream>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <tuple> 
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

#define OBJTOL 0.000001
#define TOL 0.000001
#define NOISE_TOL 0.05

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
	#define M_SQRT1_2 0.707106781186547524401
#endif


#endif
