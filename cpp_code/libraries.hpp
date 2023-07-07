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

#ifndef LIBRARIES_HPP
#define LIBRARIES_HPP

#include <algorithm>
#include <cfloat>
#include <ctime>
#include <cmath>
// Uncomment this for everything except running rFeatureAnalysis.cpp from windows
// #include <Eigen/Dense>
#include "C:\Users\nandr\Documents\eigen-3.4.0\Eigen\Dense"

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
