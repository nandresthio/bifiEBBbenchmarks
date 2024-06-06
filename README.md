# Bi-fidelity Surrogate Modelling benchmark problems

The software in this repository and data linked to it provide multiple benchmarks for Bi-Fidelity Expensive Black Box (Bf-EBB) problems, both for surrogate model building and for optimisation. It implements both classical instances from the literature, as well as the instance generating procedure presented in the paper [Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances](https://doi.org/10.1287/ijoc.2019.0934) by Andrés-Thió N, Muñoz MA, and Smith-Miles K. It also provides instance measures (here called "features") for analysis of the differences between benchmarks, and uses these features to generate objectively varied benchmark sets of different sizes. These benchmark sets can be used by the community when analysing the performance of new Bi-fidelity Expensive Black-Box problems. The paper [Characterising Harmful Data Sources When Constructing Multi-fidelity Surrogate Models](https://arxiv.org/pdf/2403.08118) by Andrés-Thió N, Muñoz MA, and Smith-Miles K details the procedure followed to calculate and standarise the feature values for all instances, and the generation of unbiased benchark test suites of different sizes.

The main benefits of using this piece of code are the following:

- A very large set of candidate functions (over 200 literature instances, and a virtually endless set of newly proposed instances).
- Multiple Latin Hypercube Sampling (LHS) plans of varying size and dimension which are well spread out and interface with the code,
as well as chosen subsets of different sizes. Both sets are locally optimised in terms of the minimum distance between pairs of points.
- A very large set of instance measures, known as features, supplied for all instances, which can be used to analyse the difference between them.
- Instance subsets of different sizes are supplied, chosen to be as varied as possible using an instance filtering procedure.

Note that no surrogate models or algorithms are supplied in this reposoritory, although a different repository which provides an implementation of Kriging and Co-Kriging among others and is run on this set of benchmarks can be found [here](https://github.com/nandresthio/bifiEBBmethods).

## Cite

To cite the usage of this software, please cite the preprint paper using the following BibTex citation:

```
@article{andres2024characterising,
  title={Characterising harmful data sources when constructing multi-fidelity surrogate models},
  author={Andr{\'e}s-Thi{\'o}, Nicolau and Mu{\~n}oz, Mario Andr{\'e}s and Smith-Miles, Kate},
  journal={arXiv preprint arXiv:2403.08118},
  year={2024}
}
} 
```

To cite the software itself, use the following DOI: [![DOI](https://zenodo.org/badge/659941002.svg)](https://zenodo.org/badge/latestdoi/659941002); below is the BibTex for citing this version of the code.

```
@misc{nandres2022,
  author = {Andr\'{e}s-Thi\'{o}, Nicolau},
  title = {Bifidelity Surrogate Modelling benchmark problems},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  note = {available for download at https://github.com/nandresthio/bifiEBBbenchmarks},
  doi = {10.5281/zenodo.8353690}
} 
```


## Description

This code implements test functions for Bi-Fidelity Expensive Black-Box problems. These problems have a high- and low-fidelity source of information denoted $f_h$ and $f_l$ respectively, which are used to either accurately model or to optimise the high-fidelity source. This repository provides implementations of a vast set of functions pairs $(f_h,f_l)$, implemented as derived classes of the class `BiFidelityFunction`. This C++ contains the functions `evaluate` and `evaluateLow` which return the objective function value of $f_h$ and $f_l$ respectively. The function `processFunctionName` is provided for ease, which takes the name of an instance as a string and returns an instanciated class of the relevant instance. The names of the instances implemented so far are given in the folder `data/availableFunctions`. The files `literatureBiSourceDimN.txt` contain the implemented literature instances of dimension $N$, the file `literatureSingleSourceAll.txt` provides the function names derived from unique high-fidelity functions, and the file`disturbanceFunctions.txt` contains the instances defined using the instance generating procedure for which features are available. For an example on how to instanciate a BiFidelityFunction class, relying on already optimised sampling plans, and querying a function pair for high and low fidelity functions, look at the `cpp_code/main.cpp` file. A python interface is also provided in `pythonUsage.ipynb` (which can be opened with a Jupyter notebook) for examples on how to use the implemented instances in a python script. Note the code needs to first be compiled as detailed below.

### Instances

A Bf-EBB instance is made up of a high- and low-fidelity source of information denoted $f_h$ and $f_l$ respectively, which are used to either accurately model or to optimise the high-fidelity source. This repository provides implementations of a vast set of function pairs $(f_h,f_l)$, implemented as derived classes of the class `BiFidelityFunction`. This C++ class contains the functions `evaluate` and `evaluateLow` which return the objective function value of $f_h$ and $f_l$ respectively. The function `processFunctionName` is provided for ease of use, which takes the name of an instance as a string and returns an instanciated class of the relevant instance.

The names of the instances implemented so far are given in the folder `data/availableFunctions`. The files `literatureBiSourceDim*.txt` contain the names of the implemented literature instances of dimension `*`, with the file `literatureBiSourceDimAll.txt` containing all literature test instances. These are generated from 38 unique $f_h$ functions; these are provided in the file `litreatureSingleSourceAll.txt`. Finally, the file `disturbanceFunctions.txt` contains over 80,000 instances generated from both the literature high fidelity functions and the COCO test suite using the disturbance instance generating procedure[^1]. Note that more instances of this type can be generated; this file only conatins a large set of examples.

In addition, code has been implemented to interface with the SOLAR simulation engine, which simulates the functioning of a solar power plant. The code implemented here queries the tenth function of the simulator, which takes a fidelity as input. The high fidelity source is taken to be this source with a fidelity of 1, and the low-fidelity soource the simulation with a fidelity specified when generating an instance of the `BiFidelityFunction` class. For instance, calling `processFunctionName` with the string `SOLAR0.70` creates an instances where calling `evaluateLow` runs the simulator with a fidelity of 0.70.

Finally the files `chosenTestSuiteN*.txt` contain an objectively varied set of instances of size `*`, constructed by using an instance filtering method to remove similar instances. For further details, the reader is referred to the paper which details the construction of these datasets. For further details on the implementation of each of the function pairs and their original publication[^2][^3][^4][^5][^6][^7][^8][^9][^10][^11][^12][^13][^14], please refer to the `Appendix` folder.


### Features


Features are measures of an instance which allows for comparison between different instances. In this repository features are provided which measure the relationship between $f_h$ and $f_l$, and the landscape of $f_h$, $f_l$ and $f_h - f_l$. The relationship features are taken from the work of Toal[^2] and that of Andrés-Thió et al.[^1], whereas the landscape features are taken from the flacco[^15] package from R. Some modifications are applied to the flacco code in the file `R_code/customFeatureCalculation` to calculate features with small samples, or in special cases where the objective function values are all identical.

The features are calculated both with a very large sample, and for a variety of smaller sample sizes. The former is done for all of the implemented function pairs (i.e. roughly 80,000 function pairs) using the file `R_code\featureCalculation.R`, and the later is done only for a selected set of 321 function pairs using the file `R_code\featureCalculationWithActualSample.R`. Note new users should not rerun this, as the results are stored in the files `data\features\features.txt` and `data\features\sampleFeatures.txt` respectively available in [figshare](https://figshare.com/projects/Bi-fidelity_Surrogate_Modelling_benchmark_problems_data/178659). The same features are "cleaned" (i.e. features with NAs or Inf are removed) and made available in `data\features\featuresClean.txt` and `data\features\sampleFeaturesClean.txt` respectively. Finally, a processing of the features is also conducted in order to make them comparable. For features which are unbounded, a box-cox transformation is applied which gives the feature values a normal distribution. The values are then standarised and bounded within [-4,4]. For features which are bound, a linear transformation is applied sp that the new bound is [-2,2]. As roughly 95% of the unbound features should lie within the [-2,2] range, this allows for comparison between features. These standarised feature values are stored in the files `data\features\featuresCleanStandarised.txt` and `data\features\sampleFeaturesCleanStandarised.txt`. The files `data\features\sampleAndRealFeatures.txt`and `data\features\sampleAndRealFeaturesStandarised.txt` contain features calculated with a large and small sample for a selected set of 321 instances.

The files `R_code\standariseFeaturesAndFilterInstances.R` and `R_code\standariseSampleFeatures.R` process the feature values and select benchmark sets of different sizes. For detailed definitions of the features, consult the `Appendix` folder.

### Compilation

Before compiling and running the code, it is recommended to download the populated data folder available from [figshare](https://figshare.com/projects/Bi-fidelity_Surrogate_Modelling_benchmark_problems_data/178659). The folder structure of the data folder in this repository is made so that all the code should run without issues, but in order to replicate results and rely on already optimised sampling plans it is recommended to replace the given data folder with the one available in figshare. 

The code related to the SOLAR simulation engine is not provided here, only the code with interfaces with it. If these type of instances are required, add the SOLAR code following these steps:
  
  - From the root directory, call `cd cpp_code` to move to this folder
  - Call `git clone https://github.com/bbopt/solar.git` to get the SOLAR code
  - Call `cd solar/src/`
  - Call `make` 
  - Call `cd ..`
  - Call `bin\solar -check`

If all goes well, SOLAR should be installed and the executable ready to be called from the software provided here.

To compile the C++ code given here, first make sure the library [Eigen has been installed](https://eigen.tuxfamily.org/index.php?title=Main_Page). Then go to the relevant makefile (i.e. `cpp_code\Makefile.linux` or `cpp_code\Makefile.windows`) and edit the `CFLAGS` to point to the right folder (which contains the installed Eigen library). To compile the C++ code, simply move to the `cpp_code` directory and call `make linux` if using Linux, or `make windows` if using windows. Note all code (both the R files using Rstudio or Rscript, or calling the c++ executable) should be called from the root folder due to the folder structure.



[^1]: Andrés-Thió N, Muñoz M A, Smith-Miles K (2022), "Bifidelity surrogate modelling: Showcasing the need for new test instances", INFORMS Journal on Computing 34 (6) 3007–3022.
[^2]: Toal DJ (2015) "Some considerations regarding the use of multi-fidelity kriging in the construction of surrogate models". Structural and Multidisciplinary Optimization 51(6):1223–1245.
[^3]: Song X, Lv L, Sun W, Zhang J (2019) "A radial basis function-based multi-fidelity surrogate model: exploring correlation between high-fidelity and low-fidelity models". Structural and Multidisciplinary Optimization 60(3):965–981.
[^4]: March A, Willcox K (2012) "Provably convergent multifidelity optimization algorithm not requiring high-fidelity derivatives". AIAA journal 50(5):1079–1089.
[^5]: Rajnarayan D, Haas A, Kroo I (2008) "A multifidelity gradient-free optimization method and application to aerodynamic design". 12th AIAA/ISSMO multidisciplinary analysis and optimization conference, 6020.
[^6]: Liu B, Koziel S, Zhang Q (2016) "A multi-fidelity surrogate-model-assisted evolutionary algorithm for computationally expensive optimization problems". Journal of computational science 12:28–37.
[^7]: Wang H, Jin Y, Doherty J (2017) "A generic test suite for evolutionary multifidelity optimization". IEEE Transactions on Evolutionary Computation 22(6):836–850.
[^8]: Liu H, Ong YS, Cai J, Wang Y (2018b) "Cope with diverse data structures in multi-fidelity modeling: a gaussian process method". Engineering Applications of Artificial Intelligence 67:211–225.
[^9]: Wu Y, Hu J, Zhou Q, Wang S, Jin P (2020) "An active learning multi-fidelity metamodeling method based on the bootstrap estimator". Aerospace Science and Technology 106:106116.
[^10]: Shi M, Lv L, Sun W, Song X (2020) "A multi-fidelity surrogate model based on support vector regression". Structural and Multidisciplinary Optimization 1–13.
[^11]: Dong H, Song B, Wang P, Huang S (2015) "Multi-fidelity information fusion based on prediction of kriging". Structural and Multidisciplinary Optimization 51(6):1267–1280.
[^12]: Xiong S, Qian PZ, Wu CJ (2013) "Sequential design and analysis of high-accuracy and low-accuracy computer codes". Technometrics 55(1):37–46.
[^13]: Park C, Haftka RT, Kim NH (2018) "Low-fidelity scale factor improves bayesian multi-fidelity prediction by reducing bumpiness of discrepancy function". Structural and Multidisciplinary Optimization 58(2):399–414.
[^14]: Hansen N, Auger A, Ros R, Mersmann O, Tuˇsar T, Brockhoff D (2020) "COCO: A platform for comparing continuous optimizers in a black-box setting". Optimization Methods and Software URL http://dx.doi.org/https://doi.org/10.1080/10556788.2020.1808977
[^15]: Kerschke, P. & Trautmann, H. (2019). Comprehensive Feature-Based Landscape Analysis of Continuous and Constrained Optimization Problems Using the R-package flacco. In: Bauer N., Ickstadt K., Lübke K., Szepannek G., Trautmann H., Vichi M. (eds.) Applications in Statistical Computing -- From Music Data Analysis to Industrial Quality Improvement, pp. 93-123, Studies in Classification, Data Analysis, and Knowledge Organization, Springer. URL: https://link.springer.com/chapter/10.1007/978-3-030-25147-5_7