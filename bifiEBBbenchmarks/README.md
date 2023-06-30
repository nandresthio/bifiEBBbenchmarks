# Bi-fidelity Surrogate Modelling benchmark problems

The software and data in this repository aim to provide a benchmark for Bi-Fidelity Expensive Black Box problems, both for surrogate model building and for optimisation. It implements both classical instances from the literature, as well as the instance generating procedure presented in the paper [Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances](https://doi.org/10.1287/ijoc.2019.0934) by Andrés-Thió N, Muñoz MA, and Smith-Miles K. It also employs instance measures in order to generate varied testing suites of different sizes, to be used by the community when analysing the performance of new Bi-fidelity Expensive Black-Box problems.

The main benefits of using this piece of code are the following:

- A very large set of candidate functions (over 200 literature instances, and a virtually endless set of newly proposed instances).
- Latin Hypercube Sampling (LHS) of varying size and dimension which are well spread out are supplied for reproducibility purposes.
- A set of instance measures, known as features, are supplied for all instances, which can be used to analyse the difference between them.
- Instance subsets of different sizes are supplied, chosen to be as varied as possible using an instance filtering procedure.

Note that no surrogate models or algorithms are supplied in this reposoritory, although a different repository which provides an implementation of Kriging and Co-Kriging among others and is run on this set of benchmarks can be found here (ADD LINK).

## Cite

To cite this software, please cite the [paper](https://doi.org/10.1287/ijoc.2019.0934) (LINK IS INCORRECT, UPDATE TO WHAT IS TO BE SUBMITTED) using its DOI and the software itself, using the following DOI. (AGAIN, LINK IS NOT CORRECT, UPDATE)

[![DOI](https://zenodo.org/badge/285853815.svg)](https://zenodo.org/badge/latestdoi/285853815)

Below is the BibTex for citing this version of the code. (UPDATE ZENODO)

```
@misc{nandres2022,
  author = {Andr\'{e}s-Thi\'{o}, Nicolau},
  title = {Bifidelity Surrogate Modelling benchmark problems},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  note = {available for download at https://github.com/nandresthio/bifiEBBbenchmarks},
  doi = {10.5281/zenodo.6208147}
} 
```


## Description

This code implements test functions for Bi-Fidelity Expensive Black-Box problems. These problems have a high and low source of information, which are used to either accurately model or to optimise the high fidelity source. This repository provides implementations of a vast set of functions pairs, provides Latin Hypercube Sampling (LHS) plans which have been spread out.

TO BE CONTINUED!!

## To do
 - Add instructions on populating folders `data/samplePlans` and `features` from figshare
 - Add instructions on adding SOLAR folder if want to use SOLAR instances
 - Deal with linking Eigen when using Sourcecpp

