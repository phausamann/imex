# imex
A MATLAB MEX to Eigen interface in C++

## Overview
A header-only library written in C++ that provides various functions for easy conversions from `mxArray` pointers to Eigen objects and back. Inputs will be typecast automatically, if possible.

## Running the example
The example demonstrates how inputs are automatically typecast depending on the requested Eigen type.
* Set the path to your Eigen installation in line 21 of `gen_mex.m`
* Run `test_imex.m`
* `C` should be equal to `A` (except for the datatype)

## Disclaimer
This is a very basic library in its early alpha stage. A lot of features are yet to be implemented (e.g. `Eigen::Array`, `Eigen::SparseVector`, support of MATLAB logicals).

So far, the library has only been tested on MATLAB R2016a x64 on Windows 8/10 with the Microsoft Visual C++ 2010 Compiler and Eigen 3.2.1.
