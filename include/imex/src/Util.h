// Util.h

#ifndef IMEX_UTIL_H
#define IMEX_UTIL_H

#include "mex.h"

#include <type_traits>
#include <complex>
#include <exception>

namespace imex {
	
	template <typename T> void validateNumericInput(const mxArray* in, int M, int N) {
		
		if (!mxIsNumeric(in)) 
			throw(std::invalid_argument("Input must be numeric"));
		
		if (mxIsComplex(in) && imex::Complexity<T>() == mxREAL)
			throw(std::invalid_argument("Complex input assigned to non-complex object"));
		
		validateArrayInput<T>(in, M, N);
		
	}
	
	template <typename T> void validateLogicalInput(const mxArray* in, int M, int N) {
		
		if (!mxIsLogical(in)) 
			throw(std::invalid_argument("Input must be logical"));
		
		validateArrayInput<T>(in, M, N);
		
	}
	
	template <typename T> void validateLogicalOrNumericInput(const mxArray* in, int M, int N) {
		
		if (!mxIsLogical(in) && !mxIsNumeric(in)) 
			throw(std::invalid_argument("Input must be logical or numeric"));
		
		if (mxIsComplex(in) && imex::Complexity<T>() == mxREAL)
			throw(std::invalid_argument("Complex input assigned to non-complex object"));
		
		validateArrayInput<T>(in, M, N);
		
	}
	
	template <typename T> void validateArrayInput(const mxArray* in, int M, int N) {
		
		if (mxGetNumberOfDimensions(in) > 2)
			throw(std::invalid_argument("Multidimensional arrays are not supported"));
		
		if (M > 0 && mxGetM(in) != M) 
			throw(std::invalid_argument("Bad number of rows"));
		
		if (N > 0 && mxGetN(in) != N) 
			throw(std::invalid_argument("Bad number of columns"));
		
	}
	
}


#endif