// Type.h

#ifndef IMEX_TYPE_H
#define IMEX_TYPE_H

#include "mex.h"

#include <type_traits>
#include <complex>

namespace imex {

	template <typename T> T* Cast(const mxArray* in, const bool& imag = false) {
		
		switch (mxGetClassID(in)) {

		case mxDOUBLE_CLASS:
			 return imex::Cast<T, double>(in, imag);

		case mxSINGLE_CLASS:
			 return imex::Cast<T, float>(in, imag);

		case mxINT64_CLASS:
			return imex::Cast<T, long long>(in, imag);

		case mxINT32_CLASS:
			return imex::Cast<T, int>(in, imag);

		case mxINT16_CLASS:
			return imex::Cast<T, short>(in, imag);

		case mxINT8_CLASS:
			return imex::Cast<T, char>(in, imag);

		case mxUINT64_CLASS:
			return imex::Cast<T, unsigned long long>(in, imag);

		case mxUINT32_CLASS:
			return imex::Cast<T, unsigned int>(in, imag);

		case mxUINT16_CLASS:
			return imex::Cast<T, unsigned short>(in, imag);

		case mxUINT8_CLASS:
			return imex::Cast<T, unsigned char>(in, imag);

		default:
			throw(std::invalid_argument("Cannot cast non-numeric array types"));
		}

		return nullptr;

	}

	template <typename T> T* CastImag(const mxArray* in) {
		return imex::Cast<T>(in, true);
	}

	template <typename T, typename S> T* Cast(const mxArray* in, const bool& imag = false) {

		int nElements = mxGetNumberOfElements(in);
		T* T_cptr = new T[nElements];
		S* S_cptr;

		if (imag)
			S_cptr = (S*) mxGetImagData(in);
		else
			S_cptr = (S*) mxGetData(in);

		for (int it=0; it<nElements; it++)
			T_cptr[it] = static_cast<T>(*S_cptr++);

		return T_cptr;
	}
	
	template <typename T> mxClassID Type(void) {
		return mxUNKNOWN_CLASS;
	}
	
	template <> mxClassID Type<bool>(void) {
		return mxLOGICAL_CLASS;
	}
	
	template <> mxClassID Type<double>(void) {
		return mxDOUBLE_CLASS;
	}
	
	template <> mxClassID Type<float>(void) {
		return mxSINGLE_CLASS;
	}
	
	template <> mxClassID Type<long long>(void) {
		return mxINT64_CLASS;
	}
	
	template <> mxClassID Type<int>(void) {
		return mxINT32_CLASS;
	}
	
	template <> mxClassID Type<short>(void) {
		return mxINT16_CLASS;
	}
	
	template <> mxClassID Type<char>(void) {
		return mxINT8_CLASS;
	}
	
	template <> mxClassID Type<unsigned long long>(void) {
		return mxUINT64_CLASS;
	}
	
	template <> mxClassID Type<unsigned int>(void) {
		return mxUINT32_CLASS;
	}
	
	template <> mxClassID Type<unsigned short>(void) {
		return mxUINT16_CLASS;
	}
	
	template <> mxClassID Type<unsigned char>(void) {
		return mxUINT8_CLASS;
	}
	
	template<typename T> struct is_complex : std::false_type {};
	
	template<typename T> struct is_complex< std::complex<T> > : std::true_type {};
	
	template <typename T> mxComplexity Complexity(typename std::enable_if< is_complex<T>::value >::type* dummy = 0) {
		return mxCOMPLEX;
	}
	
	template <typename T> mxComplexity Complexity(typename std::enable_if< !is_complex<T>::value >::type* dummy = 0) {
		return mxREAL;
	}
	
}


#endif