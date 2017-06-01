// Vector.h

#ifndef IMEX_VECTOR_H
#define IMEX_VECTOR_H

#include "mex.h"

#include <vector>
#include <complex>

#include <Eigen/Core>

namespace imex {
	
	template <typename T> class Vector : public Eigen::Matrix<T, -1, 1> {

	public:

		// input routines
		static Eigen::Matrix<T, -1, 1> InputWrapper(const mxArray* in);

		//output routines
		static mxArray* OutputWrapper(const Eigen::Matrix<T, -1, 1>& vec);
		static mxArray* OutputWrapper(const std::vector<T>& vec);

	};
	
	template <typename T> class Vector< std::complex<T> > : public Eigen::Matrix<std::complex<T>, -1, 1> {

	public:

		// input routines
		static Eigen::Matrix<std::complex<T>, -1, 1> InputWrapper(const mxArray* in);

		//output routines
		static mxArray* OutputWrapper(const Eigen::Matrix<std::complex<T>, -1, 1>& vec);
		static mxArray* OutputWrapper(const std::vector< std::complex<T> >& vec);

	};

}


// input routines
template <typename T> 
Eigen::Matrix<T, -1, 1> imex::Vector<T>::InputWrapper(const mxArray* in) {
	
	validateNumericInput<T>(in, -1, 1);
	
	Eigen::Matrix<T, -1, -1> mat;
	
	imex::Matrix<T>::InputWrapper(mat, in);
	
	return mat.col(0);
	
}

template <typename T> 
Eigen::Matrix<std::complex<T>, -1, 1> imex::Vector< std::complex<T> >::InputWrapper(const mxArray* in) {
	
	validateNumericInput< std::complex<T> >(in, -1, 1);
	
	Eigen::Matrix<std::complex<T>, -1, -1> mat;
	
	imex::Matrix< std::complex<T> >::InputWrapper(mat, in);
	
	return mat.col(0);
	
}

//output routines
template <typename T> 
mxArray* imex::Vector<T>::OutputWrapper(const Eigen::Matrix<T, -1, 1>& vec) {
	
	Eigen::Matrix<T, -1, -1> mat(vec);
	
	return imex::Matrix<T>::OutputWrapper(mat);
	
}

template <typename T> 
mxArray* imex::Vector< std::complex<T> >::OutputWrapper(const Eigen::Matrix<std::complex<T>, -1, 1>& vec) {
	
	Eigen::Matrix<std::complex<T>, -1, -1> mat(vec);
	
	return imex::Matrix< std::complex<T> >::OutputWrapper(mat);
	
}

template <typename T> 
mxArray* imex::Vector<T>::OutputWrapper(const std::vector<T>& vec) {
	
	mxArray* out = mxCreateNumericMatrix(vec.size(), 1, imex::Type<T>(), mxREAL);
	memcpy((T*) mxGetData(out), vec.data(), vec.size()*sizeof(T));
	
	return out;
	
}

template <typename T> 
mxArray* imex::Vector< std::complex<T> >::OutputWrapper(const std::vector< std::complex<T> >& vec) {
	
	std::vector<T> vecR(vec.size(), T()), vecI(vec.size(), T());
	for (int i=0; i<vec.size(); i++) {
		vecR[i] = vec[i].real();
		vecI[i] = vec[i].imag();
	}
	
	mxArray* out = mxCreateNumericMatrix(vec.size(), 1, imex::Type<T>(), mxCOMPLEX);
	memcpy((T*) mxGetData(out), vecR.data(), vecR.size()*sizeof(T));
	memcpy((T*) mxGetImagData(out), vecI.data(), vecI.size()*sizeof(T));
	
	return out;
	
}

#endif