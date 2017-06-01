// Matrix.h

#ifndef IMEX_MATRIX_H
#define IMEX_MATRIX_H

#include "mex.h"

#include <complex>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace imex {

	template <typename T> class Matrix : public Eigen::Matrix<T, -1, -1> {

	public:

		// input routines
		static Eigen::Matrix<T, -1, -1> InputWrapper(const mxArray* in);
		static void InputWrapper(Eigen::Matrix<T, -1, -1>& mat, const mxArray* in);

		//output routines
		static mxArray* OutputWrapper(const Eigen::Matrix<T, -1, -1>& mat);

	};

	template <typename T> class Matrix< std::complex<T> > : public Eigen::Matrix<std::complex<T>, -1, -1> {

	public:

		// input routines
		static Eigen::Matrix< std::complex<T>, -1, -1> InputWrapper(const mxArray* in);
		static void InputWrapper(Eigen::Matrix<std::complex<T>, -1, -1>& mat, const mxArray* in);
		
		//output routines
		static mxArray* OutputWrapper(const Eigen::Matrix<std::complex<T>, -1, -1>& mat);

	};

}


// input routines
template <typename T> 
Eigen::Matrix<T, -1, -1> imex::Matrix<T>::InputWrapper(const mxArray* in) {

	validateNumericInput<T>(in, -1, -1);

	Eigen::Matrix<T, -1, -1> mat;

	if (mxIsSparse(in))
		mat = Eigen::Matrix<T, -1, -1>(imex::SparseMatrix<T>::InputWrapper(in));
	else
		imex::Matrix<T>::InputWrapper(mat, in);

	return mat;

}

template <typename T> 
Eigen::Matrix<std::complex<T>, -1, -1> imex::Matrix< std::complex<T> >::InputWrapper(const mxArray* in) {

	validateNumericInput< std::complex<T> >(in, -1, -1);

	Eigen::Matrix<std::complex<T>, -1, -1> mat;

	if (mxIsSparse(in))
		mat = Eigen::Matrix<std::complex<T>, -1, -1>(imex::SparseMatrix< std::complex<T> >::InputWrapper(in));
	else
		imex::Matrix< std::complex<T> >::InputWrapper(mat, in);

	return mat;

}

template <typename T> 
void imex::Matrix<T>::InputWrapper(Eigen::Matrix<T, -1, -1>& mat, const mxArray* in) {

	if (imex::Type<T>() == mxGetClassID(in))
		mat = Eigen::Map<Eigen::Matrix<T, -1, -1> >((T*) mxGetData(in), mxGetM(in), mxGetN(in));
	else
		mat = Eigen::Map<Eigen::Matrix<T, -1, -1> >(imex::Cast<T>(in), mxGetM(in), mxGetN(in));

}

template <typename T> 
void imex::Matrix< std::complex<T> >::InputWrapper(Eigen::Matrix<std::complex<T>, -1, -1>& mat, const mxArray* in) {

	int M = mxGetM(in), N = mxGetN(in);

	mat = Eigen::Matrix<std::complex<T>, -1, -1>(M, N);

	if (imex::Type<T>() == mxGetClassID(in)) {

		mat.real() = Eigen::Map<Eigen::Matrix<T, -1, -1> >((T*) mxGetData(in), M, N);

		if (mxIsComplex(in))
			mat.imag() = Eigen::Map<Eigen::Matrix<T, -1, -1> >((T*) mxGetImagData(in), M, N);
		else
			mat.imag() = Eigen::Matrix<T, -1, -1>::Zero(M, N);

	} else {

		mat.real() = Eigen::Map<Eigen::Matrix<T, -1, -1> >(imex::Cast<T>(in), M, N);

		if (mxIsComplex(in))
			mat.imag() = Eigen::Map<Eigen::Matrix<T, -1, -1> >(imex::CastImag<T>(in), M, N);
		else
			mat.imag() = Eigen::Matrix<T, -1, -1>::Zero(M, N);

	}

}

//output routines
template <typename T> 
mxArray* imex::Matrix<T>::OutputWrapper(const Eigen::Matrix<T, -1, -1>& mat) {

	mxArray* out = mxCreateNumericMatrix(mat.rows(), mat.cols(), imex::Type<T>(), mxREAL);
	memcpy((T*) mxGetData(out), mat.data(), mat.size()*sizeof(T));

	return out;

}

template <typename T> 
mxArray* imex::Matrix< std::complex<T> >::OutputWrapper(const Eigen::Matrix<std::complex<T>, -1, -1>& mat) {

	Eigen::Matrix<T, -1, -1> matR = mat.real(), matI = mat.imag();

	mxArray* out = mxCreateNumericMatrix(mat.rows(), mat.cols(), imex::Type<T>(), mxCOMPLEX);
	memcpy((T*) mxGetData(out), matR.data(), matR.size()*sizeof(T));
	memcpy((T*) mxGetImagData(out), matI.data(), matI.size()*sizeof(T));

	return out;

}

#endif