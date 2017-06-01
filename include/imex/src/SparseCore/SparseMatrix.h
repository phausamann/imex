// SparseMatrix.h

#ifndef IMEX_SPARSEMATRIX_H
#define IMEX_SPARSEMATRIX_H

#include "mex.h"

#include <vector>
#include <complex>
#include <exception>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace imex {
	
	template <typename T> class SparseMatrix : public Eigen::SparseMatrix<T> {

	public:

		// input routines
		static Eigen::SparseMatrix<T> InputWrapper(const mxArray* in);
		static void InputWrapper(Eigen::SparseMatrix<T>& mat, const mxArray* in);

		//output routines
		static mxArray* OutputWrapper(Eigen::SparseMatrix<T>& mat);
		static mxArray* OutputWrapper(int M, int N, const std::vector<T>& pr, const std::vector<int>& ir, const std::vector<int>& jc);

	};
	
	template <typename T> class SparseMatrix< std::complex<T> > : public Eigen::SparseMatrix< std::complex<T> > {

	public:

		// input routines
		static Eigen::SparseMatrix< std::complex<T> > InputWrapper(const mxArray* in);
		static void InputWrapper(Eigen::SparseMatrix< std::complex<T> >& mat, const mxArray* in);

		//output routines
		static mxArray* OutputWrapper(Eigen::SparseMatrix< std::complex<T> >& mat);
		static mxArray* OutputWrapper(int M, int N, const std::vector< std::complex<T> >& p, const std::vector<int>& ir, const std::vector<int>& jc);

	};

}


// input routines
template <typename T> 
Eigen::SparseMatrix<T> imex::SparseMatrix<T>::InputWrapper(const mxArray* in) {
	
	validateNumericInput<T>(in, -1, -1);
	
	Eigen::SparseMatrix<T> mat;
	
	if (mxIsSparse(in)) {
	
		imex::SparseMatrix<T>::InputWrapper(mat, in);
		
	} else {
		
		Eigen::Matrix<T, -1, -1> dmat = imex::Matrix<T>::InputWrapper(in);
		mat = dmat.sparseView();
		
	}
	
	return mat;
	
}

template <typename T> 
Eigen::SparseMatrix< std::complex<T> > imex::SparseMatrix< std::complex<T> >::InputWrapper(const mxArray* in) {
	
	validateNumericInput< std::complex<T> >(in, -1, -1);
	
	Eigen::SparseMatrix< std::complex<T> > mat;
	
	if (mxIsSparse(in)) {
	
		imex::SparseMatrix< std::complex<T> >::InputWrapper(mat, in);
		
	} else {
		
		Eigen::Matrix<std::complex<T>, -1, -1> dmat = imex::Matrix< std::complex<T> >::InputWrapper(in);
		mat = dmat.sparseView();
		
	}
	
	return mat;
	
}

template <typename T> 
void imex::SparseMatrix<T>::InputWrapper(Eigen::SparseMatrix<T>& mat, const mxArray* in) {
		
	mat = Eigen::SparseMatrix<T>(mxGetM(in), mxGetN(in));
	
	T* pr_cptr;
	mwIndex* ir_cptr = mxGetIr(in);
	mwIndex* jc_cptr = mxGetJc(in);

	int nnz = mxGetNzmax(in);
	std::vector< Eigen::Triplet<T> > tripletList;
	tripletList.reserve(nnz);

	if (imex::Type<T>() == mxGetClassID(in))
		pr_cptr = (T*) mxGetData(in);
	else
		pr_cptr = imex::Cast<T>(in);
	
	int j=0;
	for (int k=0; k<nnz; k++) {
		while (jc_cptr[j+1] <= k && j <= mxGetN(in)) j++;
		if (jc_cptr[j+1] == 0) break; // necessary for empty sparse matrices apparently
		tripletList.push_back(Eigen::Triplet<T>(ir_cptr[k], j, pr_cptr[k]));
	}
	
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
}

template <typename T> 
void imex::SparseMatrix< std::complex<T> >::InputWrapper(Eigen::SparseMatrix< std::complex<T> >& mat, const mxArray* in) {

	T *pr_cptr, *pi_cptr;
		
	mat = Eigen::SparseMatrix< std::complex<T> >(mxGetM(in), mxGetN(in));

	mwIndex* ir_cptr = mxGetIr(in);
	mwIndex* jc_cptr = mxGetJc(in);

	int nnz = mxGetNzmax(in);
	std::vector< Eigen::Triplet< std::complex<T> > > tripletList;
	tripletList.reserve(nnz);

	if (imex::Type<T>() == mxGetClassID(in)) {
	
		pr_cptr = (T*) mxGetData(in);

		if (mxIsComplex(in))
			pi_cptr = (T*) mxGetImagData(in);
		else
			pi_cptr = new T[mxGetNumberOfElements(in)]();

	} else {
	
		pr_cptr = imex::Cast<T>(in);

		if (mxIsComplex(in))
			pi_cptr = imex::CastImag<T>(in);
		else
			pi_cptr = new T[mxGetNumberOfElements(in)]();

	}
	
	int j=0;
	for (int k=0; k<nnz; k++ && j <= mxGetN(in)) {
		while (jc_cptr[j+1] <= k) j++;
		if (jc_cptr[j+1] == 0) break; // necessary for empty sparse matrices apparently
		tripletList.push_back(Eigen::Triplet< std::complex<T> >(ir_cptr[k], j, std::complex<T>(pr_cptr[k], pi_cptr[k])));
	}
	
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
}

//output routines
template <typename T> 
mxArray* imex::SparseMatrix<T>::OutputWrapper(Eigen::SparseMatrix<T>& mat) {
	
	// TODO: make faster (memcpy for mat.valuePtr()?)
	
	std::vector<T> p;
	std::vector<int> ir, jc(mat.cols()+1, 0);
	
	p.reserve(mat.nonZeros());
	ir.reserve(mat.nonZeros());
	
	for (int k=0; k<mat.outerSize(); k++) {
		jc[k+1] = jc[k];
		for (Eigen::SparseMatrix<T>::InnerIterator it(mat,k); it; ++it) {
			p.push_back(it.value());
			ir.push_back(it.row());
			jc[it.col()+1]++;
		 }
	}
	
	return imex::SparseMatrix<T>::OutputWrapper(mat.rows(), mat.cols(), p, ir, jc);
	
}

template <typename T> 
mxArray* imex::SparseMatrix< std::complex<T> >::OutputWrapper(Eigen::SparseMatrix< std::complex<T> >& mat) {
	
	// TODO: make faster (memcpy for mat.valuePtr()?)
	
	std::vector<T> p;
	std::vector<int> ir, jc(mat.cols()+1, 0);
	
	p.reserve(mat.nonZeros());
	ir.reserve(mat.nonZeros());
	
	for (int k=0; k<mat.outerSize(); k++) {
		jc[k+1] = jc[k];
		for (Eigen::SparseMatrix< std::complex<T> >::InnerIterator it(mat,k); it; ++it) {
			p.push_back(it.value());
			ir.push_back(it.row());
			jc[it.col()+1]++;
		 }
	}
	
	return imex::SparseMatrix< std::complex<T> >::OutputWrapper(mat.rows(), mat.cols(), p, ir, jc);
	
}

template <typename T> 
mxArray* imex::SparseMatrix<T>::OutputWrapper(int M, int N, const std::vector<T>& pr, const std::vector<int>& ir, const std::vector<int>& jc) {
	
	mxArray* out = mxCreateSparse(M, N, pr.size(), mxREAL);
	
	mwIndex* ir_cptr = mxGetIr(out);
	mwIndex* jc_cptr = mxGetJc(out);
	
	memcpy((T*) mxGetData(out), pr.data(), pr.size()*sizeof(T));
	
	for (int it=0; it<ir.size(); it++) ir_cptr[it] = static_cast<mwIndex>(ir[it]);
	for (int it=0; it<jc.size(); it++) jc_cptr[it] = static_cast<mwIndex>(jc[it]);
	
	return out;
	
}

template <typename T> 
mxArray* imex::SparseMatrix< std::complex<T> >::OutputWrapper(int M, int N, const std::vector< std::complex<T> >& p, const std::vector<int>& ir, const std::vector<int>& jc) {
	
	std::vector<T> pr(p.size(), T()), pi(p.size(), T());
	for (int i=0; i<p.size(); i++) {
		pr[i] = p[i].real();
		pi[i] = p[i].imag();
	}
	
	mxArray* out = mxCreateSparse(M, N, p.size(), mxCOMPLEX);

	mwIndex* ir_cptr = mxGetIr(out);
	mwIndex* jc_cptr = mxGetJc(out);
	
	memcpy((T*) mxGetData(out), pr.data(), pr.size()*sizeof(T));
	memcpy((T*) mxGetImagData(out), pi.data(), pi.size()*sizeof(T));
	
	for (int it=0; it<ir.size(); it++) ir_cptr[it] = static_cast<mwIndex>(ir[it]);
	for (int it=0; it<jc.size(); it++) jc_cptr[it] = static_cast<mwIndex>(jc[it]);
	
	return out;
	
}

#endif