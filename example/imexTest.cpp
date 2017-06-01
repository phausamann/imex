// imexTest.cpp

#include "mex.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <imex/Core>
#include <imex/SparseCore>

#define A_IN prhs[0]
#define B_IN prhs[1]

#define C_OUT plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    typedef double T;
    typedef std::complex<double> Tc;
    
    if (mxIsComplex(A_IN) || mxIsComplex(B_IN)) {
        
        Eigen::Matrix<Tc, -1, -1> A = imex::Matrix<Tc>::InputWrapper(A_IN);
        
        if (nrhs == 1) {
            C_OUT = imex::Matrix<Tc>::OutputWrapper(A);
        } else {
            Eigen::SparseMatrix<Tc> B = imex::SparseMatrix<Tc>::InputWrapper(B_IN);
            Eigen::Matrix<Tc, -1, -1> C = A*B;
            C_OUT = imex::Matrix<Tc>::OutputWrapper(C);
        }
        
    } else {
        
        Eigen::Matrix<T, -1, -1> A = imex::Matrix<T>::InputWrapper(A_IN);
        
        if (nrhs == 1) {
            C_OUT = imex::Matrix<T>::OutputWrapper(A);
        } else {
            Eigen::SparseMatrix<T> B = imex::SparseMatrix<T>::InputWrapper(B_IN);
            Eigen::Matrix<T, -1, -1> C = A*B;
            C_OUT = imex::Matrix<T>::OutputWrapper(C);
        }
                
    }
    
}