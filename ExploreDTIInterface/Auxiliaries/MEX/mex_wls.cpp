// Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
// under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

/*
 C++ mex implementation of weighted least squares fit

 Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
            University of Antwerp
*/

#include <mex.h>
#include <lapack.h>

inline void gels(char* trans, ptrdiff_t* m, ptrdiff_t* n, ptrdiff_t* nrhs, float* a, ptrdiff_t* lda, float* b, ptrdiff_t* ldb, float* work, ptrdiff_t* lwork, ptrdiff_t* info) {
    sgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
};

inline void gels(char* trans, ptrdiff_t* m, ptrdiff_t* n, ptrdiff_t* nrhs, double* a, ptrdiff_t* lda, double* b, ptrdiff_t* ldb, double* work, ptrdiff_t* lwork, ptrdiff_t* info) {
    dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
};

template <typename T> inline void wls(T* A, ptrdiff_t A_r, ptrdiff_t A_c, T* b, ptrdiff_t b_r, ptrdiff_t b_c, T* w, T* x) {
    /* allocation */
    T* A_i = (T*) mxMalloc(A_r*A_c*sizeof(T));
    T* b_i = (T*) mxMalloc(b_r*sizeof(T));
    T* work = (T*) mxMalloc(sizeof(T));
    
    char trans = 'N'; ptrdiff_t info; ptrdiff_t one = 1;
    
    /* worksize query */
    ptrdiff_t lwork = -1;
    gels(&trans, &A_r, &A_c, &one, A_i, &A_r, b_i, &b_r, work, &lwork, &info);
    lwork = work[0]; mxFree(work); work = (T*) mxMalloc(lwork*sizeof(T));
    
    /* actual calculation */
    for (ptrdiff_t i = 0; i < b_c; i++) {
        /* multiply A_i with w_i */
        for (ptrdiff_t j = 0; j < A_c; j++) { for (ptrdiff_t k = 0; k < A_r; k++) { A_i[j*A_r+k] = A[j*A_r+k]*w[i*b_r+k]; } }
        /* multiply b_i with w_i */
        for (ptrdiff_t j = 0; j < b_r; j++) { b_i[j] = b[i*b_r+j]*w[i*b_r+j]; }
        /* calculate least squares solution x = Aw\bw */
        gels(&trans, &A_r, &A_c, &one, A_i, &A_r, b_i, &b_r, work, &lwork, &info);
        /* copy solution x_i from b_i */
        for (ptrdiff_t j = 0; j < A_c; j++) { x[i*A_c+j] = b_i[j]; }
    }
    
    /* deallocation */
    mxFree(A_i);
    mxFree(b_i);
    mxFree(work);
}

template <typename T> inline void wls_function(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    /* input */
	T* A = (T*) mxGetData(prhs[0]); ptrdiff_t A_r = mxGetM(prhs[0]); ptrdiff_t A_c = mxGetN(prhs[0]);
    T* b = (T*) mxGetData(prhs[1]); ptrdiff_t b_r = mxGetM(prhs[1]); ptrdiff_t b_c = mxGetN(prhs[1]);
    T* w = (T*) mxGetData(prhs[2]); ptrdiff_t w_r = mxGetM(prhs[2]); ptrdiff_t w_c = mxGetN(prhs[2]);
    
    /* output */
    plhs[0] = mxCreateNumericMatrix(A_c, b_c, mxGetClassID(prhs[0]), mxREAL);
    T* x = (T*) mxGetData(plhs[0]); ptrdiff_t x_r = mxGetM(plhs[0]); ptrdiff_t x_c = mxGetN(plhs[0]);
    
    /* calculation */
    wls(A, A_r, A_c, b, b_r, b_c, w, x);
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    /* check parameters */
    if (nrhs != 3) {
        mexErrMsgTxt("Usage: x = mex_wls(A,b,w)");
    }

    for (int i = 1; i < 3; i++) {
        if (mxGetClassID(prhs[0]) != mxGetClassID(prhs[i])) {
            mexErrMsgTxt("Types do not match!");
        }
        if (mxGetNumberOfDimensions(prhs[i]) != 2) {
            mexErrMsgTxt("Wrong number of dimensions!");
        }
        if (mxGetM(prhs[0]) != mxGetM(prhs[i])) {
            mexErrMsgTxt("Dimensions don't match!");
        }
    }
    if (mxGetN(prhs[1]) != mxGetN(prhs[2])) {
        mexErrMsgTxt("Dimensions don't match!");
    }
    
    /* run it for the right type */
    switch (mxGetClassID(prhs[0])) {
        case mxSINGLE_CLASS:
            wls_function<float>(nlhs, plhs, nrhs, prhs);
            break;
        case mxDOUBLE_CLASS:
            wls_function<double>(nlhs, plhs, nrhs, prhs);
            break;
        default:
            mexErrMsgTxt("Type not supported yet!");
            break;
    }
}