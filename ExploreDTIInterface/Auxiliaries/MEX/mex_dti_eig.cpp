// Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
// under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
/*
 C++ mex implementation of diffusion tensor eigenvalues
 Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
            University of Antwerp
*/

#include <mex.h>
#include <lapack.h>

inline void spevd(char* jobz, char* uplo, ptrdiff_t* n, float* ap, float* w, float* z, ptrdiff_t* ldz, float* work, ptrdiff_t *lwork, ptrdiff_t* iwork, ptrdiff_t* liwork, ptrdiff_t* info) {
    sspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
}

inline void spevd(char* jobz, char* uplo, ptrdiff_t* n, double* ap, double* w, double* z, ptrdiff_t* ldz, double* work, ptrdiff_t *lwork, ptrdiff_t* iwork, ptrdiff_t* liwork, ptrdiff_t* info) {
    dspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
}

template<class T> inline void eig_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    T* dt = (T*) mxGetData(prhs[0]);

    plhs[0] = mxCreateNumericMatrix(3,mxGetN(prhs[0]),mxGetClassID(prhs[0]),mxREAL);
    plhs[1] = mxCreateNumericMatrix(9,mxGetN(prhs[0]),mxGetClassID(prhs[0]),mxREAL);
    T* eigval = (T*) mxGetData(plhs[0]);
    T* eigvect = (T*) mxGetData(plhs[1]);
    
    char jobz = 'V'; char uplo = 'L';
    ptrdiff_t n = 3; ptrdiff_t ldz = 3;
    T* work = new T[28]; ptrdiff_t lwork = 28;
    ptrdiff_t* iwork = new ptrdiff_t[18]; ptrdiff_t liwork = 18;
    ptrdiff_t info = 0;
    for (int pos = 0; pos < mxGetN(prhs[0]); pos++) {
        spevd(&jobz, &uplo, &n, &dt[pos*6], &eigval[pos*3], &eigvect[pos*9], &ldz, work, &lwork, iwork, &liwork, &info);
    }
    delete[] work;
    delete[] iwork;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    /* check parameters */
    if (nrhs != 1 || nlhs != 2) {
        mexErrMsgTxt("Usage: [eigval, eigvect] = mex_dti_eig(dt)");
    }

    if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetM(prhs[0]) != 6) {
        mexErrMsgTxt("dt must be a 6 x n single or double matrix ");
    }
    
    /* run it for the right type */
    switch (mxGetClassID(prhs[0])) {
        case mxSINGLE_CLASS:
            eig_function<float>(nlhs, plhs, nrhs, prhs);
            break;
        case mxDOUBLE_CLASS:
            eig_function<double>(nlhs, plhs, nrhs, prhs);
            break;
        default:
            mexErrMsgTxt("Type not supported!");
            break;
    }
}
