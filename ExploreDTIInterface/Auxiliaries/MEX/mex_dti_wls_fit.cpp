#include <mex.h>
#include <cmath>
#include <lapack.h>
// Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
// under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

inline void gels(char* trans, ptrdiff_t* m, ptrdiff_t* n, ptrdiff_t* nrhs, float* a, ptrdiff_t* lda, float* b, ptrdiff_t* ldb, float* work, ptrdiff_t* lwork, ptrdiff_t* info) {
    sgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
};

inline void gels(char* trans, ptrdiff_t* m, ptrdiff_t* n, ptrdiff_t* nrhs, double* a, ptrdiff_t* lda, double* b, ptrdiff_t* ldb, double* work, ptrdiff_t* lwork, ptrdiff_t* info) {
    dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
};

template<class T> inline void pinv(T* A, T* AI, ptrdiff_t m, ptrdiff_t n) {
    ptrdiff_t ldb = m < n ? n : m;
    T* B = (T *)mxMalloc(ldb*m*sizeof(T));
	for(int i = 0; i<ldb; i++){
        for(int j = 0; j < m; j++){
            if(i == j) B[i+j*ldb] = 1.0;
            else       B[i+j*ldb] = 0.0;
        }
	}
    char trans = 'N'; ptrdiff_t lwork = -1; ptrdiff_t info;
    T* work = (T*) mxMalloc(sizeof(T));
    gels(&trans, &m, &n, &m, A, &m, B, &ldb, work, &lwork, &info);
    lwork = work[0]; mxFree(work); work = (T*) mxMalloc(lwork*sizeof(T));
    gels(&trans, &m, &n, &m, A, &m, B, &ldb, work, &lwork, &info);
    mxFree(work);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            AI[i+j*n] = B[i+j*ldb];
        }
    }
    mxFree(B);
}

template<class T> inline void tensorfit(T dt[7], int nVols, const T* w, const T* logData, const T* Xinv) {
    for(int i = 0; i < 7; i++) { dt[i] = 0; }
	for(int i = 0; i < nVols; i++) {
        T d;
        if (w == NULL) {
            d = logData[i];
        } else {
            d = w[i]*logData[i];
        }
        int j = i*7;
        dt[0] = dt[0] + d*Xinv[j  ];
        dt[1] = dt[1] + d*Xinv[j+1];
        dt[2] = dt[2] + d*Xinv[j+2];
        dt[3] = dt[3] + d*Xinv[j+3];
        dt[4] = dt[4] + d*Xinv[j+4];
        dt[5] = dt[5] + d*Xinv[j+5];
        dt[6] = dt[6] + d*Xinv[j+6];
    }
}

template<class T> inline void computeWeightDiag(T* W, T* wX, int nVols, const T* X, const T* dt) {
	for(int j = 0; j < nVols; j++){
        T pLogData = 0;
        int l = j ; pLogData = pLogData + X[l] * dt[0];
        l += nVols; pLogData = pLogData + X[l] * dt[1];
        l += nVols; pLogData = pLogData + X[l] * dt[2];
        l += nVols; pLogData = pLogData + X[l] * dt[3];
        l += nVols; pLogData = pLogData + X[l] * dt[4];
        l += nVols; pLogData = pLogData + X[l] * dt[5];
        l += nVols; pLogData = pLogData + X[l] * dt[6];
        W[j] = exp(pLogData);
        for(int l = j; l < j+7*nVols; l += nVols){
            wX[l] = X[l]*W[j];
        }
    }
}

template<class T> inline void computeHatDiag(T* H, int nVols, const T* X, const T* wXinv, const T* W) {
	for(int j = 0; j < nVols; j++){
        H[j] = 0;
        int k = j*7;
        int l = j ; H[j] += X[l] * wXinv[k++] * W[j];
        l += nVols; H[j] += X[l] * wXinv[k++] * W[j];
        l += nVols; H[j] += X[l] * wXinv[k++] * W[j];
        l += nVols; H[j] += X[l] * wXinv[k++] * W[j];
        l += nVols; H[j] += X[l] * wXinv[k++] * W[j];
        l += nVols; H[j] += X[l] * wXinv[k++] * W[j];
        l += nVols; H[j] += X[l] * wXinv[k  ] * W[j];
    }
}

template<class T> inline void function(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    /* input */
    T* y = (T*) mxGetData(prhs[0]);
    T* X = (T*) mxGetData(prhs[1]);
    int nVols = mxGetM(prhs[0]);
    int nVox = mxGetN(prhs[0]);
    
    /* output */
    plhs[0] = mxCreateNumericMatrix(    7, nVox, mxGetClassID(prhs[0]), mxREAL);
    plhs[1] = mxCreateNumericMatrix(nVols, nVox, mxGetClassID(prhs[0]), mxREAL);
    plhs[2] = mxCreateNumericMatrix(nVols, nVox, mxGetClassID(prhs[0]), mxREAL);
    T* dt = (T*) mxGetData(plhs[0]);
    T* W = (T*) mxGetData(plhs[1]);
    T* H = (T*) mxGetData(plhs[2]);
    
    /* allocate workspace */
    T* wX = (T*) mxMalloc(nVols*7*sizeof(T));
    T* wXinv = (T*) mxMalloc(nVols*7*sizeof(T));
    T* dt_ = (T*) mxMalloc(7*sizeof(T));
    T* Xinv = (T*) mxMalloc(nVols*7*sizeof(T));
    T* Xcopy = (T*) mxMalloc(nVols*7*sizeof(T));
    for (int i = 0; i < nVols*7; i++) {
        Xcopy[i] = X[i];
    }
    
    /* actual work */
    pinv(Xcopy,Xinv,nVols,7);
    for (int i = 0; i < nVox; i++) {
        tensorfit<T>(dt_, nVols, NULL, &y[i*nVols], Xinv);
        computeWeightDiag<T>(&W[i*nVols], wX, nVols, X, dt_);
        pinv(wX,wXinv,nVols,7);
        computeHatDiag<T>(&H[i*nVols], nVols, X, wXinv, &W[i*nVols]);
        tensorfit<T>(&dt[i*7], nVols, &W[i*nVols], &y[i*nVols], wXinv);
    }
    
    /* deallocate workspace */
    mxFree(wX);
    mxFree(wXinv);
    mxFree(dt_);
    mxFree(Xinv);
    mxFree(Xcopy);
}

inline void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgTxt("need 2 input arguments!");
    }
    if (mxGetClassID(prhs[0]) != mxGetClassID(prhs[1])) {
        mexErrMsgTxt("types do not match!");
    }
    
    switch (mxGetClassID(prhs[0])) {
        case mxSINGLE_CLASS:
            function<float>(nlhs,plhs,nrhs,prhs);
            break;
        case mxDOUBLE_CLASS:
            function<double>(nlhs,plhs,nrhs,prhs);
            break;
        default:
            mexErrMsgTxt("type not supported yet!");
            break;
    }
}