/*
 Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
 None of this code can be copied, reused or distributed without explicit permission from the author.
*/
#include <mex.h>
#include <cmath>
template<class T> inline void function(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    const T PI = 2*acos(0.0);
    const mwSize dt_ndim = mxGetNumberOfDimensions(prhs[0]);
	const mwSize* dt_dims = mxGetDimensions(prhs[0]);
 	mxClassID dt_classid = mxGetClassID(prhs[0]);
    T* dt = (T*) mxGetData(prhs[0]);
    
    mwSize eigval_dims[2] = {3,dt_dims[1]};
    mwSize eigvec_dims[2] = {9,dt_dims[1]};
    
    plhs[0] = mxCreateNumericMatrix(0,0,dt_classid,mxREAL);
    mxSetDimensions(plhs[0], eigval_dims, 2);
    mxSetData(plhs[0], mxMalloc(sizeof(T)*3*eigval_dims[1]));
    T* eigval = (T*) mxGetData(plhs[0]);

    T* eigvec;
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericMatrix(0,0,dt_classid,mxREAL);
        mxSetDimensions(plhs[1], eigvec_dims, 2);
        mxSetData(plhs[1], mxMalloc(sizeof(T)*9*eigvec_dims[1]));
        eigvec = (T*) mxGetData(plhs[1]);
    }

    T* dt_i;
    T* eigval_i;
    T* eigvec_i;
    T I1, I2, I3, v, s, phi, A, B, C, n;
            
    for (mwSize i = 0; i < dt_dims[1]; i++) {
        dt_i = &dt[i*dt_dims[0]];
        eigval_i = &eigval[i*eigval_dims[0]];

        I1 = dt_i[0] + dt_i[3] + dt_i[5];
        I2 = dt_i[0]*dt_i[3] + dt_i[0]*dt_i[5] + dt_i[3]*dt_i[5] - (dt_i[1]*dt_i[1] + dt_i[2]*dt_i[2] + dt_i[4]*dt_i[4]);
        I3 = dt_i[0]*dt_i[3]*dt_i[5] + 2*dt_i[1]*dt_i[2]*dt_i[4] - (dt_i[1]*dt_i[1]*dt_i[5] + dt_i[2]*dt_i[2]*dt_i[3] + dt_i[0]*dt_i[4]*dt_i[4]);

        v = I1*I1/9 - I2/3;
        s = I1*I1*I1/27 - I1*I2/6 + I3/2;
        phi = acos((s/v)*(1/sqrt(v)))/3;

        eigval_i[0] = I1/3 + 2*sqrt(v)*cos(phi);
        eigval_i[1] = I1/3 - 2*sqrt(v)*cos(PI/3+phi);
        eigval_i[2] = I1/3 - 2*sqrt(v)*cos(PI/3-phi);

        if (nlhs > 1) {
            eigvec_i = &eigvec[i*eigvec_dims[0]];
            
            A = dt_i[0] - eigval_i[0]; B = dt_i[3] - eigval_i[0]; C = dt_i[5] - eigval_i[0];
            eigvec_i[0] = (dt_i[1]*dt_i[4] - B*dt_i[2])*(dt_i[2]*dt_i[4] - C*dt_i[1]);
            eigvec_i[1] = (dt_i[2]*dt_i[4] - C*dt_i[1])*(dt_i[2]*dt_i[1] - A*dt_i[4]);
            eigvec_i[2] = (dt_i[1]*dt_i[4] - B*dt_i[2])*(dt_i[2]*dt_i[1] - A*dt_i[4]);
            n = sqrt(eigvec_i[0]*eigvec_i[0]+eigvec_i[1]*eigvec_i[1]+eigvec_i[2]*eigvec_i[2]);
            eigvec_i[0] /=n; eigvec_i[1] /=n; eigvec_i[2] /=n;

            A = dt_i[0] - eigval_i[1]; B = dt_i[3] - eigval_i[1]; C = dt_i[5] - eigval_i[1];
            eigvec_i[3] = (dt_i[1]*dt_i[4] - B*dt_i[2])*(dt_i[2]*dt_i[4] - C*dt_i[1]);
            eigvec_i[4] = (dt_i[2]*dt_i[4] - C*dt_i[1])*(dt_i[2]*dt_i[1] - A*dt_i[4]);
            eigvec_i[5] = (dt_i[1]*dt_i[4] - B*dt_i[2])*(dt_i[2]*dt_i[1] - A*dt_i[4]);
            n = sqrt(eigvec_i[3]*eigvec_i[3]+eigvec_i[4]*eigvec_i[4]+eigvec_i[5]*eigvec_i[5]);
            eigvec_i[3] /=n; eigvec_i[4] /=n; eigvec_i[5] /=n;


            eigvec_i[6] = eigvec_i[1]*eigvec_i[5]-eigvec_i[4]*eigvec_i[2];
            eigvec_i[7] = eigvec_i[2]*eigvec_i[3]-eigvec_i[0]*eigvec_i[5];
            eigvec_i[8] = eigvec_i[0]*eigvec_i[4]-eigvec_i[1]*eigvec_i[3];
    //         A = dt_i[0] - eigval_i[2]; B = dt_i[3] - eigval_i[2]; C = dt_i[5] - eigval_i[2];
    //         eigvec_i[6] = (dt_i[1]*dt_i[4] - B*dt_i[2])*(dt_i[2]*dt_i[4] - C*dt_i[1]);
    //         eigvec_i[7] = (dt_i[2]*dt_i[4] - C*dt_i[1])*(dt_i[2]*dt_i[1] - A*dt_i[4]);
    //         eigvec_i[8] = (dt_i[1]*dt_i[4] - B*dt_i[2])*(dt_i[2]*dt_i[1] - A*dt_i[4]);
    //         n = sqrt(eigvec_i[6]*eigvec_i[6]+eigvec_i[7]*eigvec_i[7]+eigvec_i[8]*eigvec_i[8]);
    //         eigvec_i[6] /=n; eigvec_i[7] /=n; eigvec_i[8] /=n;
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs != 1) {
        mexErrMsgTxt("requires one argument!");
    }
    const mwSize dt_ndim = mxGetNumberOfDimensions(prhs[0]);
	const mwSize* dt_dims = mxGetDimensions(prhs[0]);
 	mxClassID dt_classid = mxGetClassID(prhs[0]);
	if (dt_ndim != 2) {
        mexErrMsgTxt("input must be a 6 x n matrix");
    }
    if (dt_dims[0] != 6) {
        mexErrMsgTxt("input must be a 6 x n matrix");
    }
    switch (dt_classid) {
        case mxDOUBLE_CLASS:
            function<double>(nlhs,plhs,nrhs,prhs);
            break;
        case mxSINGLE_CLASS:
            function<float>(nlhs,plhs,nrhs,prhs);
            break;
        default:
            mexErrMsgTxt("type not supported yet!");
            break;
    }
}