#include <mex.h>
#include <limits>
// Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
// under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
template<class T> inline void function(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    const mwSize* func_size = mxGetDimensions(prhs[0]);
    const mwSize* points_size = mxGetDimensions(prhs[1]);
    T* func = (T*) mxGetData(prhs[0]);
    T* points = (T*) mxGetData(prhs[1]);
    
    plhs[0] = mxCreateNumericMatrix(func_size[3],points_size[1],mxGetClassID(prhs[0]),mxREAL);
    const mwSize* values_size = mxGetDimensions(plhs[0]);
    T* values = (T*) mxGetData(plhs[0]);
    
    for (int pos = 0; pos < points_size[1]; pos++) {
        T* point = &points[pos*points_size[0]];
        T* value = &values[pos*values_size[0]];
      
        bool oob = false;
        
        int min0 = int(point[0]); if (min0 < 1) oob = true;
        int min1 = int(point[1]); if (min1 < 1) oob = true;
        int min2 = int(point[2]); if (min2 < 1) oob = true;
        
        int max0 = min0 + 1; if (max0 > func_size[0]) oob = true;
        int max1 = min1 + 1; if (max1 > func_size[1]) oob = true;
        int max2 = min2 + 1; if (max2 > func_size[2]) oob = true;
        
        T d0 = point[0]-min0;
        T d1 = point[1]-min1;
        T d2 = point[2]-min2;
       
        T w0 = (1-d0)*(1-d1)*(1-d2); if (w0 < 1e-6) w0 = 0.0;
        T w1 = (1-d0)*(1-d1)*   d2 ; if (w1 < 1e-6) w1 = 0.0;
        T w2 = (1-d0)*   d1 *(1-d2); if (w2 < 1e-6) w2 = 0.0;
        T w3 =    d0 *(1-d1)*(1-d2); if (w3 < 1e-6) w3 = 0.0;
        T w4 = (1-d0)*   d1 *   d2 ; if (w4 < 1e-6) w4 = 0.0;
        T w5 =    d0 *(1-d1)*   d2 ; if (w5 < 1e-6) w5 = 0.0;
        T w6 =    d0 *   d1 *(1-d2); if (w6 < 1e-6) w6 = 0.0;
        T w7 =    d0 *   d1 *   d2 ; if (w7 < 1e-6) w7 = 0.0;
      
        min0--; min1--; min2--;
        max0--; max1--; max2--;
       
        if (!oob) {
            for (int i = 0; i < values_size[0]; i++) {
                value[i] = 
                  w0*func[min0+func_size[0]*(min1+func_size[1]*(min2+func_size[2]*i))]
                + w1*func[min0+func_size[0]*(min1+func_size[1]*(max2+func_size[2]*i))]
                + w2*func[min0+func_size[0]*(max1+func_size[1]*(min2+func_size[2]*i))]
                + w3*func[max0+func_size[0]*(min1+func_size[1]*(min2+func_size[2]*i))]
                + w4*func[min0+func_size[0]*(max1+func_size[1]*(max2+func_size[2]*i))]
                + w5*func[max0+func_size[0]*(min1+func_size[1]*(max2+func_size[2]*i))]
                + w6*func[max0+func_size[0]*(max1+func_size[1]*(min2+func_size[2]*i))]
                + w7*func[max0+func_size[0]*(max1+func_size[1]*(max2+func_size[2]*i))];
            }
        } else {
            for (int i = 0; i < values_size[0]; i++) {
                value[i] = std::numeric_limits<T>::quiet_NaN();
            }
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgTxt("types do not match!");
    }
    if ((mxGetClassID(prhs[0]) != mxGetClassID(prhs[1]))) {
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
