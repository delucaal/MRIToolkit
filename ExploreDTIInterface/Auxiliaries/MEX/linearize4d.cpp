#include <mex.h>
#include <limits>
// Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
// under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

inline int offset2(int n1, int n2, int N1) {
    return n1+N1*n2;
}

inline int offset3(int n1, int n2, int n3, int N1, int N2) {
    return n1+N1*(n2+N2*n3);
}

inline int offset4(int n1, int n2, int n3, int n4, int N1, int N2, int N3) {
    return n1+N1*(n2+N2*(n3+N3*n4));
}

template<class T> inline void function(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    const mwSize* unlin_size = mxGetDimensions(prhs[0]);
    const mwSize* mask_size = mxGetDimensions(prhs[1]);
    T* unlin = (T*) mxGetData(prhs[0]);
    bool* mask = (bool*) mxGetData(prhs[1]);
    
    int count = 0;
    for (int i = 0; i < mask_size[0]*mask_size[1]*mask_size[2]; i++) {
        if (mask[i]) {
            count++;
        }
    }
    const mwSize lin_size[] = {unlin_size[3],count};
    
    plhs[0] = (mxArray*) mxCreateNumericArray(2,lin_size,mxGetClassID(prhs[0]),mxREAL);
    T* lin = (T*) mxGetData(plhs[0]);
    
    for (int w = 0; w < unlin_size[3]; w++) {
        int count = 0;
        for (int z = 0; z < unlin_size[2]; z++) {
            for (int y = 0; y < unlin_size[1]; y++) {
                for (int x = 0; x < unlin_size[0]; x++) {
                    if (mask[offset3(x,y,z,mask_size[0],mask_size[1])]) {
                        lin[offset2(w,count,lin_size[0])] = unlin[offset4(x,y,z,w,unlin_size[0],unlin_size[1],unlin_size[2])];
                        count++;
                    }
                    
                    
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    switch (mxGetClassID(prhs[0])) {
        case mxINT8_CLASS:
            function<char>(nlhs,plhs,nrhs,prhs);
            break;
        case mxUINT8_CLASS:
            function<unsigned char>(nlhs,plhs,nrhs,prhs);
            break;
        case mxINT16_CLASS:
            function<short>(nlhs,plhs,nrhs,prhs);
            break;
        case mxUINT16_CLASS:
            function<unsigned short>(nlhs,plhs,nrhs,prhs);
            break;
        case mxINT32_CLASS:
            function<int>(nlhs,plhs,nrhs,prhs);
            break;
        case mxUINT32_CLASS:
            function<unsigned int>(nlhs,plhs,nrhs,prhs);
            break;
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