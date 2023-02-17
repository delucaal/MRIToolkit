/*
 Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
 None of this code can be copied, reused or distributed without explicit permission from the author.
*/

#include <mex.h>
// #include <limits>

template<class T> inline void function(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray* p_data = prhs[0];
    const mwSize data_ndim = mxGetNumberOfDimensions(p_data);
	const mwSize* data_dims = mxGetDimensions(p_data);
    const mwSize data_nels = mxGetNumberOfElements(p_data);
    T* data = (T*) mxGetData(p_data);
    
//     mexPrintf("%i\n",data_ndim);
//     for (mwSize i = 0; i < data_ndim; i++) {
//         mexPrintf("%i ",data_dims[i]);
//     }
//     mexPrintf("\n");
//     mexPrintf("%i\n",data_nels);
    
    const mxArray* p_mask = prhs[1];
    const mwSize mask_ndim = mxGetNumberOfDimensions(p_mask);
	const mwSize* mask_dims = mxGetDimensions(p_mask);
    const mwSize mask_nels = mxGetNumberOfElements(p_mask);
    bool* mask = (bool*) mxGetData(p_mask);
    
//     mexPrintf("%i\n",mask_ndim);
//     for (mwSize i = 0; i < mask_ndim; i++) {
//         mexPrintf("%i ",mask_dims[i]);
//     }
//     mexPrintf("\n");
//     mexPrintf("%i\n",mask_nels);
    
    const mxArray* p_output = plhs[0];
    const mwSize output_ndim = mxGetNumberOfDimensions(p_output);
	const mwSize* output_dims = mxGetDimensions(p_output);
    const mwSize output_nels = mxGetNumberOfElements(p_output);
    T* output = (T*) mxGetData(p_output);        
 
//     mexPrintf("%i\n",output_ndim);
//     for (mwSize i = 0; i < output_ndim; i++) {
//         mexPrintf("%i ",output_dims[i]);
//     }
//     mexPrintf("\n");
//     mexPrintf("%i\n",output_nels);

    mwSize output_dims3 = output_dims[3];;
    if (output_ndim == 3) {
        output_dims3 = 1;
    }
    
    mwSize ii = 0;
    for (mwSize k = 0; k < output_dims[2]; k++) {
        for (mwSize j = 0; j < output_dims[1]; j++) {
            for (mwSize i = 0; i < output_dims[0]; i++) {
                if (mask[i+mask_dims[0]*(j+mask_dims[1]*k)]) {
                    for (mwSize l = 0; l < output_dims3; l++) {
                        output[i+output_dims[0]*(j+output_dims[1]*(k+output_dims[2]*l))] = data[ii];
                        ii++;
                        if (ii == data_nels) {
                            return;
                        }
                    }
                }
            }
        }
    }
}




void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgTxt("requires two arguments!");
    }
    const mxArray* p_data = prhs[0];
    const mxArray* p_mask = prhs[1];
    
    const mwSize mask_ndim = mxGetNumberOfDimensions(p_mask);
	const mwSize* mask_dims = mxGetDimensions(p_mask);
    
    const mwSize data_ndim = mxGetNumberOfDimensions(p_data);
	const mwSize* data_dims = mxGetDimensions(p_data);
    
	mxClassID mask_classid = mxGetClassID(p_mask);
    mxClassID data_classid = mxGetClassID(p_data);
    
    mwSize new_dims[4];
    new_dims[0] = mask_dims[0];
    new_dims[1] = mask_dims[1];
    new_dims[2] = mask_dims[2];
    new_dims[3] = data_dims[0];
    
   
    if (mask_classid != mxLOGICAL_CLASS) {
        mexErrMsgTxt("mask must be a 3D logical matrix");
    }
    
    if (mask_ndim != 3) {
        mexErrMsgTxt("mask must be a 3D logical matrix");
    }
    
    if (data_ndim != 2) {
        mexErrMsgTxt("input must be a 2D matrix");
    }
    
    plhs[0] = mxCreateNumericArray(mask_ndim+1,new_dims,data_classid,mxREAL);
    switch (data_classid) {
        case mxDOUBLE_CLASS:
            function<double>(nlhs,plhs,nrhs,prhs);
            break;
        case mxSINGLE_CLASS:
            function<float>(nlhs,plhs,nrhs,prhs);
            break;
        case mxINT64_CLASS:
            function<int64_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxUINT64_CLASS:
            function<uint64_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxINT32_CLASS:
            function<int32_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxUINT32_CLASS:
            function<uint32_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxINT16_CLASS:
            function<int16_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxUINT16_CLASS:
            function<uint16_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxINT8_CLASS:
            function<int8_T>(nlhs,plhs,nrhs,prhs);
            break;            
        case mxUINT8_CLASS:
            function<uint8_T>(nlhs,plhs,nrhs,prhs);
            break;
        case mxLOGICAL_CLASS:
            function<bool>(nlhs,plhs,nrhs,prhs);
            break;  
        default:
            mexErrMsgTxt("type not supported yet!");
            break;
    }
}
