// VChooseKO.c
// VChooseKO - Permutations of K elements, (order, no repetitions) [MEX]
// VChooseKO(V, K) creates a matrix, which rows are all permutations of
// choosing K elements of the vector V with order and without repetitions.
//
// INPUT:
//   V: Array of class DOUBLE, SINGLE, (U)INT8/16/32/64, LOGICAL, CHAR.
//      V can have any shape.
//   K: Number of elements to choose. Scalar DOUBLE with integer value.
//      0 <= K <= NUMEL(V).
//
// OUTPUT:
//   Y: [N!/(N-K)!, K] matrix with N is the number of elements of V.
//      Y has the same class as the input V.
//      The rows are sorted in lexicographical order: smaller indices at first.
//      Note: N!/(N-K)! = N*(N-1)*..*(N-K+1)
//
// NOTE: A warning appears, if the output exceeds 500MB, and an error for 1GB.
//   Both limits can be adjusted according to the available RAM in the C-Mex
//   source.
//
// EXAMPLES:
//   Choose 2 elements from [1, 2, 3]:
//     VChooseKO(1:3, 2)  % ==> [1,2; 1,3; 2,1; 2,3; 3,1; 3,2]
//   For speed cast the input to integer types or SINGLE whenever possible:
//     Y = VChooseKO(uint8(0:255), 3);  % 5 times faster than:
//     Y = VChooseKO(0:255, 3);
//   To get the permutations of cell arrays, permute the index:
//     C  = {'a', 'b', 'c', 'd'};
//     C2 = C(VChooseKO(uint8(1:4), 2))
//     ==> C2 = {'a','b'; 'a','c'; 'a','d'; 'b','a'; 'b','c'; 'b','d'; ...
//               'c','a'; 'c','b'; 'c','d'; 'd','a'; 'd','b'; 'd','c'}
//   Equivalent to PERMS:
//     isequal(sortrows(perms(1:5)), VChooseKO(1:5, 5))   % TRUE
//   For an output with sorted values (not indices!), sort the input:
//     X = [3, 1, 2];  VChooseKO(sort(X), 3)
//     ==> [1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1]
//
// COMPILE:
//   Windows: mex -O VChooseKO.c
//   Linux:   mex CFLAGS="\$CFLAGS -std=C99" -O VChooseKO.c
//            (for C99 comments - thanks to Sebastiaan Breedveld)
//   Precompiled MEX: http://www.n-simon.de/mex
//   Compatibility to 64-bit machines is assumed, but not tested.
//   Please run the unit-test TestVChooseKO after compiling!
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
//         Compilers: BCC5.5, LCC2.4/3.8, Open Watcom 1.8
//         Watcom lib is 20% faster than BCC and LCC libs for INT8!
// Author: Jan Simon, Heidelberg, (C) 2010 matlab.THIS_YEAR(a)nMINUSsimonDOTde
// License: BSD (use/copy/modify on own risk, but mention author)
//
// See also: NCHOOSEK, PERMS.
// FEX: COMBINATOR, NPERMUTEK, COMBINATIONS, COMBN,
//      VCHOOSEK, VCHOOSEKR, VCHOOSEKRO.

// INSPIRATION:
// After VChooseK, VChooseKR, VChooseKRO, it was impossible not to publish
// VChooseKO - although I never needed it personally.
// At first I followed Knuth (The Art of Computer Programming):
//   function Y = perm3(V, k)
//   n  = length(V);
//   nY = prod((n - k + 1):n);
//   Y  = zeros(nY, k);
//   F  = fliplr(cumprod(1:n-1));
//   iY = 0;
//   while iY < nY
//      S = V;
//      for j = 1:n - 1
//         tempj = mod(floor(iY / F(j)), n + 1 - j);
//         if tempj
//            jt        = j + tempj;
//            temps     = S(jt);
//            S(j+1:jt) = S(j:jt-1);
//            S(j)      = temps;
//         end
//      end
//      iY       = iY + 1;
//      Y(iY, :) = S;
//   end
// (See also: FEX: FACTORADIC, 2009, Darren Rowland)
// But even with tricks (reusing already shifted chains, loop unrolling) this
// is 10 times slower in C than the trivial nested loop approach:
//   for i1 = 1:n
//      S(1) = V(i1)
//      for i2 = 1:n
//         if i1 ~= i2
//            S(2) = V(i2);
//            for i3 = 1:n
//               if i3 ~= i2 && i3 ~= i1
//                  etc;
//                  Y(iY, :) = S;
//   end; end; end; end; ... end
// Even the method of simulating the nested loops by a vector of loop indices
// for free K, this is much faster than the shifting. See VChooseKO.inc.
// And in addition the nested loop method yields to a lexicographical order!
//
// Other related publications in the FEX:
// (http://www.mathworks.com/matlabcentral/fileexchange/<number>)
//   COMBINATOR (Matt Fig) [Fast and general Matlab implementation!]: 24325
//   NPERMUTEK (Matt Fig): 11462
//   COMBINATIONS (Gautam Vallabha): 23080
//   COMBN (Jos van der Geest): 7147
//   VCHOOSEK (Jan Simon) no order, no repetitions: 26190
//   VCHOOSEKR (Jan Simon) no order, repetitions: 26277
//   VChooseKRO (Jan Simon) order, repetitions: 26242
//   VCHOOSEKO (Jan Simon) order, no repetitions: ?

/*
% $JRev: R0i V:015 Sum:2fb5EEqeopST Date:17-Jan-2010 17:22:05 $
% $File: Published\VChooseK_\VChooseKO.c $
% History:
% 001: 14-Jan-2010 13:25, First version. MEX interface taken from VChooseK.c.
*/

#include "mex.h"
#include <math.h>
#include <string.h>

// Assume 32 bit array dimensions for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Limits for warning and error:
#define WARN_LIMIT   500000000L   // 500 MB, *must* be smaller than ERROR_LIMIT
#define ERROR_LIMIT 1000000000L   // 1 GB (fair for 32&64 bit systems)

// Error messages do not contain the function name in Matlab 6.5!
#define ERR_FUNC "*** VChooseKO[MEX]: "

// Prototypes:
void BadInputTypeError(void);
mwSize GetOutputSize(mwSize n, mwSize k, double *C_double);

// Core functions: =============================================================
// Include the core function for different input types:

#define DATA_TYPE int8_T
#define FUNC_NAME(X) Elem ## X ## _Byte1
#include "VChooseKO.inc"
#undef DATA_TYPE
#undef FUNC_NAME

#define DATA_TYPE int16_T
#define FUNC_NAME(X) Elem ## X ## _Byte2
#include "VChooseKO.inc"
#undef DATA_TYPE
#undef FUNC_NAME

#define DATA_TYPE int32_T
#define FUNC_NAME(X) Elem ## X ## _Byte4
#include "VChooseKO.inc"
#undef DATA_TYPE
#undef FUNC_NAME

#define DATA_TYPE double         // or int64_T
#define FUNC_NAME(X) Elem ## X ## _Byte8
#include "VChooseKO.inc"
#undef DATA_TYPE
#undef FUNC_NAME

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize N, nY, K, ElementSize;
  const mxArray *X;
  double nY_double;
  
  // Check number of inputs and outputs:
  if (nrhs != 2) {
     mexErrMsgIdAndTxt("JSimon:VChooseKO:BadNInput",
                       ERR_FUNC "2 inputs required.");
  }
  if (nlhs > 1) {
     mexErrMsgIdAndTxt("JSimon:VChooseKO:BadNInput",
                       ERR_FUNC "1 output allowed.");
  }
  
  // Check type of input:
  X = prhs[0];
  if (!mxIsNumeric(X) && !mxIsChar(X) & !mxIsLogical(X)) {
     mexErrMsgIdAndTxt("JSimon:VChooseKO:BadNInput",
                       ERR_FUNC
                       "Input array must be numerical, CHAR or LOGICAL.");
  }
  if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) {
     mexErrMsgIdAndTxt("JSimon:VChooseKO:BadNInput",
                       ERR_FUNC "Input K must be a numerical scalar.");
  }
  
  // Get number of elements to choose:
  N = mxGetNumberOfElements(X);
  K = (mwSize) floor(mxGetScalar(prhs[1]) + 0.5);  // Explicite rounding
  
  // Empty input => empty output:
  if (N == 0) {
     plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     return;
  }
  
  // Degenerated cases:
  //   K == 1: reply input as column vector
  //   K == 0: empty output
  //   K < 0:  error
  if (K <= 1 || K > N) {
     if (K == 1) {           // Choose 1 element => output equals input:
        plhs[0] = mxDuplicateArray(X);
        mxSetM(plhs[0], N);  // Column vector
        mxSetN(plhs[0], 1);
     } else if (K == 0) {    // Empty matrix:
        plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(X), mxREAL);
     } else {
        mexErrMsgIdAndTxt("JSimon:VChooseK:BadK",
                          ERR_FUNC "Input K: 0 <= K <= N.");
     }
     return;
  }
  
  // Calculate size of output N!/(N-K)!
  nY = GetOutputSize(N, K, &nY_double);
  
  // Check limits - use nY as DOUBLE to control overflows:
  ElementSize  = mxGetElementSize(X);
  nY_double   *= ElementSize;
  if (nY_double > WARN_LIMIT) {
     if (nY_double < ERROR_LIMIT && nY_double < MWSIZE_MAX) {
        mexWarnMsgIdAndTxt("JSimon:VChooseKO:LargeOutput",
                           ERR_FUNC "Output will be large and slow.");
     } else {
        mexErrMsgIdAndTxt("JSimon:VChooseKO:TooLargeOutput",
                          ERR_FUNC "Output would be too large: "
                          "(%.8g x %d) * %d Byte", nY_double, K, ElementSize);
     }
  }
  
  // Create the output:
  // In a MEX file this stops with an error message automatically on failure.
  // For usage in an engine file, check result for NULL!
  plhs[0] = mxCreateNumericMatrix(nY, K, mxGetClassID(X), mxREAL);
  
  // Call different core functions according to the element size and K:
  switch (K) {
     case 2:
        switch (ElementSize) {
           case 8:  Elem2_Byte8(mxGetPr(X), N, mxGetPr(plhs[0]), nY);
                    break;
           case 4:  Elem2_Byte4((int32_T *) mxGetData(X), N,
                                (int32_T *) mxGetData(plhs[0]), nY);
                    break;
           case 2:  Elem2_Byte2((int16_T *) mxGetData(X), N,
                                (int16_T *) mxGetData(plhs[0]), nY);
                    break;
           case 1:  Elem2_Byte1((int8_T *) mxGetData(X), N,
                                (int8_T *) mxGetData(plhs[0]), nY);
                    break;
           default: BadInputTypeError();
        }
        break;
        
     case 3:
        switch (ElementSize) {
           case 8:  Elem3_Byte8(mxGetPr(X), N, mxGetPr(plhs[0]), nY);
                    break;
           case 4:  Elem3_Byte4((int32_T *) mxGetData(X), N,
                                (int32_T *) mxGetData(plhs[0]), nY);
                    break;
           case 2:  Elem3_Byte2((int16_T *) mxGetData(X), N,
                                (int16_T *) mxGetData(plhs[0]), nY);
                    break;
           case 1:  Elem3_Byte1((int8_T *) mxGetData(X), N,
                                (int8_T *) mxGetData(plhs[0]), nY);
                    break;
           default: BadInputTypeError();
        }
        break;
        
     case 4:
        // With LCC3, this is 20% faster than general method with free K.
        // But with the Open Watcom 1.8 both methods have equivalent speed.
        switch (ElementSize) {
           case 8:  Elem4_Byte8(mxGetPr(X), N, mxGetPr(plhs[0]), nY);
                    break;
           case 4:  Elem4_Byte4((int32_T *) mxGetData(X), N,
                                (int32_T *) mxGetData(plhs[0]), nY);
                    break;
           case 2:  Elem4_Byte2((int16_T *) mxGetData(X), N,
                                (int16_T *) mxGetData(plhs[0]), nY);
                    break;
           case 1:  Elem4_Byte1((int8_T *) mxGetData(X), N,
                                (int8_T *) mxGetData(plhs[0]), nY);
                    break;
           default: BadInputTypeError();
        }
        break;
             
     default:  // General case for free K:
        switch (ElementSize) {
           case 8:  ElemK_Byte8(mxGetPr(X), N, mxGetPr(plhs[0]), nY, K);
                    break;
           case 4:  ElemK_Byte4((int32_T *) mxGetData(X), N,
                                (int32_T *) mxGetData(plhs[0]), nY, K);
                    break;
           case 2:  ElemK_Byte2((int16_T *) mxGetData(X), N,
                                (int16_T *) mxGetData(plhs[0]), nY, K);
                    break;
           case 1:  ElemK_Byte1((int8_T *) mxGetData(X), N,
                                (int8_T *) mxGetData(plhs[0]), nY, K);
                    break;
           default: BadInputTypeError();
        }
  }
  
  return;
}

// Error about bad input type: =================================================
void BadInputTypeError(void)
{
  mexErrMsgIdAndTxt("JSimon:VChooseKO:BadInputType",
                    ERR_FUNC "Input must have 1, 2, 4 or 8 bytes per element.");
}

// Get output size: ============================================================
mwSize GetOutputSize(mwSize n, mwSize k, double *nY_double)
{
  // Number of samples for taking K elements from a set of N without repetitions
  // and with order: N!/(N-K)! = N*(N-1)*..*(N-K+1)
  // Calculate the exact value as int32 (or int64 on 64 bit machine) and the
  // eventually less precise DOUBLE value, which is less susceptible for an
  // overflow.
  
  mwSize nY = 1, i;
  double dY = 1.0;
  
  for (i = n - k + 1; i <= n; i++) {
     nY *= i;
     dY *= (double) i;
  }
  
  *nY_double = dY;
  return (nY);
}
