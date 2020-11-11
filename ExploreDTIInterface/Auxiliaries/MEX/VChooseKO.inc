// VChooseKO.inc
// VChooseKO.INC - Helper function for VChooseKO.c
// Choose K elements without repetitions and with order.
// Instantiated macros:
//   FUNC_NAME: Name containg size of elements, e.g. Elem2_Byte2
//   DATA_TYPE: Type of the elements
//
// Author: Jan Simon, Heidelberg, (C) 2010 matlab.THISYEAR(a)nMINUSsimon.de
// License: BSD (use/copy/modify on own risk, but mention author)

/*
% $JRev: R0h V:008 Sum:iv64UU99b3nc Date:16-Jan-2010 02:01:51 $
% $File: Published\VChooseK_\VChooseKO.inc $
% History:
% 001: 14-Jan-2010 11:43, Core functions for array types as macro.
%      It is much safer to edit the different subfunction once only and not
%      separately for each input type.
%      The explicite loops for K==4 are 20% faster with LCC2.4 and 3.8, but
%      BCC5.5 and OpenWatcom 1.8 the general method has the same speed.
*/

// Prototypes:
void FUNC_NAME(2)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY);
void FUNC_NAME(3)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY);
void FUNC_NAME(4)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY);
void FUNC_NAME(K)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY, mwSize K);

// =============================================================================
void FUNC_NAME(2)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY)
{
  // Called for input data with 1/2/4/8 bytes per element:
  // (U)INT8,  LOGICAL, (U)INT16, CHAR, (U)INT32, SINGLE, (U)INT64, DOUBLE.
  // K == 2
  // INPUT:
  //   X: Pointer to data array of input.
  //   N: Number of elements of input vector.
  //   Y: Pointer to output array.
  
  DATA_TYPE *Y1, *Y2, S0, *Xp;
  mwSize i1, i2;
  
  Y1 = Y;
  Y2 = Y1 + nY;
  
  // Create the output:
  for (i1 = 0; i1 < N; i1++) {
     S0 = X[i1];
   
     // 2 loops: exclude i1==i2
     Xp = X;
     for (i2 = 0; i2 < i1; i2++) {
        *Y1++ = S0;
        *Y2++ = *Xp++;
     }
     
     Xp++;  // Skip X[i1]
     
     for (i2 = i1 + 1; i2 < N; i2++) {
        *Y1++ = S0;
        *Y2++ = *Xp++;
     }
  }
  
  return;
}

// =============================================================================
void FUNC_NAME(3)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY)
{
  // Called for input data with 1/2/4/8 bytes per element:
  // (U)INT8,  LOGICAL, (U)INT16, CHAR, (U)INT32, SINGLE, (U)INT64, DOUBLE.
  // K == 3
  // INPUT:
  //   X: Pointer to data array of input.
  //   N: Number of elements of input vector.
  //   Y: Pointer to output array.
  
  DATA_TYPE *Y1, *Y2, *Y3, S0, S1;
  UINT8_T *Q, *Qp, *Qf;
  mwSize i1, i2, i3;
  
  // Q: Vector of logical flags - TRUE: Element was chosen already:
  if ((Q = (UINT8_T *) mxMalloc(N * sizeof(UINT8_T))) == NULL) {
     mexWarnMsgIdAndTxt("JSimon:VChooseKO:NoMemory",
                        "Cannot allocate memory for [Q].");
  }
  
  // Initialize Q - all values are TRUE:
  Qf = Q + N;
  for (Qp = Q; Qp < Qf; *Qp++ = 1) ;  // empty loop
  
  Y1 = Y;
  Y2 = Y1 + nY;
  Y3 = Y2 + nY;
  
  // Create the output:
  for (i1 = 0; i1 < N; i1++) {
     Q[i1] = 0;
     S0    = X[i1];
     
     for (i2 = 0; i2 < N; i2++) {
        if (Q[i2]) {
           Q[i2] = 0;
           S1    = X[i2];
      
           Qp = Q;
           for (i3 = 0; i3 < N; i3++) {
              if (*Qp++) {
                 *Y1++ = S0;
                 *Y2++ = S1;
                 *Y3++ = X[i3];
              }
           }
           Q[i2] = 1;
        }
     }
     Q[i1] = 1;
  }
  
  mxFree(Q);
  
  return;
}

// =============================================================================
void FUNC_NAME(4)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY)
{
  // Called for input data with 1/2/4/8 bytes per element:
  // (U)INT8,  LOGICAL, (U)INT16, CHAR, (U)INT32, SINGLE, (U)INT64, DOUBLE.
  // K == 4
  // INPUT:
  //   X: Pointer to data array of input.
  //   N: Number of elements of input vector.
  //   Y: Pointer to output array.
  
  DATA_TYPE *Y1, *Y2, *Y3, *Y4, S0, S1, S2;
  UINT8_T *Q, *Qp, *Qf;
  mwSize i1, i2, i3, i4;
  
  // Q: Vector of logical flags - TRUE: Element was chosen already:
  if ((Q = (UINT8_T *) mxMalloc(N * sizeof(UINT8_T))) == NULL) {
     mexWarnMsgIdAndTxt("JSimon:VChooseKO:NoMemory",
                        "Cannot allocate memory for [Q].");
  }
  
  // Initialize Q - all values are TRUE:
  Qf = Q + N;
  for (Qp = Q; Qp < Qf; *Qp++ = 1) ;  // empty loop
  
  Y1 = Y;
  Y2 = Y1 + nY;
  Y3 = Y2 + nY;
  Y4 = Y3 + nY;
  
  // Create the output:
  for (i1 = 0; i1 < N; i1++) {
     Q[i1] = 0;
     S0    = X[i1];
     
     for (i2 = 0; i2 < N; i2++) {
        if (Q[i2]) {
           Q[i2] = 0;
           S1    = X[i2];
      
           for (i3 = 0; i3 < N; i3++) {
              if (Q[i3]) {
                 Q[i3] = 0;
                 S2    = X[i3];
                 
                 Qp = Q;
                 for (i4 = 0; i4 < N; i4++) {
                    if (*Qp++) {
                       *Y1++ = S0;
                       *Y2++ = S1;
                       *Y3++ = S2;
                       *Y4++ = X[i4];
                     }
                  }
                  Q[i3] = 1;
               }
           }
           Q[i2] = 1;
        }
     }
     Q[i1] = 1;
  }
  
  mxFree(Q);
  
  return;
}

// =============================================================================
void FUNC_NAME(K)(DATA_TYPE *X, mwSize N, DATA_TYPE *Y, mwSize nY, mwSize K)
{
  // Called for input data with 1/2/4/8 bytes per element:
  // (U)INT8,  LOGICAL, (U)INT16, CHAR, (U)INT32, SINGLE, (U)INT64, DOUBLE.
  // INPUT:
  //   X:  Pointer to data array of input.
  //   N:  Number of elements of input vector.
  //   Y: Pointer to output array.
  //   K:  Number of elements to choose, 2nd dimension of output.
  //       K > 1 for this lgorithm.
  
  DATA_TYPE *Yp, *Yarrayf, *Ycolumnf, *S, *Sp, *SKm1;
  mwSize *V, c, cK, jK, iN, Kminus1, Kminus2, Nminus1;
  UINT8_T *Q, *Qp, *Qf, found;
  
  // Some constants:
  Kminus1 = K - 1;
  Kminus2 = K - 2;
  Nminus1 = N - 1;
  
  // S: Vector of currently choosen values:
  if ((S = (DATA_TYPE *) mxMalloc(K * sizeof(DATA_TYPE))) == NULL) {
     mexWarnMsgIdAndTxt("JSimon:VChooseKO:NoMemory",
                        "Cannot allocate memory for [S].");
  }
  
  // Q: Vector of logical flags - TRUE: Element was chosen already:
  if ((Q = (UINT8_T *) mxMalloc(N * sizeof(UINT8_T))) == NULL) {
     mexWarnMsgIdAndTxt("JSimon:VChooseKO:NoMemory",
                        "Cannot allocate memory for [Q].");
  }
  
  // V: Vector of loop indices - instead of cascaded loops, the loop indices
  //    are managed in a vector for the general case.
  if ((V = (mwSize *) mxMalloc(K * sizeof(mwSize))) == NULL) {
     mexWarnMsgIdAndTxt("JSimon:VChooseKO:NoMemory",
                        "Cannot allocate memory for [V].");
  }
  
  // Initialize Q - the first K-1 values are TRUE:
  Qf = Q + Kminus1;
  for (Qp = Q; Qp < Qf; *Qp++ = 0) ;  // empty loop
  Qf = Q + N;
  for (      ; Qp < Qf; *Qp++ = 1) ;  // empty loop
  
  // Initialize V and S with 0:K-1 and X(0:K-1):
  for (c = 0; c < K; c++) {
    V[c] = c;                    // Loop indices
    S[c] = X[c];                 // Corresponding values of the input
  }

  // Create the output:
  SKm1     = S + Kminus1;
  Yarrayf  = Y + K * nY;         // End of output array
  Ycolumnf = Y + nY;             // End of first colum
  while (Y < Ycolumnf) {
     // Write S to row of output Y, loop over K.th index:
     Qp = Q;
     for (iN = 0; iN < N; iN++) {
        if (*Qp++) {             // Q[iN]
           *SKm1 = X[iN];        // Insert as last element of S
           Sp    = S;            // Copy S as row to output Y:
           for (Yp = Y++; Yp < Yarrayf; Yp += nY) {
             *Yp = *Sp++;
           }
           
           // NOTE: Remarkable slower here (but faster in VChooseK.c!):
           // A = vector of output column pointers;
           // for (i=0; i<K; i++) {**A++ = *Sp++; (*A)++;}
        }
     }
     
     // Advance indices as in a cascaded loop - rightmost index which is < N:
     c     = Kminus2;            // Zero based index to index K-1
     found = 0;
     while (c >= 0) {            // Stop if even 1st index is complete
        cK    = V[c];            // Value of c.th index
        Q[cK] = 1;               // Release this element in Q list
        if (cK < Nminus1) {      // Before zero based N.th element?
           for (jK = cK + 1; jK < N; jK++) {
              if (Q[jK]) {       // Search smallest unused index:
                 V[c]  = jK;
                 S[c]  = X[jK];
                 Q[jK] = 0;
                 found = 1;      // Stop while(C >= 0) loop
                 break;          // Stop for(jK) loop
              }
           }
           if (found) {
              break;
           }
        }
        c--;                     // Check next index to the left
     }
     
     // Restart higher loops with lowest available index:
     if (c >= 0) {               // Not if calculations are complete
        for (c = c + 1; c < Kminus1; c++) {
           Qp = Q;
           for (jK = 0; jK < N; jK++) {
              if (*Qp++) {       // Q[jK] - This input is not used already
                 V[c]  = jK;     // Remember the index
                 S[c]  = X[jK];  // Insert value in S
                 Q[jK] = 0;      // Set flag: now this input is used
                 break;          // Stop for(jK)
              }
           }
        }
     }
  }  // while Y < Ycolumnf
  
  mxFree(S);
  mxFree(Q);
  mxFree(V);
  
  return;
}
