% VChooseKO - Permutations of K elements, (order, no repetitions) [MEX]
% VChooseKO(V, K) creates a matrix, which rows are all permutations of
% choosing K elements of the vector V with order and without repetitions.
%
% INPUT:
%   V: Array of class DOUBLE, SINGLE, (U)INT8/16/32/64, LOGICAL, CHAR.
%      V can have any shape.
%   K: Number of elements to choose. Scalar DOUBLE with integer value.
%      0 <= K <= NUMEL(V).
%
% OUTPUT:
%   Y: [N!/(N-K)!, K] matrix with N is the number of elements of V.
%      Y has the same class as the input V.
%      Note: N!/(N-K)! = N*(N-1)*..*(N-K+1)
%      The rows are sorted in lexicographical order: smaller indices at first.
%
% NOTE: A warning appears, if the output exceeds 500MB, and an error for 1GB.
%   Both limits can be adjusted according to the available RAM in the C-Mex
%   source.
%
% EXAMPLES:
%   Choose 2 elements from [1, 2, 3]:
%     VChooseKO(1:3, 2)  % ==> [1,2; 1,3; 2,1; 2,3; 3,1; 3,2]
%   For speed cast the input to integer types or SINGLE whenever possible:
%     Y = VChooseKO(uint8(1:100), 3);  % 5 times faster than:
%     Y = VChooseKO(1:100, 3);
%   To get the permutations of cell arrays, permute the index:
%     C  = {'a', 'b', 'c', 'd'};
%     C2 = C(VChooseKO(uint8(1:4), 2))
%     ==> C2 = {'a','b'; 'a','c'; 'a','d'; 'b','a'; 'b','c'; 'b','d'; ...
%               'c','a'; 'c','b'; 'c','d'; 'd','a'; 'd','b'; 'd','c'}
%   Equivalent to PERMS:
%     isequal(sortrows(perms(1:5)), VChooseKO(1:5, 5))   % TRUE
%   For an output with sorted values (not indices!), sort the input:
%     X = [3, 1, 2];  VChooseKO(sort(X), 3)
%     ==> [1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1]
%
% COMPILE:
%   Windows: mex -O VChooseKO.c
%   Linux:   mex CFLAGS="\$CFLAGS -std=C99" -O VChooseKO.c
%   Please run the unit-test TestVChooseKO after compiling!
%   More help in: VChooseKO.c
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
%         Compilers: BCC5.5, LCC2.4/3.8, Open Watcom 1.8
% Author: Jan Simon, Heidelberg, (C) 2010 matlab.THIS_YEAR(a)nMINUSsimonDOTde
% License: BSD (use/copy/modify on own risk, but mention author)
%
% See also: NCHOOSEK, PERMS, PERMUTE.
% FEX: COMBINATOR, NPERMUTEK, COMBINATIONS, COMBN,
%      VCHOOSEK, VCHOOSEKR, VCHOOSEKRO.

% $JRev: R0m V:012 Sum:Hl7laSwVXeHN Date:16-Jan-2010 02:02:10 $
% $File: Published\VChooseK_\VChooseKO.m $
% History:
% 001: 14-Jan-2010 15:39, Finally the last of the 4 choosing programs.

% This is a dummy file to support the help command.
% The unit-test function TestVChooseKO contains a Matlab implementation of the
% used algorithm.
% More help and notes: VChooseKO.c
