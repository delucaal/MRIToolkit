%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



%-fanDTasia ToolBox------------------------------------------------------------------
% This Matlab script is part of the fanDTasia ToolBox: a Matlab library for Diffusion 
% Weighted MRI (DW-MRI) Processing, Diffusion Tensor (DTI) Estimation, High-order 
% Diffusion Tensor Analysis, Tensor ODF estimation, Visualization and more.
%
% A Matlab Tutorial on DW-MRI can be found in:
% http://www.cise.ufl.edu/~abarmpou/lab/fanDTasia/tutorial.php
%
%----------------------------------------------
%writeFDT
%This function writes an image volume in FDT file format.
%example1: writeFDT('DTI.fdt', A)
%A can be 1D matrix for 1D data, 2D for a 2D image, 3D for a volume, 4D for a set of volume images.
%
%example2: writeFDT('DTI.fdt', A, B)
%A can be 1D matrix for 1D data, 2D for a 2D image, 3D for a volume, 4D for a set of volume images.
%B is a Nx4 matrix that contains the gradient orientations and bvalues.
%B is saved to a .txt file with the same name as the fdt file. (In example2: DTI.txt)
%
%----------------------------------------------
%
%Downloaded from Angelos Barmpoutis' web-page.
%
%This program is freely distributable without 
%licensing fees and is provided without guarantee 
%or warrantee expressed or implied. This program 
%is NOT in the public domain. 
%
%Copyright (c) 2007 Angelos Barmpoutis. 
%
%VERSION : 20110611
%
%-----------------------------------------------
function writeFDT(name,data,bvalue)
sz=size(data);
sz1=sz(1);

if length(sz)>=2
    sz2=sz(2);
else
    sz2=1;
end

if length(sz)>=3
    sz3=sz(3);
else
    sz3=1;
end

if length(sz)>=4
    sz4=sz(4);
else
    sz4=1;
end

fid=fopen(name,'w','b');
sum=1;

fwrite(fid,sz1,'int');
fwrite(fid,sz2,'int');
fwrite(fid,sz3,'int');
fwrite(fid,sz4,'int');


for i=1:sz4
   for j=1:sz3
      for x=1:sz2
          fwrite(fid,data(:,x,j,i),'float');
      end
   end
end

fclose(fid);

if nargin==3
    bname=name;
    bname(length(name))='t';
    bname(length(name)-1)='x';
    bname(length(name)-2)='t';
    write_TXT(bname,bvalue);
end

function write_TXT(fnm,A);

sz=size(A);

if length(sz)~=2
    fprintf(1,'ERROR: Matrix A must be 2D. You gave a %dDmatrix.\n',length(sz));
    return;
end
    
f=fopen(fnm,'wt');
if f~=0
    for i=1:sz(1)
        for j=1:sz(2)
            fprintf(f,'%.8f ',A(i,j));
        end
        fprintf(f,'\n');
    end
    fclose(f);
else
    fprintf(1,'ERROR: The file cannot be opened.\n');
    return;
end