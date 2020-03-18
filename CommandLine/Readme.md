<a href="https://github.com/delucaal/MRIToolkit"> 
<span style='align:center'> <img src="img/MRIToolkitLogo.png" style="width:200px;  display: block;  margin-left: auto;  margin-right: auto;"/> </span>
 </a> 
 
 # MRIToolkit - Command Line Tools [update 14-03-2020] 

- This is a collection of utilities that will be distributed as compiled executables.  
- They are designed to work with Niftis, so Nifti in - Nifti out, to maximize interoperability

## Tools description
- **mrtd_coordsys_fix**: 
```matlab
This tools can be used to check whether spatial dimension flip / permute are needed to properly process the given diffusion MRI data

usage: mrtd_coordsys_fix -nii file.nii -bval file.bval -bvec file.bvec -out screenshot_file_no_extension (other_options)

-perm: how to permute the spatial dimensions,e.g. "1 2 3" or "2 1 3" etc.
-flip: how to flip the sign of the diffusion gradients, e.g. "0 0 0" or "0 1 0" etc.
```
- **mrtd_test_gradients**: 
```matlab
This tools can be used to check whether gradients flip / permute are needed to properly process the given diffusion MRI data

usage: mrtd_test_gradients -nii file.nii -bval file.bval -bvec file.bvec -out screenshot_file_no_extension (other_options)

-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x] 4=[x z y] =[y z x] 6=[z x y]
-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]
```
- **mrtd_moco_epi**: 
```matlab
This tool performs motion / eddy currents / EPI correction of diffusion data.

usage: mrtd_moco_epi -nii file.nii -bval file.bval -bvec file.bvec -out corrected_file (other_options)

-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x] 4=[x z y] =[y z x] 6=[z x y]
-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]
-epi: .nii file to perform EPI correction (T1 or T2 image)
-epi_constraint: allow deformation on specific axis. Input between quotas "1 1 1"
-epi_reg_mode: image to use for registration to structural. One between "fa" (default), "b0", "dwis"
-epi_normcorr: 0 or 1. Use normalized correlation in place of mutual information to drive the registration.
```
- **mrtd_deconv_fod**: 
```matlab
This tool computed a fiber orientation distribution (FOD) with spherical deconvolution methods.

usage: mrtd_deconv_fod -method chosen_method -nii file.nii -bval file.bval -bvec file.bvec -out corrected_file.nii (other_options)

-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x] 4=[x z y] =[y z x] 6=[z x y]
-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]
-method: one in "csd","mscsd","grl","mfod"
-t1_seg: T1 segmentation file derived from "fsl_anat" (mscsd only)
-lmax: spherical harmonic' order

GRL - mFOD specific options

-aniso_rf_dti: Add an anisotropic RF using the DTI model. Specify the eigenvalues in ms2/um2 as "1.7 0.2 0.2"
-aniso_rf_dki: Add an anisotropic RF using the DKI model. Specify the eigenvalues in ms2/um2 and the mean kurtosis as "1.7 0.2 0.2 1" or "auto" to estimate from the data
-aniso_rf_noddi: Add an anisotropic RF using the NODDI model. Specify 6 parameters  "intra-cellular-volume free-diffusivity*10^9 watson-concentration isotropic-fraction isotropic-diffusivity b0-amplitude"
-iso_rf: Add an isotropic RF using the ADC model. Specify the diffusivity in um2/ms, as 0.7 for GM and 3 for CSF
-shell_weight: Inner-shell weighting factor (GRL only)

GRL default usage: mrtd_deconv_fod -method grl -nii file.nii -bval file.bval -bvec file.bvec -out corrected_file.nii -aniso_rf_dti "2.1 0 0" -iso_rf 0.7 -iso_rf 3
mFOD default usage: mrtd_deconv_fod -method mfod -nii file.nii -bval file.bval -bvec file.bvec -out corrected_file.nii -aniso_rf_dki "auto" -aniso_rf_noddi "0.4 1.7 1 0 3 1" -iso_rf 3
```