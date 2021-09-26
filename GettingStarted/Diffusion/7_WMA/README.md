<p align="center">
<a href="https://github.com/delucaal/MRIToolkit"> 
<img src="../../../img/MRIToolkitLogo.png" height="150"/> 
 </a> 
 </p>

# Fiber tractography segmentation with WMA [update 14-03-2020] 
The White Matter Analysis (WMA) package [Zhang et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29920375/) allows to automatically reconstruct major white matter bundles from whole brain fiber tractography, as in the following example:

<p align="center">
<img src="../../../img/ExampleTractsWMA" width="600"/> 
 </p>

# How does it practically work?
WMA is coded in Python and works with VTK files (from Slicer3D). MRIToolkit allows to run the WMA pipeline on whole brain fiber tractography obtained with MRIToolkit, ExploreDTI or mrtrix3 (see below) with a single call:
```matlab
MRTTrack.WMATractographyClustering('mat_file','target_mat_file','tract_file','TractographyInMatFormat.mat','output','output_folder_for_WMA')
```

After some computation time, fiber tractography of the major bundles will be stored in the folder output_folder_for_WMA/AnatomicalTracts_EDTI.

To run the call above, WMA should be installed in your system, as detailed in the next section.

# Installing WMA for MRIToolkit
- Step1: if you do not have Anaconda/Miniconda installed, please download and install it [miniconda3](https://docs.conda.io/en/latest/miniconda.html). There are many guides for this step available online. 
- Step2: in terminal, type the following lines:
```bash
conda create -n wma_mritoolkit python
conda activate wma_mritoolkit
pip install git+https://github.com/delucaal/whitematteranalysis.git
```
As you might notice, I forked the [original repository](https://github.com/SlicerDMRI/whitematteranalysis) to ensure the installed version is compatible with MRIToolkit. 
- Step3: download and unpack/install [Slicer3D v4.10.2](https://slicer-packages.kitware.com/#collection/5f4474d0e1d8c75dfc70547e/folder/60add36eae4540bf6a89be73). Note: the code will likely not run correctly with versions other than 4.10.2.
- Step4: Download and unpack the [WMA atlas](https://www.dropbox.com/s/beju3c0g9jqw5uj/WMA_tutorial_data.zip?dl=0)
- Step5: include the following lines in your MRIToolkitDefineLocalVars.m script and adjust according to your system: 
```matlab
MRIToolkit.SlicerPath = 'path to slicer'; % On MacOS, it will be something like: /Applications/Slicer3D.app/Contents/MacOS
MRIToolkit.WMAAtlasFolder = 'path to the white matter atlas'; % Point to the root folder of the WMA atlas
MRIToolkit.Miniconda3Path = 'path to the bin folder of your miniconda3 installation'; % For example: /opt/bin/miniconda3/bin/
MRIToolkit.Miniconda3Env = 'Name of the Conda3 environment with the WMA package. Empty for the base environment'; % If you followed the steps above, type 'wma_mritoolkit' here
```

# Can I combine this with tractography from mrtrix3?
Yes! Just convert your tractography from tck to MAT format, run WMA and convert the output back to tck, as explained [here](https://github.com/delucaal/MRIToolkit/tree/master/GettingStarted/Diffusion/4_MRTrix_Dipy). I realise command line tools would probably be more efficient for interfacing with mrtrix3 - they are on my todo list. 