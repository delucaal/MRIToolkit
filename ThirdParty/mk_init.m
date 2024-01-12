%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



% A. De Luca
% Sub-toolbox specific initialization
% 24/09/2018: creation - v1.1
% 12/01/2024: v1.2 Added ShellAndPython.m
global MRIToolkit;
MRIToolkit.thirdparty_version = 1.1;

addpath(get_executed_file_path())
addpath(fullfile(get_executed_file_path(),'EPG'))
addpath(fullfile(get_executed_file_path(),'TCK'))
addpath(fullfile(get_executed_file_path(),'TRK'))
