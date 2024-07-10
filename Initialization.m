function DefaultOptions = Initialization
%INITIALIZATION Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 31-Mar-2023 NL

clc
clear
close all

% cd(fileparts(which('migration.m')))
addpath(genpath('mfiles'))
load('DefaultOptions.mat')

% Check Bio-Formats and initialize logging
[status, version] = bfCheckJavaPath();
disp(['Bioformat version: ' version])
bfInitLogging('OFF');
bfUpgradeCheck(true);

try
    pythobj = pyenv('Version','C:\env-cellpose-qupath\Scripts\python.exe','ExecutionMode','InProcess');
end

try
    parpool(20);
end
disp('System ready...')
end

