clearvars
clc

%% Initialize Parameters
parameters.degree = 1;
parameters.nsub = [60 20];
parameters.vol_frac = 0.49;
parameters.rmin = 2;
parameters.change_min = 1e-4;
parameters.iter_max = 100;

%Domain and Material properties
parameters.freq = 300;
parameters.omega = parameters.freq*2*pi;
parameters.YOUNG = 210e9;
parameters.YOUNG_MIN = 1e-3;
parameters.RHO = 7860;
parameters.RHO_MIN = 1e-3;
parameters.alpha_ = 0;
parameters.beta_ = 0.1/parameters.omega;

% Optimization Options
%%%% Possible objective_function strings %%%%%
% "compliance" - static compliance
% "scaled compliance" - 1 to 100 compliance, scaled from the initial result
% "dB compliance" - 100 +10*log10(Cs)/100 +10*log10(Cs0)
% AIP - Active Input Power
% dB AIP - Active Inpute Power in dB scaling
% mixed - neta*(dB AIP) + (1-neta)*(dB compliance)
parameters.objective_function = "mixed";
parameters.neta = 0.05; % If objective_function = "mixed", neta = the mixing coefficient

parameters.philter = "density"; % "simple" for sensitivity, "density" for density
initialize_cantilever(parameters);

% Load initialized parameters
load('init.mat');

%% Optimization Strategy
% OC
% xval = OC(f1, 'init.mat', filter_options);

% MMA
% xval = MMA(f1,'init.mat', filter_options); 

% GCMMA
xval = GCMMA(f1,f2,'init.mat', filter_options);

