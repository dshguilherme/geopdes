clearvars
clc


%% Initialize Parameters
% Mesh, filter and algorithm parameters
parameters.degree = 1;
parameters.nsub = [60 20];
parameters.vol_frac = 0.49;
parameters.rmin = 2;
parameters.change_min = 1e-4;
parameters.iter_max = 100;

%Domain and Material properties
parameters.freq = 1;
parameters.omega = parameters.freq*2*pi;
parameters.YOUNG = 1;
parameters.YOUNG_MIN = 1e-3;
parameters.RHO = 1;
parameters.RHO_MIN = 1e-3;
parameters.alpha_ = 0;
parameters.beta_ = 0.1/parameters.omega;

% Optimization Options
parameters.philter = "density"; % "simple" for sensitivity, "density" for density
parameters.objective_function = "compliance";
parameters.neta = eps; % If objective_function = "mixed", neta = the mixing coefficient
initialize_mbb(parameters);

%% Load initialized parameters
load('init.mat');


%% Optimization Strategy
% OC
% xval = OC(f1, 'init.mat', filter_options);

% MMA
% xval = MMA(f1,'init.mat', filter_options);

% GCMMA
xval = GCMMA(f1,f2,'init.mat', filter_options);




