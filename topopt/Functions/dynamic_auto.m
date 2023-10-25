function [fobj, xval, fres, x_history] =  dynamic_auto(degree, nsub, freq, rmin, objective, neta)
%% Initialize Parameters
% Geometry
L = 1;
h = 0.5;
radius = 1;
problem_data = cantilever_beam(L,h);

% Mesh parameters
parameters.degree = degree;
parameters.nsub = nsub;

% Domain and Material
parameters.freq = freq;
parameters.omega = parameters.freq*2*pi;
parameters.YOUNG = 210e9;
parameters.YOUNG_MIN = 1e-3;
parameters.POISSON = 0.3;
parameters.RHO = 7860;
parameters.RHO_MIN = 1e-3;
parameters.alpha_ = 1.2*parameters.omega;
parameters.beta_ = 0.1/parameters.omega;

% Optimization Parameters
parameters.vol_frac = 0.49;
parameters.rmin = rmin;
parameters.change_min = 1e-4;
parameters.iter_max = 50;
parameters.philter = "density"; % "simple" for sensitivity, "density" for density
parameters.neta = neta;
parameters.objective_function = objective;

%% Load/initialize parameters
initialize_problem(parameters, problem_data);
load('init.mat');
[xval, fobj, fres, x_history] = GCMMA(f1, f2, 'init.mat', filter_options);
end