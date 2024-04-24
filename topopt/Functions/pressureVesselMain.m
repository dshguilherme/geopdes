clearvars
close all
clc
%% Processing options
CALC_INITIAL_SPECTRA = true;
initial_filename = 'espectros_iniciais_vessel';
% mesh options
    initial_degree = 2;
    initial_nsub = [1 1];
% Spectral options
    startFreq = 1;
    stopFreq = 1000;
    discretization = 1;
    sAlpha = 1e-5;
    seBeta = 1e-7;

%% Mesh, domain and optimization options
% Mesh
parameters.degree = 2;
parameters.nsub = [2 2];
% Domain
parameters.freq = [100];
parameters.RHO = 7850; % kg/m3
parameters.YOUNG = 180e6; % Pa
parameters.POISSON = 0.33;
parameters.alpha_ = 1e-5; % Rayleigh damping parameters. Mass proportional
parameters.beta_ = 1e-7; % Rayleigh damping parameters. Stiffness proportional
parameters.pressure = 1000; % 1kPa of oscilation
% Optimization
parameters.thickness = 3e-3; % mean thickness of the vessel
parameters.min_thickness = 1e-3; % minimum thickness
parameters.max_thickness = 10e-3; % maximum thickness
parameters.max_add = 10; % Max added mass = max_add*initial_mass
parameters.max_take = 1/3; % Minimum mass = max_take*initial_mass
parameters.objective_function = 'v2_rms';
parameters.initial_guess = false;
parameters.iter_max = 1000;
parameters.change_min = 1e-3;

%% Problem declaration
parameters.A = .47;
parameters.B = 1.775;
problem_data = pressure_vessel_problem(parameters.A,parameters.B);

%% Initial step and optimization
io = firstStepPressureVessel(parameters, problem_data);
[xval, fobj, fres, x_history] = fastGCMMA_noplots(io);

%% Save paraview and Results
grafo = nrbplot(io.geometry.nurbs, [sqrt(length(xval)), sqrt(length(xval))]);
grafo.CData = reshape(xval,sqrt(length(xval)),sqrt(length(xval)));
colorbar
colormap(jet)

%% Calc Spectra
