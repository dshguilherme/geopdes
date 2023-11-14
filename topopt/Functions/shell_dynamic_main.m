clearvars
clc
close all

%% Initialization

% problem_data = square_shell_problem;
% problem_data = scordelis_problem;
problem_data = hemispherical_shell_problem(10, 10);
% Mesh parameters
parameters.degree = 2;
parameters.nsub = [50 50];

% Domain and Material properties
parameters.freq = 60;
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 0.01;
parameters.YOUNG = 6.825e7;
parameters.POISSON = 0.3;
parameters.alpha_ = 1e-4;
parameters.beta_ = 0; %0.1/parameters.omega;
parameters.Fmag = 1; % Force magnitude

% Optimization parameters
parameters.thickness = 0.1;
parameters.min_thickness = parameters.thickness/2;
parameters.max_thickness = parameters.thickness*2;
parameters.maximum_to_add = .80;
parameters.maximum_to_take = .50;

parameters.rmin = 4;
parameters.change_min = 1e-1;
parameters.iter_max = 50;
parameters.philter = "simple"; % 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.7;
parameters.objective_function = "v2_mix";
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function

%% Solve for initial step
initial_step_KL_shell(parameters, problem_data);
load('init_shell.mat')

%% Optimization
figure(1)
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_shell.mat', filter_options);

%% Plots

% FRF Functions
discretization = 4;
final_frequency = 4000;
frequencies = linspace(0,final_frequency,(final_frequency)/discretization +1);
x_init = factor*eeen;
FRF_init = CalculateFRF(frequencies, x_init);
FRF_final = CalculateFRF(frequencies, xval);

% Plots
figure(2)
ystr = plotFRFs(parameters,frequencies,FRF_init,FRF_final);

figure(3)
plotObjectiveAndRestrictions(fobj,fres,ystr)


