clearvars
clc
close all

%% Initialization

problem_data = scordelis_problem;

% Mesh parameters
parameters.degree = 3;
parameters.nsub = [30 30];

% Domain and Material properties
parameters.freq = 0;
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 7860;
parameters.YOUNG = 210e9;
parameters.POISSON = 0.3;
parameters.alpha_ = 0;
parameters.beta_ = 0.1/parameters.omega;

% Optimization parameters
parameters.min_thickness = 1e-3;
parameters.max_thickness = 1e-2;

parameters.max_added_mass = 1;
parameters.rmin = 1;
parameters.change_min = 1e-4;
parameters.iter_max = 50;
parameters.philter = "density"; % 'simple' or 'density'
parameters.neta = 0.9;
parameters.objective_function = "compliance";

%% Solve for initial step
initial_step_scordelis(parameters, problem_data);
load('init_scordelis.mat')

%% Optimization
figure(1)
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_scordelis.mat', filter_options);

%% Plots

% History plots
objective_function = 'History';
save('init_scordelis.mat', '-append', 'objective_function');
W_history = zeros(length(fobj),1);
Cs_history = W_history;
for i=1:length(fobj)
    x = x_history(:,i);
    [WC, ~] = f2(x);
    W_history(i) = WC(1);
    Cs_history(i) = WC(2);
end
