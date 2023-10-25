clearvars
clc
close all

%% Initialization

% problem_data = square_shell_problem;
problem_data = scordelis_problem;
% Mesh parameters
parameters.degree = 4;
parameters.nsub = [36 36];

% Domain and Material properties
parameters.freq = 50;
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 7860;
parameters.YOUNG = 210e9;
parameters.POISSON = 0.3;
parameters.alpha_ = 1e-8;
parameters.beta_ = 1e-8; %0.1/parameters.omega;

% Optimization parameters
parameters.thickness = 0.25;
parameters.min_thickness = parameters.thickness/2;
parameters.max_thickness = parameters.thickness*2;
parameters.maximum_to_add = .25;
parameters.maximum_to_take = .25;

parameters.rmin = 4;
parameters.change_min = 1e-4;
parameters.iter_max = 150;
parameters.philter = "simple"; % 'simple' or 'density'
parameters.neta = 0.9;
parameters.objective_function = "mixed";
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

% History plots
objective_function = 'History';
save('init_shell.mat', '-append', 'objective_function');
W_history = zeros(length(fobj),1);
Cs_history = W_history;
for i=1:length(fobj)
    x = x_history(:,i);
    [WC, ~] = f2(x);
    W_history(i) = WC(1);
    Cs_history(i) = WC(2);
end