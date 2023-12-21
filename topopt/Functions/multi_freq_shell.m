clearvars
clc
close all
%% Initialization

freqs = [168 327 547];
names = {'multi_big_R0'};
problem_data = square_shell_problem;
parameters.degree = 2;
parameters.nsub = [40 40];
parameters.freq = freqs;
parameters.omega = 2*pi.*parameters.freq;
parameters.RHO = 2700;
parameters.YOUNG = 6.9e13;
parameters.POISSON = 0.3;
parameters.alpha_ = @(omega) 0;
parameters.beta_ = @(omega) 1e-3;
parameters.Fmag = 1e6;

% Optimization parameters
parameters.thickness = 0.1;
parameters.min_thickness = parameters.thickness/2;
parameters.max_thickness = parameters.thickness*2;
parameters.maximum_to_add = .05;
parameters.maximum_to_take = .05;

parameters.rmin = 4;
parameters.change_min = 1e-1;
parameters.iter_max = 400;
parameters.philter = "simple"; % 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.9;
parameters.objective_function = "v2_rms";

%% Solve for initial step
    initial_step_KL_shell_multi(parameters, problem_data);
    load('init_shell_multi.mat')
%% Optimization
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_shell.mat', filter_options);
clearvars -except freqs ell names xval fobj fres x_history
save(names{1})