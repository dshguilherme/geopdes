clearvars
clc
close all

%% Initialization

freqs = [1350];
names = {'1350_big_R0'};
for ell=1:length(freqs)
% for ell=6:8
    problem_data = square_shell_problem;
    % problem_data = scordelis_problem;
%     problem_data = hemispherical_shell_problem(10, 10);
    % Mesh parameters
    parameters.degree = 2;
    parameters.nsub = [40 40];

    % Domain and Material properties
    parameters.freq = freqs(ell);
    parameters.omega = parameters.freq*2*pi;
    parameters.RHO = 2700;
    parameters.YOUNG = 6.9e13;
    parameters.POISSON = 0.3;
%     aalpha = [1.2*parameters.omega 1.2*parameters.omega 1.2*parameters.omega 1.2*parameters.omega 0 0 0 0 1.2*parameters.omega 1.2*parameters.omega 1.2*parameters.omega 1.2*parameters.omega];
    parameters.alpha_ = 0; %1e-3;
%     bbeta = [0 0 0 0 0.1/parameters.omega 0.1/parameters.omega 0.1/parameters.omega 0.1/parameters.omega 0.1/parameters.omega 0.1/parameters.omega 0.1/parameters.omega 0.1/parameters.omega];
    parameters.beta_ = 0.1/parameters.omega;
    parameters.Fmag = 1e6; % Force magnitude

    % Optimization parameters
    parameters.thickness = 0.1;
    parameters.min_thickness = parameters.thickness/2;
    parameters.max_thickness = parameters.thickness*2;
    parameters.maximum_to_add = .80;
    parameters.maximum_to_take = .50;

    parameters.rmin = 4;
    parameters.change_min = 1e-1;
    parameters.iter_max = 400;
    parameters.philter = "simple"; % 'simple' or 'density'
    parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
    parameters.neta = 0.9;
    parameters.objective_function = "v2_rms";
    %% Solve for initial step
    initial_step_KL_shell(parameters, problem_data);
    load('init_shell.mat')
    %% Optimization
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_shell.mat', filter_options);
clearvars -except freqs ell names xval fobj fres x_history
save(names{ell})
end
