clearvars
clc
close all

%% Initialization

freqs = [15 60 102 125];
names = {'square_AIP_beta_f15.mat','square_AIP_beta_f60.mat', 'square_AIP_beta_f102.mat', 'square_AIP_beta_f125.mat'};
for ell=4:length(freqs)
    problem_data = square_shell_problem;
    % problem_data = scordelis_problem;
    % problem_data = hemispherical_shell_problem(10, 10);
    % Mesh parameters
    parameters.degree = 2;
    parameters.nsub = [80 80];

    % Domain and Material properties
    parameters.freq = freqs(ell);
    parameters.omega = parameters.freq*2*pi;
    parameters.RHO = 2700;
    parameters.YOUNG = 69e12;
    parameters.POISSON = 0.3;
    parameters.alpha_ = 0; %1.2*parameters.omega;
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
    parameters.objective_function = "AIP";
    %% Solve for initial step
    initial_step_KL_shell(parameters, problem_data);
    load('init_shell.mat')
    %% Optimization
[xval, fobj, fres, x_history] = ...
    GCMMA_no_plot(f1, f2, 'init_shell.mat', filter_options);
clearvars -except freqs ell names xval fobj fres x_history
save(names{ell})
end
