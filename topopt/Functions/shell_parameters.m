function [parameters, problem_data] = shell_parameters(problem_handle, degree, nsub, freq, objective_function)
problem_data = problem_handle;
% problem_data = square_shell_problem;
% problem_data = scordelis_problem;
% problem_data = hemispherical_shell_problem(10, 10);
% Mesh parameters
parameters.degree = degree;
parameters.nsub = nsub;

% Domain and Material properties
parameters.freq = freq;
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
parameters.neta = 0.9;
parameters.objective_function = objective_function;
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function
end