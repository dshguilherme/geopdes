function [mesh, material, optimization] = standardParameters
mesh.degree = 2;
mesh.nsub = [134 134];

% Domain and Material properties
material.freq = 350;
material.rho = 2700;
material.young = 6.9e13;
material.POISSON = 0.3;
material.alpha = 1e-5;
material.beta = 1e-7;

% Optimization parameters
optimization.thickness = 1e-3;
optimization.min_thickness = 1e-3/2;
optimization.max_thickness = 1e-3*2;
optimization.max_add = 1;
optimization.max_take = 0.5;


optimization.obj_function = "AIP"; % v2_rms or AIP
optimization.initial_guess = 0; % 1 for initial guess with eigenvalue decompposion, 0 for constant initial guess
end