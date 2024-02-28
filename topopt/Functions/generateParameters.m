function parameters = generateParameters(mesh, material, optimization)
% Mesh parameters
parameters.degree = mesh.degree;
parameters.nsub = mesh.nsub;

% Domain and Material properties
parameters.freq = material.freq;
parameters.nfreq = length(parameters.freq);
parameters.omega = parameters.freq*2*pi;
parameters.RHO = material.rho;
parameters.YOUNG = material.young;
parameters.POISSON = 0.3;
parameters.alpha_ = material.alpha;
parameters.beta_ = material.beta;
parameters.Fmag = 1000; % Force magnitude

% Optimization parameters
parameters.thickness = optimization.thickness;
parameters.min_thickness = optimization.min_thickness;
parameters.max_thickness = optimization.max_thickness;
parameters.maximum_to_add = optimization.max_add;
parameters.maximum_to_take = optimization.max_take;

parameters.rmin = 1;
parameters.change_min = 1e-3;
parameters.iter_max = 500;
parameters.philter = "simple"; % 'none', 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.9;
parameters.objective_function = optimization.obj_function; % v2_rms or AIP

parameters.initial_guess = optimization.initial_guess; % 1 for initial guess with eigenvalue decompposion, 0 for constant initial guess
end