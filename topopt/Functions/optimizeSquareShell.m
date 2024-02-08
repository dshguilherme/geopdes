function optimizeSquareShell(degree, nsub, freq, ...
    initial_thickness, max_add, max_take, iter_max, filename)
%% Initialization
problem_data = square_shell_problem;

% Mesh parameters
parameters.degree = degree;
parameters.nsub = nsub;

% Domain and Material properties
parameters.freq = freq;
parameters.nfreq = length(parameters.freq);
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 2700;
parameters.YOUNG = 6.9e13;
parameters.POISSON = 0.3;
parameters.alpha_ = 1.2.*parameters.omega;
parameters.beta_ = 0.1./parameters.omega;
parameters.Fmag = 1e6; % Force magnitude

% Optimization parameters
parameters.thickness = initial_thickness;
parameters.min_thickness = parameters.thickness/2;
parameters.max_thickness = parameters.thickness*2;
parameters.maximum_to_add = max_add;
parameters.maximum_to_take = max_take;

parameters.rmin = 1;
parameters.change_min = 1e-3;
parameters.iter_max = iter_max;
parameters.philter = "simple"; % 'none', 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.9;
parameters.objective_function = "v2_rms";

%% Solve initial step
io = firstStepParallelKLShell(parameters,problem_data);

%% Optimize
[xval, fobj, fres, x_history] = fastGCMMA_noplots(io);

%% Save Paraview and Results
grafo = nrbplot(io.geometry.nurbs,io.nsub); view(0,90);
[connectivity, coordinates, element_vals, ~] = makeParaviewData(xval, io.nsub, grafo);
close all
save(strcat(filename,'.txt'),'connectivity','coordinates','element_vals','-ascii');
save(strcat(filename,'.mat'),'xval','x_history','fobj','fres');
end

