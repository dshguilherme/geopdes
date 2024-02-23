%% optimizeSquareShell
% A function that optimizes the harmonic problem for a square shell using
% Kirchhoff-Love shell elements from Isogeometric Analysis.
% INPUTS
% degree - degree of the NURBS basis functions
% nsub - [nsub_x nsub_y] array of mesh subdivisions
% freq - frequency of optimization. [f1 f2] for multiopt
% Fmag - total norm of the Force Vector (magnitude)
% force_type - can be 1, 2, 3, 4 or 5
%            - 1- point force, centered
%            - 2- distributed force over all plate
%            - 3- point force, offcenter x=1/4 y=3/4
%            - 4- distributed force, centered on a 0.2x0.2 area
%            - 5- distributed force, off-centered on a 0.2x0.2 area x=1/4 y=3/4
% fixed_sides - array with up to 4 numbers, for each boundary [1 2 3 4]
% initial_thickness - initial mean thickness of the plate
%  max_add - maximum percentage of inital mass to add. Positive value
% max_take - maximum percentage of initial mass to take. Positive value
%  iter_max - maximum number of iterations for the GCMMA algorithm
%  filename - name of the files to save. Will save two files:
%           - filename.txt (Paraview data)
%           - filename.mat (optimization data)

function optimizeSquareShell(degree, nsub, freq, Fmag, force_type, obj_function, ...
    fixed_sides, initial_thickness, max_add, max_take, iter_max, filename)
%% Initialization
problem_data = square_shell_problem(force_type, fixed_sides);

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
parameters.alpha_ = 0.12.*parameters.omega;
parameters.beta_ = 0.01./parameters.omega;
parameters.Fmag = Fmag; % Force magnitude

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
parameters.objective_function = obj_function; % v2_rms or AIP

parameters.initial_guess = 1; % 1 for initial guess with eigenvalue decompposion, 0 for constant initial guess

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

