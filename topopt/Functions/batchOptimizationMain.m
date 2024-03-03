clearvars
close all
clc
%% Processing options
CALC_INITIAL_SPECTRA = true; % True or False
% mesh options
    initial_degree = 2;
    initial_nsub = [60 60];
% Spectral options
    startFreq = 1;
    stopFreq = 1000;
    discretization = 1;
    sAlpha = 1e-5; % Alpha Damping for plotting
    sBeta = 1e-7; % Beta Damping for plotting


%%  Mesh, material, optimization and problem batch variables
mesh_options = []; % 1 - degree, 2 - nsub
domain_options = [1 5]; % 1 - freq, 2 - rho, 3 - young, 4 - POISSON, 
                      %5 - alpha, 6 - beta, 7 - proportional
optimization_options = []; % 1 - thickness, 2 - min_thickness, 3 - max_thickness
                           % 4 - max_add, 5 - max_take, 6 - obj_function
                           % 7 - initial_guess
problem_options = [2]; % 1 - type, 2 - force_type, 3 - boundaries

variable_fields = generateVariables(mesh_options, domain_options, ...
    optimization_options, problem_options);

%% String names and variable definition
% Mesh
mesh_variables = cell(2,2);

% Domain
domain_variables = cell(7,2);
domain_variables{1,1} = [130 350 443 1310]; % optimization frequencies
domain_variables{1,2} = string({'f0130', 'f0350', 'f0443', 'f1310'});
domain_variables{5,1} = [1e-7, 1e-6, 1e-5]; % alpha
domain_variables{5,2} = string({'a-7', 'a-6', 'a-5'});

% Optimization
optimization_variables = cell(7,2);

% Problem
problem_variables = cell(3,2);
problem_variables{2,1} = [1 2];
problem_variables{2,2} =  string({'fCentered', 'fDistributed'});

% variables = {mesh_variables, domain_variables, optimization_variables, problem_variables};

variables = [mesh_variables(:,1); domain_variables(:,1); ...
    optimization_variables(:,1); problem_variables(:,1)];
variables = variables(~cellfun('isempty',variables(:)));

variable_names = [mesh_variables(:,2); domain_variables(:,2); ...
    optimization_variables(:,2); problem_variables(:,2)];
variable_names = variable_names(~cellfun('isempty',variable_names(:)));


%% Generate filenames
prefix = string('L1_134x134_t0010_');
[all_names, indexMatrix] = generateNameList(prefix, variable_names);
[missing_names, midx] = missingFileNames({'L1'}, all_names);

% Generate struct for batch processing
p = generateBatchStruct(midx, indexMatrix, variable_fields, variables);

% Batch processes
parfor i=1:length(midx)
    filename = missing_names(i);
    batchOptimizeKLShell(p(i).problem_data, p(i), filename);
end

%% Save plots of solutions and spectra
if CALC_INITIAL_SPECTRA
    calculateInitialFRFs(initial_degree, initial_nsub, startFreq, stopFreq, ...
        discretization, sAlpha, sBeta, initial_filename)
end
% Open file with the original spectra
clearvars -except all_names problem_variables startFreq stopFreq discretization sAlpha sBeta
initial_spectra = load(initial_filename+string('.mat'));
calculateBatchSpectras(all_names, startFreq, stopFreq, discretization, ...
    sAlpha, sBeta, initial_spectra.frequency_array, initial_spectra.AIP);


