clearvars
close all
clc
%% Processing options
CALC_INITIAL_SPECTRA = true; % True or False
CALC_SPECTRAS = true; % True or False
initial_filename = 'espectros_iniciais_n3p5';
% mesh options
    initial_degree = 5;
    initial_nsub = [66 66];
% Spectral options
    startFreq = 1;
    stopFreq = 1000;
    discretization = 1;
    sAlpha = 1e-5; % Alpha Damping for plotting
    sBeta = 1e-7; % Beta Damping for plotting


%%  Mesh, material, optimization and problem batch variables
mesh_options = []; % 1 - degree, 2 - nsub
domain_options = [1 5 6]; % 1 - freq, 2 - rho, 3 - young, 4 - POISSON, 
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
% mesh_variables{1,1} = [3 4 5];
% mesh_variables{1,2} = string({'p3','p4','p5'});
% mesh_variables{2,1} = [44 66 88 110 132];
% mesh_variables{2,2} = string({'44x44','66x66','88x88', '110x110'});

% Domain
domain_variables = cell(7,2);
domain_variables{1,1} = [130 350 443]; % optimization frequencies
domain_variables{1,2} = string({'f0130', 'f0350', 'f0443'});
domain_variables{5,1} = [0 1e-4 1e-3 1e-2 1e-1]; % alpha
domain_variables{5,2} = string({'a0', 'a-4', 'a-3', 'a-2', 'a-1'});
domain_variables{6,1} = [1e-4 1e-5 1e-6 1e-7 0]; % beta
domain_variables{6,2} = string({'b-4', 'b-5', 'b-6', 'b-7', 'b0'}); 

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
prefix = string('L1_t0010_');
[all_names, indexMatrix] = generateNameList(prefix, variable_names);
[missing_names, midx] = missingFileNames({'L1'}, all_names);

% Generate struct for batch processing
p = generateBatchStruct(midx, indexMatrix, variable_fields, variables);

tic
fprintf('Initiating parallel optimization of %i cases \n', length(missing_names));
% Batch processes
parfor i=1:length(missing_names)
    filename = missing_names(i);
    batchOptimizeKLShell(p(i).problem_data, p(i), filename);
end
timerVal = toc;
hours = fix(timerVal/3600);
minutes = fix(timerVal/60) -hours*60;
seconds = round(rem(timerVal,60));
fprintf('Optimization of %i cases finished in %i hours, %i minutes and %i seconds \n',length(missing_names),hours,minutes,seconds);
%% Save plots of solutions and spectra
if CALC_INITIAL_SPECTRA
    fprintf('Calculating Initial Spectras ... \n');
    calculateInitialFRFs(initial_degree, initial_nsub, startFreq, stopFreq, ...
        discretization, sAlpha, sBeta, initial_filename)
    fprintf('Initial Spectras calculated! \n');
end
% Open file with the original spectra
if CALC_SPECTRAS
    clearvars -except all_names problem_variables startFreq stopFreq discretization sAlpha sBeta initial_filename
    initial_spectra = load(initial_filename+string('.mat'));
    tic
    fprintf('Calculating spectras and plotting results... \n');
    calculateBatchSpectras(all_names, startFreq, stopFreq, discretization, ...
        sAlpha, sBeta, initial_spectra.frequency_array, initial_spectra.AIP, initial_spectra.F);
    timerVal = toc;
    hours = fix(timerVal/3600);
    minutes = fix(timerVal/60) -hours*60;
    seconds = round(rem(timerVal,60));
    if hours > 0
        fprintf('Finished calculating %i spectras in %i hours, %i minutes and %i seconds  \n',length(all_names),hours,minutes,seconds);
    else
        fprintf('Finished calculating %i spectras in  %i minutes and %i seconds  \n',length(all_names),minutes,seconds);
    end
end
