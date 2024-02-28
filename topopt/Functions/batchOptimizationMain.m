%%  Mesh, material, optimization and problem batch variables
mesh_options = []; % 1 - degree, 2 - nsub
domain_options = [1]; % 1 - freq, 2 - rho, 3 - young, 4 - POISSON, 
                      %5 - alpha, 6 - beta
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
domain_variables = cell(6,2);
domain_variables{1,1} = [130 350 443 1310]; % optimization frequencies
domain_variables{1,2} = string({'f0130', 'f0350', 'f0443', 'f1310'});

% Optimization
optimization_variables = cell(7,2);

% Problem
problem_variables = cell(3,2);
problem_variables{2,1} = 6:11;
problem_variables{2,2} =  string({'FLine','FSenoidalLine','FRandom', ...
    'FAngled15','FAngled30','FAngled45'});

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

% Batch process
for i=1:length(midx)
    filename = missing_names(i);
    batchOptimizeKLShell(p(i).problem_data, p(i), filename);
end

