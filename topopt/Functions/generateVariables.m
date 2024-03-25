function variable_fields = generateVariables(mesh_options, ...
    domain_options, optimization_options, problem_options)
% Mesh
mesh_fields = {'degree', 'nsub'}; % degree, nsub
mesh_variables = mesh_fields(mesh_options);

% Domain
domain_fields = {'freq', 'rho', 'young', 'POISSON', 'alpha_', ...
    'beta_', 'proportional'}; %freq, rho, young, POISSON, alpha, beta, proportional
domain_variables = domain_fields(domain_options);

% Optimization
optimization_fields = {'thickness', 'min_thickness', 'max_thickness', ...
    'max_add','max_take', 'obj_function', 'initial_guess'};
optimization_variables = optimization_fields(optimization_options);

% Problem batch variables
problem_fields = {'type','force_type', 'boundaries'};
problem_options = [2];
problem_variables = problem_fields(problem_options);

variable_fields = [mesh_variables, domain_variables, ...
    optimization_variables, problem_variables];
end