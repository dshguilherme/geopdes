prefix = string('L1_134x134_t0010_');
variable_fields = {'freq', 'force_type'};

freqs = string({'f0100', 'f1000', 'f1500'});
freq_num = [100 1000 1500];

forces = string({'angled15', 'angled30'});
force_num = [9 10];

variables = {freq_num, force_num};
input = {freqs, forces};

[all_names, indexMatrix] = generateNameList(prefix, input);
[missing_names, midx] = missingFileNames({'L1'}, all_names);

[mesh, material, optimization] = standardParameters;
parameters = generateParameters(mesh, material, optimization);
problem_data = generateProblem("square", 1, [1 2 3 4]);

for i=1:length(midx)
    idx = indexMatrix(midx,:);
    p = parameters;
    for j=1:length(variable_fields)
        p.(variable_fields{j}) = variables{j}(idx(j));
    end
    filename = missing_names(i);
    batchOptimizeKLShell(problem_data, p, filename);
end

