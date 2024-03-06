function p = generateBatchStruct(midx, iM, variable_fields, variables)

[mesh, material, optimization] = standardParameters;

problem_set = {'type','force_type','boundaries'};
problem_variables = intersect(problem_set,variable_fields);
count = numel(problem_variables);

for i=1:length(midx)
    p(i) = generateParameters(mesh, material, optimization);
end
for i=1:length(midx)
    idx = iM(midx(i),:);
    for j=1:length(idx)
        p(i).(variable_fields{j}) = variables{j}(idx(j));
        if strcmp(variable_fields{j},'freq')
            if iscell(variables{j})
                tmp = variables{j}(idx(j));
                p(i).('freq') = tmp{1};
            end
            p(i).('nfreq') = length(p(i).('freq'));
            p(i).('omega') = 2*pi*p(i).('freq');
        end
        if strcmp(variable_fields{j},'proportional')
            damping = p(i).('proportional'){1};
            p(i).('alpha_') = damping(1)*p(i).('freq');
            p(i).('beta_') = damping(2)/p(i).('freq');
        end
        if strcmp(variable_fields{j},'nsub')
            p(i).(variable_fields{j}) = repmat(variables{j}(idx(j)),1,2);
        end
    end
    tf = strcmp(problem_variables, problem_set);
    arg = {"square", [1], [1 2 3 4]};
    arg{tf} = variables{length(idx)-count+1:length(idx)}(idx(end-count+1:end));   
    p(i).problem_data = generateProblem(arg{1}, arg{2}, arg{3});   
end