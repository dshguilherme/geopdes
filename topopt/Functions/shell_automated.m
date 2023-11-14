clearvars
close all
clc

%% Simulation parameters
degree = [2 3 4];
nsub = {[36 36] [50 50] [80 80]};
freq = [60 140 155 183] ;
objective_function = "v2_rms";
X_hist = cell(3,3,4);
F_res = X_hist;
F_obj = X_hist;

for i=1:3
    for j=1:3
        for k=1:4
[parameters, problem_data] = ...
    shell_parameters(square_shell_problem, degree(i), nsub{j}, freq(k), objective_function);

%% Initial Step and Load
initial_step_KL_shell(parameters, problem_data);
load('init_shell.mat')


%% GCMMA solve
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_shell.mat', filter_options);
        X_hist{i,j,k} = x_history;
        F_res{i,j,k} = fres;
        F_obj{i,j,k} = fobj;
        end
    end
end

%% Save results
clearvars -except X_hist F_res F_obj
save('square_shell_GCMMA_v2.mat');