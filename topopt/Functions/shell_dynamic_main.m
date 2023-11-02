clearvars
clc
close all

%% Initialization

problem_data = square_shell_problem;
% problem_data = scordelis_problem;
% problem_data = hemispherical_shell_problem(10, 10);
% Mesh parameters
parameters.degree = 2;
parameters.nsub = [10 10];

% Domain and Material properties
parameters.freq = 50;
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 2700;
parameters.YOUNG = 69e12;
parameters.POISSON = 0.3;
parameters.alpha_ = 1e-4;
parameters.beta_ = 1e-8; %0.1/parameters.omega;
parameters.Fmag = 1e6; % Force magnitude

% Optimization parameters
parameters.thickness = 0.1;
parameters.min_thickness = parameters.thickness/5;
parameters.max_thickness = parameters.thickness*5;
parameters.maximum_to_add = .5;
parameters.maximum_to_take = .49;

parameters.rmin = 4;
parameters.change_min = 1e-4;
parameters.iter_max = 150;
parameters.philter = "simple"; % 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.9;
parameters.objective_function = "v2_rms";
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function

%% Solve for initial step
initial_step_KL_shell(parameters, problem_data);
load('init_shell.mat')

%% Optimization
figure(1)
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_shell.mat', filter_options);

%% Plots

% History plots
objective_function = 'History';
save('init_shell.mat', '-append', 'objective_function');
%        f0val = [V2_rms, aW, Cd, vals(1), vals(2)-vals(1)];

% FRF Functions
discretization = 1;
final_frequency = 200;
frequencies = linspace(0,final_frequency,(final_frequency)/discretization +1);
x_init = thickness*eeen;
FRF_init = CalculateFRF(frequencies, x_init);
FRF_final = CalculateFRF(frequencies, xval);


% Plots
figure(2)
switch parameters.objective_function
     case "v2_rms"
         idx = 1;
         ystr = 'Quadratic Mean Velocity [m/s]';
    case "AIP"
        idx = 2;
        ystr = 'Active Input Power [W]';
    case "dynamic compliance"
        idx = 3;
        ystr = 'Dynamic Compliance [m/N]';
end
    

loglog(frequencies,FRF_init(:,idx),'b','LineWidth',2)
hold on
loglog(frequencies,FRF_final(:,idx),'g','LineWidth',2)
grid on
xline(parameters.freq,'--k',{'Excitation Frequency'},'FontWeight','bold','FontSize',18)
legend('Initial Design','Optimized Design')
set(gca,'FontSize',20)
xlabel('Frequency [Hz]','FontSize',24,'FontWeight','bold')
ylabel(ystr,'FontSize',24,'FontWeight','bold')
title('Initial vs Optimized Design','FontWeight','bold','FontSize',24)

figure(3)
    hold on
    subplot(1,2,1)
    semilogy(fobj,'LineWidth',2)
    grid on
    set(gca,'FontSize',20)
    xlabel('Iteration','FontSize',24,'FontWeight','bold')
    title('Objective Function History','FontWeight','bold','FontSize',24)
    ylabel(ystr,'FontSize',24,'FontWeight','bold')
    subplot(1,2,2)
    hold on
    semilogy(fres(:,1),'LineWidth',2)
    semilogy(fres(:,2),'LineWidth',2)
    grid on
    set(gca,'FontSize',20)
    legend('Restriction 1','Restriction 2')
    xlabel('Iteration','FontSize',24,'FontWeight','bold')
    title('Restriction Function History','FontWeight','bold','FontSize',24)
    ylabel('Mass Restriction','FontSize',24,'FontWeight','bold')


