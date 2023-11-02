clearvars
clc
close all

%% Initialize Parameters
% Geometry
L = 1;
h = 0.5;
radius = 1;
problem_data = compliant_mechanism(L,h);

% Mesh Parameters
parameters.degree = 2;
parameters.nsub = [100 50];

% Domain and Material
parameters.freq = 5;
parameters.omega = parameters.freq*2*pi;
parameters.YOUNG = 1; %210e9;
parameters.YOUNG_MIN = 1e-3;
parameters.POISSON = 0.3;
parameters.RHO = 1; %7860;
parameters.RHO_MIN = 1e-3;
parameters.alpha_ = 1.2*parameters.omega;
parameters.beta_ = 0.1/parameters.omega;

% Optimization Parameters
parameters.vol_frac = 0.49;
parameters.rmin = 2;
parameters.change_min = 1e-4;
parameters.iter_max = 200;
parameters.philter = "simple"; % "simple" for sensitivity, "density" for density
parameters.neta = 0.9;
parameters.objective_function = "compliance";
%%%% Possible objective_function strings %%%%%
% "compliance" - static compliance
% "scaled compliance" - 1 to 100 compliance, scaled from the initial result
% "dB compliance" - 100 +10*log10(Cs)/100 +10*log10(Cs0)
% AIP - Active Input Power
% dB AIP - Active Inpute Power in dB scaling
% mixed - neta*(dB AIP) + (1-neta)*(dB compliance)

%% Load/initialize parameters
initialize_problem(parameters, problem_data);
% initialize_biharmonic_problem(parameters, problem_data);
load('init.mat');

%% Optimize
figure(1)
[xval, fobj, fres, x_history] = GCMMA(f1, f2, 'init.mat', filter_options);

%% Plotting

% Calc W and Cs from each time step
objective_function = 'History';
save('init.mat','-append','objective_function');
W_history = zeros(length(fobj),1);
Cs_history = W_history;
for i=1:length(fobj)
    x = x_history(:,i);
    [WC, ~] = f2(x);
    W_history(i) = WC(1);
    Cs_history(i) = WC(2);
end

%% Plots
x_init = vol_frac*ones(size(xval));
[freq, W0] = plot_W_FRF(x_init);
[~, W] = plot_W_FRF(xval);
figure(2)
semilogy(freq,W0,'LineWidth',2)
hold on
grid on
semilogy(freq,W,'LineWidth',2)
xline(parameters.freq,'--', 'LineWidth',2)
set(gca,'FontSize',18)
legend('Initial', 'Min W', 'Target Frequency')
title('Active Input Power','FontWeight','bold','FontSize',20)
xlabel('Frequency [Hz]','FontWeight','bold','FontSize',20)

objective_function = 'eigenvalues';
save('init.mat','-append','objective_function');
eig_history = zeros(length(fobj),10);
for i=1:length(fobj)
eig_history(i,:) = f2(x_history(:,i))
end

% 
% figure(3)
% plot(fobj,'LineWidth',2)
% hold on
% plot(parameters.neta*W_history,'LineWidth',2)
% plot((1-parameters.neta)*Cs_history,'LineWidth',2)
% set(gca,'FontSize',18)
% legend('f_{obj}','W_{dB}','C_s')
% title('Convergence','FontSize',20,'FontWeight','bold')
% grid on
% xlabel('Iteration','FontSize',20,'FontWeight','bold')
% ylabel('Objective Function','FontSize',20,'FontWeight','bold')

