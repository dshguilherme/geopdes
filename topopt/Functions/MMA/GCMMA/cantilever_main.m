clearvars
clc
close all

%% Initialize Parameters
parameters.degree = 1;
parameters.nsub = [60 20];
parameters.vol_frac = 0.49;
parameters.rmin = 2;
parameters.change_min = 1e-4;
parameters.iter_max = 200;

%Domain and Material properties
parameters.freq = 5;
parameters.omega = parameters.freq*2*pi;
parameters.YOUNG = 210e9;
parameters.YOUNG_MIN = 1e-3;
parameters.RHO = 7860;
parameters.RHO_MIN = 1e-3;
parameters.alpha_ = 0; %1.2*parameters.omega;
parameters.beta_ = 0.1/parameters.omega;

% Optimization Options
%%%% Possible objective_function strings %%%%%
% "compliance" - static compliance
% "scaled compliance" - 1 to 100 compliance, scaled from the initial result
% "dB compliance" - 100 +10*log10(Cs)/100 +10*log10(Cs0)
% AIP - Active Input Power
% dB AIP - Active Inpute Power in dB scaling
% mixed - neta*(dB AIP) + (1-neta)*(dB compliance)
parameters.objective_function = "mixed";
parameters.neta = 0.9; % If objective_function = "mixed", neta = the mixing coefficient

parameters.philter = "density"; % "simple" for sensitivity, "density" for density
initialize_cantilever(parameters);

% Load initialized parameters
load('init.mat');

%% Optimization Strategy
% OC
% xval = OC(f1, 'init.mat', filter_options);

% MMA
% xval = MMA(f1,'init.mat', filter_options); 

% GCMMA
figure(1)
[xval, fobj, fres, x_history] = GCMMA(f1,f2,'init.mat', filter_options);

%% Calc W and Cs from each time step
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


figure(3)
plot(fobj,'LineWidth',2)
hold on
plot(parameters.neta*W_history,'LineWidth',2)
plot((1-parameters.neta)*Cs_history,'LineWidth',2)
set(gca,'FontSize',18)
legend('f_{obj}','W_{dB}','C_s')
title('Convergence','FontSize',20,'FontWeight','bold')
grid on
xlabel('Iteration','FontSize',20,'FontWeight','bold')
ylabel('Objective Function','FontSize',20,'FontWeight','bold')


