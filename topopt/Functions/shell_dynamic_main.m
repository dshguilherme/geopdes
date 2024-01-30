clearvars
clc
close all

%% Initialization

problem_data = square_shell_problem;
% problem_data = scordelis_problem;
% problem_data = hemispherical_shell_problem(10, 10);
% Mesh parameters
parameters.degree = 2;
parameters.nsub = [30 30];

% Domain and Material properties
parameters.freq = 1;
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 2700;
parameters.YOUNG = 6.9e13;
parameters.POISSON = 0.3;
parameters.alpha_ = 0;
parameters.beta_ = 0.1/parameters.omega;
parameters.Fmag = 1e6; % Force magnitude

% Optimization parameters
parameters.thickness = 0.1;
parameters.min_thickness = parameters.thickness/2;
parameters.max_thickness = parameters.thickness*2;
parameters.maximum_to_add = .05;
parameters.maximum_to_take = .05;

parameters.rmin = 4;
parameters.change_min = 1e-1;
parameters.iter_max = 44;
parameters.philter = "none"; % 'none', 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.9;
parameters.objective_function = "v2_db";
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function

%% Solve for initial step
initial_step_KL_shell(parameters, problem_data);
load('init_shell.mat')

%% Optimization
figure(1)
[xval, fobj, fres, x_history] = ...
    GCMMA(f1, f2, 'init_shell.mat', filter_options);

%% Modal Extraction
eigen_amount = 10;
frequencies = 1:5:100;

% Reconstruct Ks and M from xval
tval = (tmax-tmin)*xval/100 +tmin;
t = apply_x_filter(filter_options,tval);

Ks2 = shellStiffnessFromElements(Bke, Ske, lm, t, t, YOUNG, modo);
M2 = shellMassFromElements(Me, lm, t, t, RHO, modo);

autovec = zeros(length(Ks2),eigen_amount);
autovec2 = autovec;
[V, W] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),eigen_amount,'smallestabs');
[V2, W2] = eigs(Ks2(free_dofs,free_dofs),M2(free_dofs,free_dofs),eigen_amount,'smallestabs');

autovec(free_dofs,:) = V;
autovec2(free_dofs,:) = V2;

w_n = sqrt(diag(W));
w_n2 = sqrt(diag(W2));
csi_r = alpha_./(2*w_n) +beta_.*w_n/2;
csi_r = zeros(size(csi_r));
csi_r2 = alpha_./(2*w_n2) +beta_.*w_n2/2;
csi_r2 = csi_r;

X = zeros(length(F),length(frequencies));
X2 = zeros(length(F),length(frequencies));

v2rms = zeros(length(frequencies),1);
v2rms2 = zeros(length(frequencies),1);
bigR = blkdiag(R0,R0,R0);
for i=1:length(frequencies)
    omega = 2*pi*frequencies(i);
    damped_nats = w_n.^2 -omega^2 +1i*2*csi_r*omega.*w_n;
    damped_nats2 = w_n2.^2 -omega^2 +1i*2*csi_r2*omega.*w_n2;
    bigW = diag(1./damped_nats);
    bigW2 = diag(1./damped_nats2);
    receptance = autovec*bigW*(autovec.');
    receptance2 = autovec2*bigW2*(autovec2.');
    X(:,i) = receptance*F;
    X2(:,i) = receptance2*F;
    velocity = 1i*omega*X(:,i);
    velocity2 = 1i*omega*X2(:,i);
    v2rms(i) = velocity'*bigR*velocity;
    v2rms2(i) = velocity2'*bigR*velocity2;
end


% FRF Functions
% discretization = 2;
% final_frequency = parameters.freq +50;
% frequencies = linspace(final_frequency-100,final_frequency,(final_frequency)/discretization +1);
% x_init = factor*eeen;
% FRF_init = CalculateFRF(frequencies, x_init);
% FRF_final = CalculateFRF(frequencies, xval);
% % % 
% % Plots
% figure(2)
% ystr = plotFRFs(parameters,frequencies,FRF_init,FRF_final);
% 
% figure(3)
% plotObjectiveAndRestrictions(fobj,fres,ystr)
% clearvars -except xval fobj fres x_history
% save('square_f15_AIP_unmixed.mat')

