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
parameters.freq = 500;
parameters.omega = parameters.freq*2*pi;
parameters.RHO = 2700;
parameters.YOUNG = 6.9e13;
parameters.POISSON = 0.3;
parameters.alpha_ = 1.2*parameters.omega;
parameters.beta_ = 0.1/parameters.omega;
parameters.Fmag = 1e6; % Force magnitude

% Optimization parameters
parameters.thickness = 0.002;
parameters.min_thickness = parameters.thickness/2;
parameters.max_thickness = parameters.thickness*2;
parameters.maximum_to_add = .05;
parameters.maximum_to_take = .05;

parameters.rmin = 1;
parameters.change_min = 1e-1;
parameters.iter_max = 200;
parameters.philter = "simple"; % 'none', 'simple' or 'density'
parameters.modo = "Continuous"; % 'SIMP' or 'Continuous'
parameters.neta = 0.9;
parameters.objective_function = "v2_rms";

%% Solve for initial step
io = firstStepKLShell(parameters,problem_data);

%% Optimize
figure(1)
[xval, fobj, fres, x_history] = fastGCMMA(io);

%% FRF
frequencies = 1:1:1000;

tmax = io.tmax;
tmin = io.tmin;
tval = (tmax-tmin)*xval/100 +tmin;
t = apply_x_filter(io.filter_options,tval);
t_init = parameters.thickness*ones(size(t));

Ks = shellStiffnessFromElements(io.Bke, io.Ske, io.lm, t, t, parameters.YOUNG, parameters.modo);
Ksi = shellStiffnessFromElements(io.Bke, io.Ske, io.lm, t_init, t_init, parameters.YOUNG, parameters.modo);
M = shellMassFromElements(io.Me, io.lm, t, t, parameters.RHO, parameters.modo);
Mi = shellMassFromElements(io.Me, io.lm, t_init, t_init, parameters.RHO, parameters.modo);

aalpha = 1e-7;
bbeta = 1e-8;
C = aalpha*M +bbeta*Ks;
Ci = aalpha*Mi +bbeta*Ksi;

for i=1:length(frequencies)
    omega = 2*pi*frequencies(i);
    Kd = Ks -omega*omega*M +1j*omega*C;
    Kdi = Ksi -omega*omega*Mi +1j*omega*Ci;
    u = SolveDirichletSystem(Kd, io.F, io.dr_dofs, io.free_dofs, io.dr_values);
    ui = SolveDirichletSystem(Kdi, io.F, io.dr_dofs, io.free_dofs, io.dr_values);
    velocity = -1j*omega*u;
    veli = -1j*omega*ui;
    v2rms(i) = real(velocity'*io.R0*velocity);
    v2rmsi(i) = real(veli'*io.R0*veli);
end


figure(2)
semilogy(frequencies,v2rms)
hold on
semilogy(frequencies,v2rmsi)
xline(parameters.freq,'-k','Optimization Frequency');
legend('Optimized shell','Initial shell')
%% Modal Extraction
% eigen_amount = 100;
% dofs = io.free_dofs;
% [V W] = eigs(Ks(dofs,dofs),M(dofs,dofs),eigen_amount,'smallestabs');
% eigenV = zeros(length(Ks),eigen_amount);
% eigenV(dofs,:) = V;
% V = eigenV;
% Mme = V'*M*V; Km = V'*Ks*V; Cm = V'*C*V; %Modal matrices
% zeta = Cm/(2*sqrt(Km*Mme)); % Modal damping ratio
% 
% % Modal response
% v2rms = zeros(length(frequencies),1);
% for i=1:length(frequencies)
%     omega = 2*pi*frequencies(i);
%     H = inv(Km -omega*omega*Mme +1j*omega*Cm);
%     X = (V*H*V')*io.F;
%     v2 = 1j*omega*X;
%     v2rms(i) = v2'*v2;
% end



