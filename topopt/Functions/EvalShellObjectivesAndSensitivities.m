function [f0val, df0dx, fval, dfdx] = EvalShellObjectivesAndSensitivities(xval)
%% Filtering
t = xval;
load('init_shell.mat');
t = apply_x_filter(filter_options,t);

%% Assembly
Ks = shellStiffnessFromElements(Bke, Ske, lm, t);
M = shellMassFromElements(Me, lm, t);
C = alpha_*M +beta_*Ks;
Kd = Ks +1j*C -omega*omega*M;

%% Solve problems
dr_values = zeros(length(dr_dofs),1);
u = SolveDirichletSystem(Kd, F, dr_dofs, free_dofs, dr_values);
us = SolveDirichletSystem(Ks, F, dr_dofs, free_dofs, dr_values);

%% Objective Functions and Constraints

% Compliance
Cs = F'*us; % Static Compliance
Cs_scaled = 100*Cs/Cs0; % Scaled Compliance

% Mean quadratic velocity
V2_rms = 0;
Area = sum(Ae);
for i=1:msh.nel
    dofs = lm(i,:);
    V2_rms = V2_rms + Ae(i)*(u(dofs)'*u(dofs));
end
V2_rms = (omega*omega/Area)*V2_rms;

V2_scaled = 100*V2_rms/V0; % Scaled quadratic velocity
V2_db = 100 +10*log10(V2_rms); % dB mean quadratic velocity
V0_db = 100 +10*log10(V0); % V0 dB
V2_db_scaled = 100*V2_db/V0_db; % Scaled dB mean quadratic velocity

M_max = RHO*sum(Ae.*thickness)*(1+maximum_to_add); % Maximum added mass
M_min = RHO*sum(Ae.*thickness)*(1-maximum_to_take); % Minimum added mass
Mass = RHO*sum(Ae.*t); % Current mass
h1 = Mass -M_max;
h2 = -Mass +M_min;

fval = [h1; h2];

%% Sensitivities

% Adjoint Problem
UR = real(u);
UI = imag(u);
lambda_ = Kd\((-2*omega*omega/Area).*(UR-UI));

% Velocity Sensitivity
m_shape = [sp.nsh_max, sp.nsh_max];
dv2dt = zeros(msh.nel,1);
dcs = dv2dt;
for i=1:msh.nel
    dofs = lm(i,:);
    ell = lambda_(dofs).'/Ae(i);
    you = u(dofs);
    Bend = 3*t(i)*t(i)*reshape(Bke(i,:),m_shape);
    Stress = reshape(Ske(i,:),m_shape);
    dk = Bend+Stress;
    dm = reshape(Me(i,:),m_shape);
    dc = alpha_*dm +beta_*dk;
    dkd = dk+1j*dc-omega*omega*dm;
    dv2dt(i) = real(ell*dkd*you);
    dcs(i) = -us(dofs)'*dk*us(dofs);
end

% Restrictions
dh1dt = RHO.*Ae';
dh2dt = -RHO.*Ae';

% Chain Rules
dc_scaled = 100*dcs/Cs0;

dv_scaled = 100*dv2dt/V0;
dv_db_constant = 10/(log(10)*V2_rms);
dv_db = 100*dv_db_constant*dv2dt/V0_db;

%% Output
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function

switch objective_function
    case "compliance"
        f0val = Cs;
        df0dx = dcs;
    case "scaled compliance"
        f0val = Cs_scaled;
        df0dx = dc_scaled;
    case "v2_rms"
        f0val = V2_rms;
        df0dx = dv2dt;
    case "v2_scaled"
        f0val = V2_scaled;
        df0dx = dv_scaled;
    case "v2_db"
        f0val = V2_db_scaled;
        df0dx = dv_db;
    case "mixed"
        f0val = neta*V2_db_scaled +(1-neta)*Cs_scaled;
        df0dx = neta*dv_db +(1-neta)*dc_scaled;
end
dfdx = [dh1dt; dh2dt]; 
end