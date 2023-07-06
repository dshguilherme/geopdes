function [f0val, df0dx, fval, dfdx] = EvalSIMPObjectivesAndSensitivities(xval)
%% Filtering
x = xval;
load('init.mat');
x = apply_x_filter(filter_options,x);

%% Assembly
[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x);
Kd = Ks +1j*omega*C -omega*omega*M;

%% Solve problem
dr_values = zeros(length(dr_dofs),1);
u = SolveDirichletSystem(Kd,F,dr_dofs,free_dofs,dr_values);
us = SolveDirichletSystem(Ks,F,dr_dofs,free_dofs,dr_values);

%% Objective Functions and Constraints
Cs = F'*us; % Static Compliance
Cs_scaled = 100*Cs/Cs0; % Scaled Compliance
Cs_db = 100 +10*log10(Cs); % dB Compliance
Cs0_db = 100+10*log10(Cs0); % Cs0 dB
Cs_db_scaled = 100*Cs_db/Cs0_db; % Scaled dB Compliance

aW = 0.5*omega*real(1j*F'*u); % Active Input Power
aW_db = 100 +10*log10(aW); % dB Active Input Power
aW0_db = 100 +10*log10(aW0); % W0 dB
aW_db_scaled = 100*aW_db/aW0_db; % Scaled dB Active Input Power

rW = 0.5*omega*real((u')*(Ks - omega*omega*M)*u); % Reactive Input Power

Ep = real((u')*Ks*u); % Elastic/Potential Energy
Ek = omega*omega*real((u')*M*u); % Kinetic Energy
R = Ep/Ek; % Energy Quotient

Vc = sum(x.*Ve)-vol_frac*sum(Ve); % Volume Constraint

%% Sensitivities
dcs = StiffnessSensitivities(-us, us, x, msh.nel, lm, Ke, YOUNG, YOUNG_MIN); % Static Compliance
% lambda_ = SolveDirichletSystem(Kd,-0.5*1j*omega*conj(F),dr_dofs,free_dofs,dr_values);
% dk = StiffnessSensitivities(lambda_, u, x, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
% dm = MassSensitivities(lambda_,u, x, msh.nel, lm, Me, RHO, RHO_MIN);
% dc = alpha_*dm +beta_*dk;
% dkd = real(dk +1j*omega*dc -omega*omega*dm); % Active Input Power
lambda_ = -0.5*1j*omega*u;
dkd = DynamicStiffnessSensitivities(lambda_, omega, u, x, msh.nel, lm, ...
    Ke, Me, YOUNG, YOUNG_MIN, RHO, RHO_MIN, alpha_, beta_);
dkd = dkd';
dv = Ve'; % Volume Restriction

% Chain-Rules

dc_scaled = 100*dcs/Cs0;
dc_db_constant = 10/(log(10)*Cs);
dc_db = 100*dc_db_constant*dcs/Cs0_db;

dW_db_constant = 10/(log(10)*aW);
dW_db = 100*dW_db_constant*dkd/aW0_db;

%% Output
switch objective_function
    case "compliance"
        f0val = Cs;
        df0dx = dcs;
    case "scaled compliance"
        f0val = Cs_scaled;
        df0dx = dc_scaled;
    case "dB compliance"
        f0val = Cs_db_scaled;
        df0dx = dc_db;
    case "AIP"
        f0val = aW;
        df0dx = dkd;
    case "dB AIP"
        f0val = aW_db_scaled;
        df0dx = dW_db;
    case "mixed"
        f0val = neta*aW_db_scaled +(1-neta)*Cs_scaled;
        df0dx = neta*dW_db +(1-neta)*dc_scaled;
end
fval = Vc;
dfdx = dv;
% Filter chain-rule

[df0dx, dfdx] = apply_sensi_filter(filter_options, x, df0dx, dfdx);
end