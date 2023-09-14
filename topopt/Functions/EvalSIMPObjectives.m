function [f0val, fval] = EvalSIMPObjectives(xval)
%% Filtering
x = xval;
load('init.mat');
% x = apply_x_filter(filter_options,x);

%% Assembly
[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x);
% [Ks, C, M] = SIMPMatrices_WQ(msh, sp, geometry, YOUNG, POISSON, RHO, alpha_, beta_, reshape(x,filter_options.subshape));
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

%% Output
switch objective_function
    case "compliance"
        f0val = Cs;
    case "scaled compliance"
        f0val = Cs_scaled;
    case "dB compliance"
        f0val = Cs_db_scaled;
    case "AIP"
        f0val = aW;
    case "dB AIP"
        f0val = aW_db_scaled;
    case "mixed"
        f0val = neta*aW_db_scaled +(1-neta)*Cs_scaled;
    case "History"
        f0val = [aW_db_scaled, Cs_scaled];
    case "Initial"
        f0val = [Cs, aW];
end
fval = Vc;
end