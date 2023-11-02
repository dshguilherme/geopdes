function [f0val, df0dx, fval, dfdx] = EvalSIMPObjectivesAndSensitivities(xval)
%% Filtering
x = xval;
load('init.mat');
x = apply_x_filter(filter_options,x);

%% Assembly
[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x);
% [Ks2, C, M2] = SIMPMatrices_WQ(msh, sp, geometry, YOUNG, POISSON, RHO, alpha_, beta_, reshape(x,filter_options.subshape));
Kd = Ks +1j*omega*C -omega*omega*M;


%% Solve problem
dr_values = zeros(length(dr_dofs),1);
u = SolveDirichletSystem(Kd,F,dr_dofs,free_dofs,dr_values);
us = SolveDirichletSystem(Ks,F,dr_dofs,free_dofs,dr_values);
% us2 = SolveDirichletSystem(Ks2,F,dr_dofs,free_dofs,dr_values);

%% Objective Functions and Constraints
% Cs = F'*us; % Static Compliance
Cs = LL'*us; % Static Compliance force inverter

Cs_scaled = 100*Cs/Cs0; % Scaled Compliance
Cs_db = 100 +10*log10(Cs); % dB Compliance
Cs0_db = 100+10*log10(Cs0); % Cs0 dB
Cs_db_scaled = 100*Cs_db/Cs0_db; % Scaled dB Compliance

% aW = 0.5*omega*real(1j*F'*u); % Active Input Power
aW = 0.5*omega*real(1j*LL'*u); % Active Input Power force inverter
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
% dcs = fastStiffnessSentivities(ke,lm,-us,us,x,YOUNG); % Static Compliance sensitivity
% lambda_ = SolveDirichletSystem(Kd,-0.5*1j*omega*conj(F),dr_dofs,free_dofs,dr_values);
% dk = StiffnessSensitivities(lambda_, u, x, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
% dm = MassSensitivities(lambda_,u, x, msh.nel, lm, Me, RHO, RHO_MIN);
% dc = alpha_*dm +beta_*dk;
% dkd = real(dk +1j*omega*dc -omega*omega*dm); % Active Input Power
% lambda_ = -0.5*1j*omega*u; % cantilever beam
lambda_ = Kd\(-0.5*1j*omega*LL); % force inverter
dkd = DynamicStiffnessSensitivities(lambda_, omega, u, x, msh.nel, lm, ...
    Ke, Me, YOUNG, YOUNG_MIN, RHO, RHO_MIN, alpha_, beta_);
% dkd = fastDynamicStiffnessSensitivities(ke,me,lm,lambda_,u,omega,x,YOUNG,RHO,alpha_,beta_);
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
    case "dynamic compliance"
        f0val = real(F'*u);
        df0dx = dkd;
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
    case "eigenmax"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),1,'sm');
        vals = diag(vals);
        new_vec = zeros(length(free_dofs)+length(dr_dofs),1);
        new_vec(free_dofs,1) = vec(:,1);
        vec1 = new_vec(:,1);
        f0val = -vals(1);
        dke1 = StiffnessSensitivities(vec1, vec1, x, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
        dme1 = MassSensitivities(vec1, vec1, x, msh.nel, lm, Me, RHO, RHO_MIN);
        df0dx = -dke1 +vals(1)*dme1;
    case "eigengap"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        new_vec = zeros(length(free_dofs)+length(dr_dofs),2);
        new_vec(free_dofs,1) = vec(:,1);
        new_vec(free_dofs,2) = vec(:,2);
        vec1 = new_vec(:,1);
        vec2 = new_vec(:,2);
        f0val = vals(1) - vals(2);
        dke1 = StiffnessSensitivities(vec1, vec1, x, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
        dke2 = StiffnessSensitivities(vec2, vec2, x, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
        dme1 = MassSensitivities(vec1, vec1, x, msh.nel, lm, Me, RHO, RHO_MIN);
        dme2 = MassSensitivities(vec2, vec2, x, msh.nel, lm, Me, RHO, RHO_MIN);
        df0dx = (dke1 -vals(1)*dme1) -(dke2 -vals(2)*dme2);
    case "eigenvalues"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),10,'sm');
        vals = diag(vals);
        f0val = vals;
end
fval = Vc;
dfdx = dv;
% Filter chain-rule

[df0dx, dfdx] = apply_sensi_filter(filter_options, x, df0dx, dfdx);
end