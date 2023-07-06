clearvars
load('init.mat')

diffs = zeros(100,6);
for i=1:100
freq = i*5;
omega = 2*pi*freq;
dx = 0.01;
x_f = xval;
x_f(20) = x_f(20) +dx;
x_b = xval;
x_b(20) = x_b(20) -dx;


[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, xval);
Kd = Ks +1j*omega*C -omega*omega*M;
dr_values = zeros(length(dr_dofs),1);
u = SolveDirichletSystem(Kd,F,dr_dofs,free_dofs,dr_values);
% us = SolveDirichletSystem(Ks,F,dr_dofs,free_dofs,dr_values);

[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x_f);
Kd = Ks +1j*omega*C -omega*omega*M;
u_f = SolveDirichletSystem(Kd,F,dr_dofs,free_dofs,dr_values);
% us_f = SolveDirichletSystem(Ks,F,dr_dofs,free_dofs,dr_values);

[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x_b);
Kd = Ks +1j*omega*C -omega*omega*M;
u_b = SolveDirichletSystem(Kd,F,dr_dofs,free_dofs,dr_values);
% us_b = SolveDirichletSystem(Ks,F,dr_dofs,free_dofs,dr_values);


% Cs0 = F'*us;
% Cs_f = F'*us_f;
% Cs_b = F'*us_b;
% Cs_fd = (Cs_f-Cs0)/dx;
% Cs_bd = (Cs0 -Cs_b)/dx;
% Cs_cd = (Cs_f -Cs_b)/2/dx;
% dcs = StiffnessSensitivities(-us, us, xval, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
% dcs_f = StiffnessSensitivities(-us_f, us_f, x_f, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);
% dcs_b = StiffnessSensitivities(-us_b, us_b, x_b, msh.nel, lm, Ke, YOUNG, YOUNG_MIN);

W0 = 0.5*omega*real(1j*F'*u);
% W0 = real(F'*u);
W_f = 0.5*omega*real(1j*F'*u_f);
% W_f = real(F'*u);
W_b = 0.5*omega*real(1j*F'*u_b);
% W_b = real(F'*u);
W_fd = (W_f-W0)/dx;
W_bd = (W0 -W_b)/dx;
W_cd = (W_f -W_b)/2/dx;
lambda_ = -0.5*1j*omega*u;
% lambda_ = -u;
dkd = DynamicStiffnessSensitivities(lambda_, omega, u, xval, msh.nel, lm, ...
    Ke, Me, YOUNG, YOUNG_MIN, RHO, RHO_MIN, alpha_, beta_);
lambda_ = -0.5*1j*omega*u_f;
% lambda_ = -u_f;
dkd_f = DynamicStiffnessSensitivities(lambda_, omega, u_f, xval, msh.nel, lm, ...
    Ke, Me, YOUNG, YOUNG_MIN, RHO, RHO_MIN, alpha_, beta_);
lambda_ = -0.5*1j*omega*u_b;
% lambda_ = -u_b;
dkd_b = DynamicStiffnessSensitivities(lambda_, omega, u_b, xval, msh.nel, lm, ...
    Ke, Me, YOUNG, YOUNG_MIN, RHO, RHO_MIN, alpha_, beta_);

% diffs(i,1) = dcs(20)-Cs_cd;
% diffs(i,2) = dcs_f(20)-Cs_fd;
% diffs(i,3) = dcs_b(20)-Cs_bd;
diffs(i,4) = dkd(20)-W_cd;
diffs(i,5) = dkd_f(20)-W_fd;
diffs(i,6) = dkd_b(20)-W_bd;


% fprintf('Compliance Sensitivity \n')
% fprintf('Central diff: %3.2e | Forward diff: %3.2e | Backwards diff: %3.2e \n', diffs(i,1),diffs(i,2),diffs(i,3))
% fprintf('Active Power Input Sensitivity \n')
fprintf('Freq: %3.0i | Central diff: %3.2e | Forward diff: %3.2e | Backwards diff: %3.2e \n', freq, diffs(i,4),diffs(i,5),diffs(i,6))

end

