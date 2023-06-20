function [f0val, df0dx, fval, dfdx] = mbb_functions_and_sensitivities(xval)

%% Filtering
x = xval;
load('init.mat')
x = apply_x_filter(filter_options, x);

%% Assembly
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
for e=1:msh.nel
    k_e = (Emin +(x(e)^3)*(E - Emin))*squeeze(Ke(e,:,:));
    idx = lm(e,:)';
    K(idx,idx) = K(idx,idx) +k_e;
end
K = sparse(K);

%% Solution
u(dr_dofs) = 0;
F(free_dofs) = F(free_dofs) -K(free_dofs, dr_dofs)*u(dr_dofs);
u(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);

%% Objective Functions and Constraints

f0val = F'*u; % Static Compliance
fval = sum(x.*Ve)-vol_frac*sum(Ve); % Volume Constraint


%% Sensitivities
df0dx = zeros(msh.nel,1);

for e=1:msh.nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    df0dx(e) = -3*(x(e)^2)*(E-Emin)*(u(dofs)')*k_e*u(dofs);
end
dfdx = Ve';

[df0dx, dfdx] = apply_sensi_filter(filter_options, x, df0dx, dfdx);
end