function plot_W_FRF(xval)

load('init.mat');

% Probe dof
y_dofs = sp.comp_dofs{2};
b_dofs = sp.boundary(2).dofs;
by_dofs = intersect(y_dofs,b_dofs);
probe_dof = median(by_dofs);


[K, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_, YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, xval);
F = builForce(sp,msh,problem_data);
F = F/sum(F);

W = zeros(50,1);
freq = zeros(50,1);
for i=1:50
    freq(i) = 10*(i-1);
    omega = 2*pi*freq(i);
    Kd = K +1j*omega*C -omega*omega*M;
    u(dr_dofs) = 0;
    F(free_dofs) = F(free_dofs) -Kd(free_dofs,dr_dofs)*u(dr_dofs);
    u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
    W(i) = 0.5*omega*omega*real((u')*C*u); % Active input power
end

    
end