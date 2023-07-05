function [freq, W] = plot_W_FRF(x_val)

load('init.mat');

[K, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, 0, 1e-8, YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x_val);

W = zeros(200,1);
freq = zeros(200,1);
u = zeros(sp.ndof,1);
for i=1:200
    freq(i) = 10*(i-1);
    omega = 2*pi*freq(i);
    Kd = K +1j*omega*C -omega*omega*M;
    u(dr_dofs) = 0;
    F(free_dofs) = F(free_dofs) -Kd(free_dofs,dr_dofs)*u(dr_dofs);
    u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
    W(i) = 0.5*omega*omega*real((u')*C*u); % Active input power
end

    
end