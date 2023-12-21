function FRF = CalculateFRF(frequency_array, t_val)
tval = t_val;
load('init_shell.mat')
tval = (tmax-tmin)*tval/100 +tmin;
t_val = apply_x_filter(filter_options,tval);
Ks = shellStiffnessFromElements(Bke, Ske, lm, t_val, t_val, YOUNG, modo);
M = shellMassFromElements(Me, lm, t_val, t_val, RHO, modo);
alpha_ = 1e-5;
beta_ = 1e-4;
C = alpha_*M +beta_*Ks;

FRF = zeros(numel(frequency_array),4);
u = zeros(sp.ndof,1);
for i=1:length(frequency_array)
    f = frequency_array(i);
    omega = 2*pi*f;
    Kd = Ks +1j*omega*C -omega*omega*M;
    u(dr_dofs) = 0;
    F(free_dofs) = F(free_dofs) - Kd(free_dofs,dr_dofs)*u(dr_dofs);
    u(free_dofs) = Kd(free_dofs, free_dofs)\F(free_dofs);
    velocity = 1j*omega*u;
    FRF(i,1) = velocity'*velocity; % (V2 rms)
    FRF(i,2) = 0.5*omega*omega*real((u')*C*u); % Active Input Power
    FRF(i,3) = abs(F'*u); % Dynamic Stiffness
    FRF(i,4) = velocity'*blkdiag(R0,R0,R0)*velocity;
end
end