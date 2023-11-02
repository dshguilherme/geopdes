function FRF = CalculateFRF(frequency_array, t_val)
load('init_shell.mat')
Ks = shellStiffnessFromElements(Bke, Ske, lm, t_val, t_val, YOUNG, modo);
M = shellMassFromElements(Me, lm, t_val, t_val, RHO, modo);
C = alpha_*M +beta_*Ks;

FRF = zeros(numel(frequency_array),3);
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
    FRF(i,3) = real(F'*u); % Dynamic Stiffness
end
end