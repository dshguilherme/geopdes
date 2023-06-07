m = 1; % Number of Restrictions
n = msh.nel; % Number of variables
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = vol_frac*eeen;
xold1   = xval;
xold2   = xval;
xmin    = 0*eeen;
xmax    = eeen;
low     = xmin;
upp     = xmax;
c       = 1000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
raa0    = 0.01;
raa     = 0.01*eeem;
raa0eps = 0.000001;
raaeps  = 0.000001*eeem;
outeriter = 0;
maxoutit  = 1;
kkttol  = 0;

[free_dofs, dr_dofs] = grab_cantilever_dofs(sp);
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
M = K;
x = xval;
for e=1:msh.nel
    k_e = (1e-3 +(x(e)^3)*(YOUNG-1e-3))*squeeze(Ke(e,:,:));
    if x(e) > 0.1
        m_e = (1e-3 +(x(e))*(RHO -1e-3))*squeeze(Me(e,:,:));
    elseif x(e) <= 0.1
       m_e = (1e-3 +(x(e)^9)*(RHO -1e-3))*squeeze(Me(e,:,:));
    end
        idx = lm(e,:)';
        K(idx,idx) = K(idx,idx) +k_e;
        M(idx,idx) = M(idx,idx) +m_e;
end
K = sparse(K); M = sparse(M);
C = alpha*M +beta*K;
% Solving
Kd = K+1i*omega*C-(omega^2)*M;
u(dr_dofs) = 0;
u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
W0 = abs(F'*u);
W0 = 100 +10*log10(W0);
