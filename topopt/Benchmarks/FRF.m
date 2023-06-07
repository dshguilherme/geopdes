clearvars
clc
%%
degree = 3;
nsub = [60 20];
vol_frac = 0.49;
rmin = 2;
change_max = 1e-4;
max_iter = 1000;
freq = 1;
omega = 2*pi*freq;
csi = 1.2;
kappa = 6; % Penalty of the Static Compliance
restrictions = [1 2]; % Which restrictions to apply. 1=Vol, 2=Static Compliance, 3=R_db 
mm = numel(restrictions);

%% Material data
L = 1;
hh = 0.5;
V0 = vol_frac*hh*L;
YOUNG = 210e9;
RHO = 7860;
alpha = 0; 
beta = 1e-8;
c0 = 100;

%%
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; Emin = 1e-3; rho0=1; rhomin=1e-3;
problem_data.E = E;
problem_data.Emin = Emin;
nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @forceCantileverCentered;
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
% problem_data.rho = @(x, y) rho0*ones(2, size (x, 1), size (x, 2));
problem_data.rho = @(x, y) rho0*ones(size(x));
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data,method_data);

% Pre-calculate element volume, stiffness, mass and dofs
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
M = K;
Ke = zeros(msh.nel,sp.nsh_max,sp.nsh_max);
Me = Ke;
lm = zeros(msh.nel,sp.nsh_max);
for i=1:msh.nel
    [k, m, dofs] = elementaryMatrices(i, sp, sp, msh, ...
        problem_data.lambda_lame, problem_data.mu_lame, problem_data.rho);
    Ke(i,:,:) = 0.49*YOUNG*k;
    Me(i,:,:) = 0.49*RHO*m;
    lm(i,:) = dofs;
end

for e=1:msh.nel
    k_e = squeeze(Ke(e,:,:));
    m_e = squeeze(Me(e,:,:));
    idx = lm(e,:)';
    K(idx,idx) = K(idx,idx) +k_e;
    M(idx,idx) = M(idx,idx) +m_e;
end
F = buildForce(sp,msh,problem_data);
F = F/sum(F);

C = alpha*M +beta*K;


% Probe dof
y_dofs = sp.comp_dofs{2};
b_dofs = sp.boundary(2).dofs;
by_dofs = intersect(y_dofs,b_dofs);
probe_dof = median(by_dofs);

% Dirichlet (clamp)
drch_dofs = sp.boundary(1).dofs; % Index of the dofs of side 1
free_dofs = setdiff(1:sp.ndof, drch_dofs);
u_drch = zeros(numel(drch_dofs),1);

% Dirichlet Boundary Conditions
u(drch_dofs) = u_drch;

% % Applying lifting
% F(free_dofs) = F(free_dofs) -Kd(free_dofs,drch_dofs)*u_drch;

% Solving
H = zeros(100);
freq = zeros(100);
for i=1:50
    freq(i) = 10*(i-1);
    omega = 2*pi*freq(i);
    Kd = K +1*omega*1j*C -omega*omega*M;
    u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs); % Dynamic
    H(i) = u(probe_dof);
end

mag = omega.*abs(H);
phase = atan2(imag(H), real(H));

figure(1)
subplot(2,1,1)
semilogy(freq,mag)
subplot(2,1,2)
plot(freq,phase)
ylim([-pi pi])






