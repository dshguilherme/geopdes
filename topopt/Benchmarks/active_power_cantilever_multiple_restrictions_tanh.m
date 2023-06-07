clearvars
clc
%%
% 1st mode = 370 1st antires = 1110 2nde mode = 1400

degree = 1;
nsub = [75 25];
vol_frac = 0.49;
rmin = 2;
change_max = 1e-4;
max_iter = 1000;
freq = 2.5;
omega = 2*pi*freq;
update = 15;
kappa = 20; % Penalty of the Static Compliance From 1 to 40
restrictions = [1 2]; % Which restrictions to apply. 1=Vol, 2=Static Compliance, 3=R_db 
mm = numel(restrictions);
filters = 2; % 1 for linear density, 2 for heaviside + linear density

%% Material data
L = 1;
hh = 0.5;
V0 = vol_frac*hh*L;
YOUNG = 210e9;
RHO = 7860;
alpha = 0; 
beta = 0.1/freq;
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

% 3) Initialize TopOpt parameters
n = msh.nel_dir;
x = vol_frac*ones(prod(n),1);
density = x;
xx = linspace(0,L,n(1));
yy = linspace(0,hh,n(2));

% Density filtering parameters
[dy, dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(method_data.nsub),h,'same');

% Pre-calculate element volume, stiffness, mass and dofs
Ke = zeros(msh.nel,sp.nsh_max,sp.nsh_max);
Me = Ke;
lm = zeros(msh.nel,sp.nsh_max);
for i=1:msh.nel
    [k, m, dofs] = elementaryMatrices(i, sp, sp, msh, ...
        problem_data.lambda_lame, problem_data.mu_lame, problem_data.rho);
    Ke(i,:,:) = k;
    Me(i,:,:) = m;
    lm(i,:) = dofs;
end

tmp_msh = msh_precompute(msh);
Ve = (tmp_msh.element_size.^2)';
V_ = sum(Ve)*vol_frac;
clear tmp_msh
%% Boundary Conditions
% Force
F = buildForce(sp,msh,problem_data);
F = 1000*F/sum(F);

% Dirichlet (clamp)
drch_dofs = sp.boundary(1).dofs; % Index of the dofs of side 1
free_dofs = setdiff(1:sp.ndof, drch_dofs);
u_drch = zeros(numel(drch_dofs),1);

%% TopOpt Initial Conditions
loop = 0;
change = 1;
% Initializing MMA Optimizer
m = mm;
nn = numel(density); % Number of variables to be optimized
xmin = zeros(nn,1); % lower bounds for the variables
xmax = ones(nn,1); % upper bounds for the variables
xold1 = x(:); % xval one iteration ago if iter > 1
xold2 = x(:); % xval two iterations ago if iter > 2
low = ones(nn,1); % Vector with the lower asymptotes from the previous iter
upp = ones(nn,1); % Vector with the upper asymptotes from the previous iter
a0 = 1; % Constants a_0 in the term a_0*z
aa = zeros(m,1); % column vector with the constants a_i 
c_MMA = 1000*ones(m,1); %column vector with the constants c_i*y_i
d = zeros(m,1); %column vector with the constants 0.5*d_i*(y_i)^2


W0 = 0;
csi = 0.1; % Initial csi

while change > change_max && loop < max_iter
    loop = loop+1;

%% Density Filter
    
% Calculate linear filtered density
    xFil = reshape(density,n);
    xFil = conv2(xFil,h,'same')./Hs;
    xFil = xFil(:);
    reset = xFil;

% Calculate volume-preserving physical density using binary search
if filters == 2
    eta1 = 0; eta2 = 1;
    while abs(eta2-eta1) > 1e-4
        xFil = reset;
        eta = 0.5*(eta1+eta2);
        
        xFil = (tanh(csi*eta) +tanh(csi*(xFil-eta)));
        xFil = xFil/(tanh(csi*eta) +tanh(csi*(1-eta)));
        
        if sum(density.*Ve -xFil.*Ve ) > 0
            eta2 = eta;
        elseif sum(density.*Ve -xFil.*Ve ) < 0
            eta1 = eta;
        end
    end
end
xPhys = xFil(:);
% xPhys = reset;
%% Finite-Element Routines

% Pre-alocating vectors/solutions
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
M = K;

% Build Mass, Stiffness and Damping
x = xPhys;
for e=1:msh.nel
    k_e = (Emin +(x(e)^3)*(E*YOUNG - Emin))*squeeze(Ke(e,:,:));
    if x(e) > 0.1
        m_e = (rhomin +(x(e))*(rho0*RHO-rhomin))*squeeze(Me(e,:,:));
    elseif x(e) <= 0.1
        m_e = (rhomin +(x(e)^9)*(rho0*RHO-rhomin))*squeeze(Me(e,:,:));
    end
    idx = lm(e,:)';
    K(idx,idx) = K(idx,idx) +k_e;
    M(idx,idx) = M(idx,idx) +m_e;
end
K = sparse(K); M = sparse(M);
C = alpha*M +beta*K;
Kd = K +1*j*omega*C -omega*omega*M;

% Dirichlet Boundary Conditions
u(drch_dofs) = u_drch;
us = u;
FF = F;

% Applying lifting
F(free_dofs) = F(free_dofs) -Kd(free_dofs,drch_dofs)*u_drch;
FF(free_dofs) = FF(free_dofs) -K(free_dofs,drch_dofs)*u_drch;

% Solving
u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs); % Dynamic
us(free_dofs) = K(free_dofs,free_dofs)\FF(free_dofs); % Static

%% Objective Function and Restrictions

Ek = real(0.5*omega*omega*(u')*M*u); % Kinetic Energy
Ep = real(0.5*(u')*K*u); % Potential Energy

% Objective Function
W = 2*alpha*Ek +2*omega*omega*Ep;

if loop == 1
    W0_db = c0 +10*log10(W);
end

W_db = c0 +10*log10(W);
W_scaled = 100*W_db/W0_db;

f0val = W_scaled;
% f0val = W_db;

% Restrictions

R = Ep/Ek;
R_db = 100 +10*log(R);

Cs = 0.5*us'*K*us;
if loop == 1
    Cs0 = Cs;
end

V_ = sum(x.*Ve);

g1 = 100*(V_/sum(Ve) -vol_frac);
g2 = 100*(Cs/Cs0) -100*kappa;
g3 = -R_db +100;
% g1 = V_/sum(Ve) - vol_frac;
% g2 = Cs - kappa*Cs0;

fval = [g1; g2; g3];
fval = fval(restrictions);

%% Sensitivities

% Adjoint problems
LL = -(1/(2*Ek))*(K -omega*omega*R*M)*conj(u); % Adjoint analysis for R restriction
lambda_ = Kd\LL;

% Element sensitivities
dW = zeros(prod(n),1);
dR = dW;
dCs = dW;
for e=1:msh.nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    m_e = squeeze(Me(e,:,:));
    
    dk = (3*x(e)^2)*(E*YOUNG -Emin)*k_e;
    if x(e) > 0.1
        dm = (rho0*RHO -rhomin)*m_e;
    elseif x(e) <= 0.1
        dm = 9*(x(e)^8)*(rho0*RHO -rhomin)*m_e;
    end
    dc = alpha*dm +beta*dk;
    
    dkd = (dk +1j*dc -omega*omega*dm);
    
    dW(e) = -0.5*omega*real(1j*(u(dofs)')*dkd*u(dofs));
    
    dR(e) = real((0.25/Ek)*(u(dofs)'*(dk -omega*omega*R*dm)*u(dofs)) + ...
        (lambda_(dofs)')*dkd*u(dofs));
    
    dCs(e) = -real((us(dofs)'*dk*us(dofs)));
end

dV = Ve;


% Heaviside Filter Chain-Rule
if filters == 2
    d_rho = csi*(sech(csi*(x-eta).^2));
    d_rho = d_rho/(tanh(csi*eta) +tanh(csi*(1-eta)));

    dW = dW.*d_rho;
    dR = dR.*d_rho;
    dCs = dCs.*d_rho;
    dV = dV.*d_rho;
end

% Log transform Chain-Rule
dW = (10/(log(10)*W))*dW;
dR = -(10/log(10)*R)*dR;

% Linear filter Chain-Rule
dW = reshape(dW,n);
dW = conv2(dW./Hs,h,'same');
dW = dW(:);

dR = reshape(dR,n);
dR = conv2(dR./Hs,h,'same');
dR = dR(:);

dCs = reshape(dCs,n);
dCs = conv2(dCs./Hs,h,'same');
dCs = dCs(:);

dV = reshape(dV,n);
dV = conv2(dV./Hs,h,'same');
dV = dV(:);

% Constants Multiplication
dW = 100*dW/W0_db;
dV = (100/sum(Ve))*dV;
dCs = (50/Cs0)*dCs;

%% Method of Moving Asymptotes
xval = x(:);
df0dx = dW;
df0dx2 = zeros(size(df0dx));
dfdx = [dV'; dCs'; dR'];
dfdx = dfdx(restrictions,:);
dfdx2 = zeros(m,nn);
[xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m, nn, loop, xval, xmin, xmax, ...
    xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, ...
    a0, aa, c_MMA, d);
density = xmma(:);
xold2 = xold1;
xold1 = xval;
change = max(abs(xmma(:)-xold1(:)));
volum = mean(xmma(:));
x_plot = conv2(reshape(xmma,n),h,'same')./Hs;
fprintf(' Iteration.:%5i | Objective.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
    loop, f0val, volum,change);
colormap(gray); imagesc(xx,yy,1-rot90(x_plot)); caxis([0 1]); axis equal; axis off; drawnow;
if (mod(loop,update) == 1 || change < 0.01) && csi < 128
    csi = 1.2*csi;
    change = 0.5;
end
end
