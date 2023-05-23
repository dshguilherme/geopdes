L = 1;
hh = 0.5;
YOUNG = 210e9;
RHO = 7860;

problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; Emin = 1e-3; rho0=1; alpha =1e-3; beta=1e-8; rhomin=1e-3;
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
m = 2; % Number of general constraints
nn = numel(density); % Number of variables to be optimized
xmin = zeros(nn,1); % lower bounds for the variables
xmax = ones(nn,1); % upper bounds for the variables
xold1 = x(:); % xval one iteration ago if iter > 1
xold2 = x(:); % xval two iterations ago if iter > 2
low = ones(nn,1); % Vector with the lower asymptotes from the previous iter
upp = ones(nn,1); % Vector with the upper asymptotes from the previous iter
a0 = 1; % Constants a_0 in the term a_0*z
aa = zeros(m,1); % column vector with the constants a_i 
c_MMA = 10000*ones(m,1); %column vector with the constants c_i*y_i
d = zeros(m,1); %column vector with the constants 0.5*d_i*(y_i)^2

while change > change_max && loop < max_iter
loop = loop+1;

% Pre-alocating vectors/solutions
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
M = K;

% Build Mass, Stiffness and Damping
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

% Solving
Kd = K +1i*omega*C -(omega^2)*M;
u(drch_dofs) = u_drch;
F(free_dofs) = F(free_dofs) -Kd(free_dofs,drch_dofs)*u_drch;
u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);

Ek = 0.5*omega*omega*(u')*M*u;
Ep = 0.5*(u')*K*u;
Pa = real(2*alpha*Ek +2*omega*omega*Ep);
R = Ep/Ek;

% Adjoint problem for R
lambda_ = Kd\((-1/(2*Ek))*(K-omega*omega*R*M)*conj(u));

dc = zeros(length(Ve),1);
d2c = dc;
dR = dc;
d2R = dc;

for e=1:msh.nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    m_e = squeeze(Me(e,:,:));
    dkc = -3*(x(e)^2)*(E*YOUNG-Emin)* (u(dofs)')*k_e*u(dofs);
    dlc = -3*(x(e)^2)*(E*YOUNG-Emin)* (lambda_(dofs)')*k_e*u(dofs);
    d2kc = -6*x(e)*(E*YOUNG-Emin)*(u(dofs)')*k_e*u(dofs);
    d2lc = -6*x(e)*(E*YOUNG-Emin)*(lambda_(dofs)')*k_e*u(dofs);
    if x(e) > 0.1
        dmc = -(rho0*RHO-rhomin)*(u(dofs)')*m_e*u(dofs);
        dlmc = -(rho0*RHO-rhomin)*(lambda_(dofs)')*m_e*u(dofs);
        d2mc = 0;
        d2lmc = 0;
    elseif x(e) <= 0.1
        dmc = -9*(x(e)^8)*(rho0*RHO-rhomin)*(u(dofs)')*m_e*u(dofs);
        dlmc = -9*(x(e)^8)*(rho0*RHO-rhomin)*(lambda_(dofs)')*m_e*u(dofs);
        d2mc = -72*(x(e)^7)*(rho0*RHO-rhomin)*(u(dofs)')*m_e*u(dofs);
        d2lmc = -72*(x(e)^7)*(rho0*RHO-rhomin)*(lambda_(dofs)')*m_e*u(dofs);
    end
    dcc = alpha*dmc +beta*dkc;
    dlcc = alpha*dlc +beta*dlmc;
    d2cc = alpha*d2mc +beta*d2kc;
    d2lcc = alpha*d2lc +beta*d2lmc;
    
    dc(e) = -0.5*omega*real(1i*(dkc +1i*omega*dcc -(omega^2)*dmc));
    d2c(e) = -0.5*omega*real(1i*(d2kc +1i*omega*d2cc -(omega^2)*d2mc));

    dR(e) = (0.25/Ek)*(dkc -omega*omega*R*dmc) +real(dlc+ 1i*omega*dlcc -omega*omega*dlmc);
    d2R(e) = (0.25/Ek)*(d2kc -omega*omega*R*d2mc) +real(d2lc+j*omega*d2lcc+omega*omega*d2lmc);
end
v_xe = x.*Ve;
Re = dc./Ve;

% Density filter
xPhys = reshape(density,n);
dc_xy = reshape(dc,n);
d2c_xy = reshape(d2c,n);
dc = conv2(dc_xy./Hs,h,'same');
d2c = conv2(d2c_xy./Hs,h,'same');
dv = reshape(Ve,n);
dv = conv2(dv./Hs,h,'same');

% Update x value through Method of Moving Asymptotes
    xval = xPhys(:);
%     m = 1;
    %nn = length(xval), iter = loop, xmin, xmax, xold1, xold2
    f0val = 100 +10*log10(Pa);
    df0dx = (10/(Pa*log(10)))*dc(:);
    df0dx2 = (10/(Pa*Pa*log(10)))*(f0val*d2c(:) -dc(:).^2);
    fval = real([100*(sum(xPhys(:))/(vol_frac*nn) -1)]);% (100 +10*log10(R))-100]);
    dfdx = real([100*dv(:)'/(vol_frac*nn)]);% (10/(R*log(10)))*dR(:)']);
    dfdx2 = real([zeros(size(dv(:)'))]);% (10/(R*R*log(10)))*(fval(2)*d2R(:)' -dfdx(2,:).^2)]);
    [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m, nn, loop, xval, xmin, xmax, ...
        xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, ...
        a0, aa, c_MMA, d);
    xPhys = conv2(reshape(xmma,n),h,'same')./Hs;
    density = xmma(:);
    xold2 = xold1;
    xold1 = xval;
    change = max(abs(xmma(:)-x(:)));
    volum = (mean(xmma(:)));

x = xPhys(:);
fprintf(' Iteration.:%5i | Compliance.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
    loop, f0val, volum,change);
colormap(gray); imagesc(xx,yy,1-rot90(xPhys)); caxis([0 1]); axis equal; axis off; drawnow;

end