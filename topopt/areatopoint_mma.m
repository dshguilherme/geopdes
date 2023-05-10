% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
L = 1;
hh = .5;
d = L/5;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions for each side of the domain
problem_data.symm_sides = [4];
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.press_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
E = 1; Emin = 1e-3;

% Source and boundary terms
problem_data.f = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [1 1];       % Degree of the splines
method_data.regularity = [0 0];       % Regularity of the splines
method_data.nsub       = [30 10];       % Number of subdivisions
method_data.nquad      = [2 2];       % Points for the Gaussian quadrature rule
rmin = 3;

% Build spaces
[geometry, msh, ~] = buildSpaces(problem_data, method_data);
space = sp_nurbs(geometry.nurbs, msh);
sp = space;

% 3) Initialize TopOpt parameters
vol_frac = 0.5;
n = msh.nel_dir;
x = vol_frac*ones(prod(n),1);
density = x;
xx = linspace(0,L,n(1));
yy = linspace(0, hh, n(2));

% Density filtering parameters

[dy, dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(method_data.nsub),h,'same');

% Pre-calculate element volume, stiffness and dofs
Ke = zeros(msh.nel,sp.nsh_max,sp.nsh_max);
lm = zeros(msh.nel,sp.nsh_max);
for i=1:msh.nel
    [k, dofs] = conductionStiffness(i,sp,sp,msh,problem_data.c_diff);
    Ke(i,:,:) = k;
    lm(i,:) = dofs;
end
tmp_msh = msh_precompute(msh);
Ve = (tmp_msh.element_size.^2)';

%% Boundary Conditions
F = buildForce(sp,msh,problem_data);
F = F/sum(F);

% Dirichlet

b1_dofs = sp.boundary(1).dofs;
dy_dof = hh/msh.nel_dir(2);
ny1 = ceil(d/dy_dof); % # of dofs under Dirichlet BCs;

drch_dofs = b1_dofs(end+1-ny1:end);

free_dofs = setdiff(1:sp.ndof, [drch_dofs]);
u_drchlt = zeros(numel(drch_dofs),1);

%% TopOpt Initial Conditions
loop = 0;
change = 1;
% Initializing MMA Optimizer
m = 1; % Number of general constraints
nn = numel(density); % Number of variables to be optimized
xmin = zeros(nn,1); % lower bounds for the variables
xmax = ones(nn,1); % upper bounds for the variables
xold1 = x(:); % xval one iteration ago if iter > 1
xold2 = x(:); % xval two iterations ago if iter > 2
low = ones(nn,1); % Vector with the lower asymptotes from the previous iter
upp = ones(nn,1); % Vector with the upper asymptotes from the previous iter
a0 = 1; % Constants a_0 in the term a_0*z
a = zeros(m,1); % column vector with the constants a_i 
c_MMA = 10000*ones(m,1); %column vector with the constants c_i*y_i
d = zeros(m,1); %column vector with the constants 0.5*d_i*(y_i)^2


while change > 0.01 && loop < 200
    loop = loop+1;
    % Pre-alocate arrays
    u = zeros(sp.ndof,1);
    K = zeros(sp.ndof);
    % Build stiffness
    for e=1:msh.nel
        k_e = (Emin +(x(e)^3)*(E-Emin))*squeeze(Ke(e,:,:));
        idx = lm(e,:)';
        K(idx,idx) = K(idx,idx)+k_e;
    end
    % Solving
    u(drch_dofs) = u_drchlt;
    F(free_dofs) = F(free_dofs) -K(free_dofs, drch_dofs)*u_drchlt;
    u(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);
    
    compliance = F'*u;
    
    dc = zeros(msh.nel,1);
    d2c = zeros(msh.nel,1);
    for e=1:msh.nel
        dofs = lm(e,:)';
        k_e = squeeze(Ke(e,:,:));
        dc(e) = -3*(x(e)^2)*(E-Emin)*(u(dofs)'*k_e*u(dofs));
        d2c(e) = -6*(x(e))*(E-Emin)*(u(dofs)'*k_e*u(dofs));
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
m = 1;
%nn = length(xval), iter = loop, xmin, xmax, xold1, xold2
f0val = compliance;
df0dx = dc(:);
df0dx2 = d2c(:);
fval = sum(xPhys(:))/(vol_frac*nn) -1;
dfdx = dv(:)'/(vol_frac*nn);
dfdx2 = zeros(m,nn);
[xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m, nn, loop, xval, xmin, xmax, ...
    xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, ...
    a0, a, c_MMA, d);
xPhys = conv2(reshape(xmma,n),h,'same')./Hs;
density = xmma(:);
xold2 = xold1;
xold1 = xval;
change = max(abs(xmma(:)-x(:)));
x = xPhys(:);

fprintf(' Iteration.:%5i | Compliance.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
    loop, compliance, mean(xmma(:)),change);
colormap(gray); imagesc(xx,yy,1-rot90(xPhys)); caxis([0 1]); axis equal; axis off; drawnow;

end