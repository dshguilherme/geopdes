function mbb_beam(degree, nsub, vol_frac, rmin, method, change_max, max_iter)
L = 60;
hh = 20;
d = 0.05*L;

problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [4];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [1];

% Physical parameters
E  =  1; Emin = 1e-3;
problem_data.E = E;
problem_data.Emin = Emin;
nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @forceMBBDistributed;
% problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
% problem_data.g = @forceMBB;
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));


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

% Pre-calculate element volume, stiffness and dofs
Ke = zeros(msh.nel,sp.nsh_max,sp.nsh_max);
lm = zeros(msh.nel,sp.nsh_max);
for i=1:msh.nel
    [k, dofs] = elementaryStiffness(i, sp, sp, msh, problem_data.lambda_lame, problem_data.mu_lame);
    Ke(i,:,:) = k;
    lm(i,:) = dofs;
end

tmp_msh = msh_precompute(msh);
Ve = (tmp_msh.element_size.^2)';
clear tmp_msh
%% Boundary Conditions
% Force
F = buildForce(sp,msh,problem_data);
F = F/sum(F);

% Dirichlet and Symmetry
symm_dofs = findSymm(sp,msh,problem_data);

y_dofs = sp.comp_dofs{2}; % Index of y dofs
b2y_dofs = sp.boundary(2).dofs; % dofs of side 2 of the boundary
b3y_dofs = sp.boundary(3).dofs; % dofs of side 3 of the boundary

dx_dof = L/msh.nel_dir(1);
dy_dof = hh/msh.nel_dir(2);

n_drch_x = ceil(d/dx_dof);
n_drch_y = ceil(d/dy_dof);



drchlt_dofs2 = intersect(y_dofs, b2y_dofs);
drchlt_dofs3 = intersect(y_dofs, b3y_dofs);

drchlt_dofs = unique([drchlt_dofs2(1:n_drch_y); drchlt_dofs3(end+1-n_drch_x:end)]);

free_dofs = setdiff(1:sp.ndof, [drchlt_dofs; symm_dofs]);

u_drchlt = zeros(numel(drchlt_dofs),1); % For this case, homogeneous BCs.

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


while change > change_max && loop < max_iter
loop = loop+1;
% Pre-alocating vectors/solutions
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
% Build stiffness
for e=1:msh.nel
    k_e = (Emin +(x(e)^3)*(E - Emin))*squeeze(Ke(e,:,:));
    idx = lm(e,:)';
    K(idx,idx) = K(idx,idx) + k_e;
end
% Solving

u(drchlt_dofs) = u_drchlt;
F(free_dofs) = F(free_dofs) -K(free_dofs, drchlt_dofs)*u_drchlt;
u(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);


compliance = F'*u;
dc = zeros(msh.nel,1);
d2c = zeros(msh.nel,1);
for e=1:msh.nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    dc(e) = -3*(x(e)^2)*(E-Emin)* (u(dofs)')*k_e*u(dofs);
    d2c(e) = -6*x(e)*(E-Emin)*(u(dofs)')*k_e*u(dofs);
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
if strcmp(method,'MMA')
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
    volum = (mean(xmma(:)));

elseif strcmp(method,'OC')
    ell1 = 0; ell2 = 1e9; move = 0.2;
    while(ell2-ell1)/(ell1+ell2) > 1e-3
        mid = 0.5*(ell2+ell1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc(:)./dv(:)/mid)))));
        xPhys = conv2(reshape(xnew,n),h,'same')./Hs;
        if sum(xPhys(:)) > vol_frac*length(x)
            ell1 = mid;
        else
            ell2 = mid;
        end
    end
    change = max(abs(xnew(:)-x(:)));
    volum = mean(xnew(:));
end
x = xPhys(:);

fprintf(' Iteration.:%5i | Compliance.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
    loop, compliance, volum,change);
colormap(gray); imagesc(xx,yy,1-rot90(xPhys)); caxis([0 1]); axis equal; axis off; drawnow;

end
end