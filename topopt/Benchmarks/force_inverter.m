function force_inverter(degree, nsub, vol_frac, rmin, method, change_max, max_iter)
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
L = 1;
hh = L/2;
d = L/100;
k_in = 1;
k_out = k_in/100;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
problem_data.symm_sides   = [4];

% Physical parameters
E  =  1; nu = .3; Emin = 1e-9;
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule
rmin = 2;
% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data, method_data);

% 3) Initialize TopOpt parameters

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
Ke = zeros(msh.nel, sp.nsh_max, sp.nsh_max);
lm = zeros(msh.nel, sp.nsh_max);
for i=1:msh.nel
    [k, dofs] = elementaryStiffness(i, sp, sp, msh, problem_data.lambda_lame, problem_data.mu_lame);
    Ke(i,:,:) = k;
    lm(i,:) = dofs;
end

tmp_msh = msh_precompute(msh);
Ve = (tmp_msh.element_size.^2)';
clear tmp_msh

%% Boundary Conditions
x_dofs = sp.comp_dofs{1}; % Index of x dofs
y_dofs = sp.comp_dofs{2}; % Index of y dofs

% Symmetry
symm_dofs = findSymm(sp, msh, problem_data);

% Boundary 1
b1_dofs = sp.boundary(1).dofs;
b1x = intersect(x_dofs,b1_dofs); % x dofs of boundary 1 
b1y = intersect(y_dofs,b1_dofs); % y dofs of boundary 1

dy_dof = 0.5*L/msh.nel_dir(2);
ny1 = ceil(d/dy_dof); % # of dofs of the clamp

drch_b1 = [b1x(1:ny1), b1y(1:ny1)]';

% Boundary 4
b4_dofs = sp.boundary(4).dofs;

b4x = intersect(x_dofs, b4_dofs);

k1_dof = intersect(b1x, b4x); % dof of spring 1 and f_in

% Boundary 2
b2_dofs = sp.boundary(2).dofs;
b2x = intersect(x_dofs, b2_dofs);

k2_dof = intersect(b2x, b4x); % dof of the spring 2 and u_out


% Spring constants
K = sparse(zeros(sp.ndof));

% Force
F = zeros(sp.ndof,1);
LL = F;

LL(k2_dof) = 1;
F(k1_dof) = 1;

% Dirichlet and Symmetry BCs
drchlt_dofs = drch_b1;
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
        k_e = (Emin +(x(e)^3)*(E -Emin))*squeeze(Ke(e,:,:));
        idx = lm(e,:)';
        K(idx,idx) = K(idx,idx) +k_e;
    end
    
    K(k1_dof,k1_dof) = K(k1_dof,k1_dof) +k_in; % Adding spring to k1
    K(k2_dof,k2_dof) = K(k2_dof,k2_dof) +k_out; % Adding spring to k2
    % Solving
    u(drchlt_dofs) = u_drchlt;
    F(free_dofs) = F(free_dofs) -K(free_dofs, drchlt_dofs)*u_drchlt;
    u(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);
    
    compliance = LL'*u;
    adjoint_lambda = zeros(sp.ndof,1);
    adjoint_lambda(free_dofs) = K(free_dofs,free_dofs)\LL(free_dofs);
%     adjoint_lambda = K\-LL;
    dc = zeros(msh.nel,1);
      d2c = zeros(msh.nel,1);
    for e=1:msh.nel
        dofs = lm(e,:)';
        k_e = squeeze(Ke(e,:,:));
        lambda_ = adjoint_lambda(dofs);
        dc(e) = -3*(x(e)^2)*(E-Emin)*(lambda_'*k_e*u(dofs));
        d2c(e) = -6*(x(e))*(E-Emin)*(lambda_'*k_e*u(dofs));
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
    
    if strcmp(method, 'MMA')
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
        volum = mean(xmma(:));
    elseif strcmp(method, 'OC')
     % Calc Lagrange Multiplier and update x
        ell1 = 0; ell2 = 1e9; move = 0.1;
        while(ell2-ell1)/(ell1+ell2) > 1e-4 && ell2 > 1e-40
            mid = 0.5*(ell2+ell1);
            xnew = max(0,max(x-move,min(1,min(x+move,x.*max(1e-10,-dc(:)./dv(:)/mid).^.3))));
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