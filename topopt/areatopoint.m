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
method_data.nsub       = [120 40];       % Number of subdivisions
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
xold1 = x;
xold2 = x;
low = ones(size(x));
upp = ones(size(x));
while change > 0.01
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
    for e=1:msh.nel
        dofs = lm(e,:)';
        k_e = squeeze(Ke(e,:,:));
        dc(e) = -3*(x(e)^2)*(E-Emin)*(u(dofs)'*k_e*u(dofs));
    end
    v_xe = x.*Ve;
    Re = dc./Ve;
    xPhys = reshape(density,n);
    % Density filter
    dc_xy = reshape(dc,n);
    dc = conv2(dc_xy./Hs,h,'same');
    dv = reshape(Ve,n);
    dv = conv2(dv./Hs,h,'same');
    % Calc Lagrange Multiplier and update x
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
density = xnew(:);
x = xPhys(:);
fprintf(' Iteration.:%5i | Compliance.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
    loop, compliance, mean(xnew(:)),change);
colormap(gray); imagesc(xx,yy,1-rot90(xPhys)); caxis([0 1]); axis equal; axis off; drawnow;

end