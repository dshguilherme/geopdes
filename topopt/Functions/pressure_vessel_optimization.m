A = 470e-3;
B = 1775e-3;
YOUNG = 180e6; % Steel young's modulus 180MPa
RHO = 7850; % Steel density 7850 kg/m3
problem_data = pressure_vessel_problem(A,B);
iter_max = 500;

method_data.degree = [2 2]; % Degree of the BSplines
method_data.regularity = [1 1];
method_data.nsub = [1 1];
method_data.nquad = [3 3]; % number of quadrature points

% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data, method_data);

%% Pre-calc element volume, stiffness, mass and force

msh1 = msh_precompute(msh);
sp1 = sp_precompute_param(sp,msh1,'value', true,'gradient', true, 'hessian', true);
for idim = 1:msh1.rdim
    x{idim} = reshape(msh1.geo_map(idim,:,:),msh1.nqn,msh1.nel);
end

% Stiffness Matrix
[Bke, Ske] = op_KL_shells_elements(sp1, sp1, msh1, problem_data.E_coeff(x{:}), problem_data.nu_coeff(x{:}));

% Mass
Me = op_u_v_elements(sp1, sp1, msh1, RHO);
lm = sp1.connectivity';

% Restriction
Ve = (msh1.element_size.^2); % Element Area

% Force vector
F = op_pn_v(sp1,msh1,problem_data.p(x{:})); % Pressure

%% get MS Relations
% degen_dofs = [sp.boundary(3).dofs sp.boundary(4).dofs];
% master = degen_dofs(1:sp.boundary(3).ndof_dir:end);
% slaves = setdiff(degen_dofs,master);

side_slaves = sp.boundary(2).dofs;
side_master = sp.boundary(1).dofs;

% [~, ~, ib] = intersect(slaves,side_slaves);
% side_slaves(ib) = [];
% side_master(ib) = [];

% ms = zeros(length(slaves),2);

% ms(:,2) = slaves;
% ms(:,1) = sort(repmat(master',length(slaves)/6,1));
% ms = [ms; side_master(:), side_slaves(:)];
% ms = sort(ms);
% ms = zeros(size(ms));
ms = [side_master(:), side_slaves(:)];

%% Generate Bending Strip
coefs = geometry.nurbs.coefs(:,[1 end-1 end-2],:);
knots = {[0 0 0 1 1 1], geometry.nurbs.knots{2}};
strip = nrbmak(coefs,knots);
strip_dofs = sp.boundary(1).dofs;
nn = strip.number(2);
strip_dofs = [strip_dofs', strip_dofs'+nn-1, strip_dofs'+nn-2];

[bK, b_dofs] = bending_strip_stiffness(strip);

%% Set No displacement conditions
cpoints = reshape(geometry.nurbs.coefs,4,numel(geometry.nurbs.coefs)/4);
% Find all controlpoints with y below -A/2
idx_y = find(cpoints(2,:) <= -A/2);

% Find all controlpoints B/4 away from the origin, with a 10% margin
idx_x1 = find(cpoints(1,:) <= B/4+B/40);
idx_x2 = find(cpoints(1,:) >= B/4-B/40);
idx_x = intersect(idx_x1,idx_x2);

idx_x3 = find(cpoints(1,:) >= -B/4-B/40);
idx_x4 = find(cpoints(1,:) <= -B/4+B/40);
neg_idx_x = intersect(idx_x3,idx_x4);

idx_x = union(idx_x,neg_idx_x);

% Intersection of x and y points
idx = intersect(idx_x,idx_y);

dr_dofs = [idx, idx+sp.cumsum_ndof(2), idx+sp.cumsum_ndof(3)];
free_dofs = setdiff(1:sp.ndof,dr_dofs);
free_dofs = setdiff(free_dofs,ms(:,2));

%% Thickness
tmax = 6e-3; %maximum thickness allowed
thickness = 3.5e-3; % mean thickness
tmin = 1e-3; %min thickness allowed
factor = 100*(thickness-tmin)/(tmax-tmin);
dfactor = (tmax-tmin)/100;

%% Initialize MMA constants
m = 2; % Number of Restrictions
n = msh.nel; % Number of variables
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = factor*eeen;
xold1   = xval;
xold2   = xval;
xmin    = 0*eeen;
xmax    = 100*eeen;
low     = xmin;
upp     = xmax;
c       = 10000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
raa0    = 0.01;
raa     = 0.01*eeem;
raa0eps = 0.000001;
raaeps  = 0.000001*eeem;
outeriter = 0;
maxoutit  = iter_max;
kkttol  = 0;

%% Initial Solution
t = (tmax-tmin)*xval/100 +tmin;
[Ks, M] = shellMatricesFromElements(Bke, Ske, Me, lm, t, YOUNG, RHO);
Ks(strip_dofs, strip_dofs) = Ks(strip_dofs,strip_dofs) + bK(b_dofs,b_dofs);
C = alpha_*M +beta_*Ks;
Kd = Ks +1j*omega*C -omega*omega*M;


u = zeros(length(Ks),1);
u(free_dofs) = V(:,1);
u(ms(:,2)) = u(ms(:,1));
u(dr_dofs) = 0;

deformed_vessel = geo_deform(u,sp,geometry);
nrbplot(deformed_vessel.nurbs,[100 100])