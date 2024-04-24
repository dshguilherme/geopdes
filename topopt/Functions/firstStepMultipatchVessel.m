function initial_output = firstStepMultipatchVessel(paramaters, problem_data)
data_names = fieldnames(parameters);
for iopt = 1:numel(data_names)
    eval([data_names{iopt} '=parameters.(data_names{iopt});']);
end
clear method_data
method_data.degree = [2 degree];
method_data.regularity = [1 degree-1];
method_data.nsub = [1 nsub];
method_data.nquad = [2 degree+1];
% Build Spaces
[geo1, msh1, sp1] = buildSpaces(problem_data(1),method_data);
[geo2, msh2, sp2] = buildSpaces(problem_data(2),method_data);
[geo3, msh3, sp3] = buildSpaces(problem_data(3),method_data);

%% Pre-calc element volume, stiffness, mass and force

pmsh1 = msh_precompute(msh1);
psp1 = sp_precompute_param(sp1,pmsh1,'value', true,'gradient', true, 'hessian', true);
pmsh2 = msh_precompute(msh2);
psp2 = sp_precompute_param(sp2,pmsh2,'value', true,'gradient', true, 'hessian', true);
pmsh3 = msh_precompute(msh3);
psp3 = sp_precompute_param(sp3,pmsh3,'value', true,'gradient', true, 'hessian', true);
for idim = 1:pmsh1.rdim
    x1{idim} = reshape(pmsh1.geo_map(idim,:,:),pmsh1.nqn,pmsh1.nel);
    x2{idim} = reshape(pmsh2.geo_map(idim,:,:),pmsh2.nqn,pmsh2.nel);
    x3{idim} = reshape(pmsh3.geo_map(idim,:,:),pmsh3.nqn,pmsh3.nel);
end

% Stiffness Matrices
[Bke1, Ske1] = op_KL_shells_elements(psp1, psp1, pmsh1, problem_data(1).E_coeff(x1{:}), problem_data(1).nu_coeff(x1{:}));
[Bke2, Ske2] = op_KL_shells_elements(psp2, psp2, pmsh2, problem_data(2).E_coeff(x2{:}), problem_data(2).nu_coeff(x2{:}));
[Bke3, Ske3] = op_KL_shells_elements(psp3, psp3, pmsh3, problem_data(3).E_coeff(x3{:}), problem_data(3).nu_coeff(x3{:}));

% Mass
Me1 = op_u_v_elements(psp1, psp1, pmsh1, RHO);
Me2 = op_u_v_elements(psp2, psp2, pmsh2, RHO);
Me3 = op_u_v_elements(psp3, psp3, pmsh3, RHO);
lm1 = psp1.connectivity';
lm2 = psp2.connectivity';
lm3 = psp3.connectivity';

% Restrictions
Ve1 = (pmsh1.element_size.^2); % Element Area
Ve2 = (pmsh2.element_size.^2); % Element Area
Ve3 = (pmsh3.element_size.^2); % Element Area

% Force vector
F1 = pressure*op_pn_v(psp1,pmsh1,problem_data(1).p(x1{:})); % Pressure
F2 = pressure*op_pn_v(psp2,pmsh2,problem_data(2).p(x2{:})); % Pressure
F3 = pressure*op_pn_v(psp3,pmsh3,problem_data(3).p(x3{:})); % Pressure

%% Generate Bending Strips
PP1 = [0 0 0];
PP2 = [0 0 0];
line = nrbline(PP1,PP2);
[~, sz1, sz2] = size(geo1.nurbs.coefs);
interval = linspace(0,1,sz1);
strip = nrbkntins(line,interval(2:end-1));
strip.number = [strip.number 3];
strip.knots = {strip.knots, [0 0 0 1 1 1]};
strip.order = [2 3];

strip1 = strip;
coefs = zeros(4,length(interval),3);
coefs(:,:,1) = geo1.nurbs.coefs(:,:,end-1);
coefs(:,:,2) = geo1.nurbs.coefs(:,:,end);
coefs(:,:,3) = geo3.nurbs.coefs(:,:,2);
strip1.coefs = coefs;

strip2 = strip;
coefs = zeros(4,length(interval),3);
coefs(:,:,1) = geo2.nurbs.coefs(:,:,end-1);
coefs(:,:,2) = geo2.nurbs.coefs(:,:,end);
coefs(:,:,3) = geo3.nurbs.coefs(:,:,end-1);
strip2.coefs = coefs;

[bK1, ssp1] = bending_strip_stiffness(strip1);
[bK2, ssp2] = bending_strip_stiffness(strip2);

%% DOF mapping from strip to patch, patch to global dofs
% Boundary DOFs of the bending strip -> patch1 / patch2
% setdiff(strip_dofs, boundary_dofs) -> common boundary
strip1_dofs = cell2mat(ssp1.comp_dofs(1:3));
sb1_dofs = union(ssp1.boundary(3).dofs, ssp1.boundary(4).dofs);
common_bdr1 = setdiff(strip1_dofs, sb1_dofs);
strip1_dofs = [ssp1.boundary(3).dofs(:), common_bdr1(:), ssp1.boundary(4).dofs(:)];
strip2_dofs = cell2mat(ssp2.comp_dofs(1:3));
sb2_dofs = union(ssp2.boundary(3).dofs, ssp2.boundary(4).dofs);
common_bdr2 = setdiff(strip2_dofs, sb2_dofs);
strip2_dofs = [ssp2.boundary(3).dofs(:), common_bdr2(:), ssp2.boundary(4).dofs(:)];

%boundary 4 = calota's circle
%boundary 3 = cylinder circle (positive)
%boundary 4 = cylinder circle (negative)
%calota1 = positive, calota2 = negative
% format = [cap, common, cylinder]
patch_strip1_dofs = [sp1.boundary(4).dofs(:)-geo1.nurbs.number(1),sp1.boundary(4).dofs(:),sp3.boundary(3).dofs(:)+geo3.nurbs.number(1)];
patch_strip2_dofs = [sp2.boundary(4).dofs(:)-geo2.nurbs.number(1),sp2.boundary(4).dofs(:),sp3.boundary(4).dofs(:)-geo3.nurbs.number(1)];

cap1_pad = 0;
cap2_pad = sp1.ndof;
cylinder_pad = sp1.ndof+sp2.ndof;

global_strip1_dofs = patch_strip1_dofs+[cap1_pad, cap1_pad, cylinder_pad];
gs1_dofs = global_strip1_dofs(:);
global_strip2_dofs = patch_strip2_dofs+[cap2_pad, cap2_pad, cylinder_pad];
gs2_dofs = global_strip2_dofs(:);

%% Dirichlet Boundary Conditions
% Degenerate points due to revolving the geometry
degen1_dofs = sp1.boundary(3).dofs;
master_degen1 = sp1.boundary(3).dofs(1)+sp1.cumsum_ndof(1:3);
degen1_dofs = setdiff(degen1_dofs,master_degen1);
degen2_dofs = sp2.boundary(3).dofs+cap2_pad;
master_degen2 = sp2.boundary(3).dofs(1)+sp2.cumsum_ndof(1:3)+cap2_pad;
degen2_dofs = setdiff(degen2_dofs,master_degen2);

ms = zeros(length(degen1_dofs)+length(degen2_dofs),2);
tmp = repmat(master_degen1,length(degen1_dofs)/3,1);
ms(1:length(degen1_dofs),1) = sort(tmp(:));
ms(1:length(degen1_dofs),2) = degen1_dofs(:);

tmp = repmat(master_degen2, length(degen2_dofs)/3,1);
ms(length(degen1_dofs)+1:end,1) = sort(tmp(:));
ms(length(degen1_dofs)+1:end,2) = degen2_dofs;

% Degenerate boundaries due to revolving the geometry
degen1_line_dofs = sp1.boundary(2).dofs;
degen1_line_master = sp1.boundary(1).dofs;
degen2_line_dofs = sp2.boundary(2).dofs+cap2_pad;
degen2_line_master = sp2.boundary(1).dofs+cap2_pad;
degen3_line_dofs = sp3.boundary(2).dofs+cylinder_pad;
degen3_line_master = sp3.boundary(1).dofs+cylinder_pad;

ms = [ms; degen1_line_master(:), degen1_line_dofs(:); ...
    degen2_line_master(:), degen2_line_dofs(:); ...
    degen3_line_master(:), degen3_line_dofs(:)];

% Find duplicate rows in ms and delete them
[~, I, ~] = unique(ms, 'rows', 'first');
duplicates = setdiff(1:size(ms,1), I);
ms(duplicates,:) = [];
ms = sort(ms,1);

% 0 Displacement lines at the supports
x_coords = squeeze(geo3.nurbs.coefs(1,:,:));
y_coords = squeeze(geo3.nurbs.coefs(2,:,:));

% Find all points with y below -A/2:
idx_y = find(y_coords <=-A/2);

% Find all points B/4 away from the origin, with a 10% margin
idx_x1 = find(x_coords <= B/4 +B/40);
idx_x2 = find(x_coords >= B/4 -B/40);
idx_x = intersect(idx_x1, idx_x2);

idx_x3 = find(x_coords >= -B/4-B/40);
idx_x4 = find(x_coords <= -B/4+B/40);
idx_x_neg = intersect(idx_x3, idx_x4);

idx_x = union(idx_x, idx_x_neg);

% Find the points that apply to all conditions
idx = intersect(idx_x,idx_y);

dr_dofs = idx+cylinder_pad +sp3.cumsum_ndof(1:3);
dr_dofs = dr_dofs(:);

all_dofs = 1:sp1.ndof+sp2.ndof+sp3.ndof;

restricted_dofs = unique([degen1_dofs(:); degen2_dofs(:); degen1_line_dofs(:); ...
    degen2_line_dofs(:); degen3_line_dofs(:); dr_dofs]);
free_dofs = setdiff(all_dofs,restricted_dofs);

%% Initialize MMA constants
tmax = max_thickness;
tmin = min_thickness;
factor = 100*(thickness-tmin)/(tmax-tmin);
dfactor = (tmax-tmin)/100;
m = 2; % Number of Restrictions
n = msh1.nel+msh2.nel+msh3.nel; % Number of variables
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
t1 = t(1:msh1.nel);
t2 = t(msh1.nel+1:msh1.nel+msh2.nel);
t3 = t(msh1.nel+msh2.nel+1:end);
% Matrices
[Ks1, M1] = shellMatricesFromElements(Bke1, Ske1, Me1, lm1, t1, YOUNG, RHO);
[Ks2, M2] = shellMatricesFromElements(Bke2, Ske2, Me2, lm2, t2, YOUNG, RHO);
[Ks3, M3] = shellMatricesFromElements(Bke3, Ske3, Me3, lm3, t3, YOUNG, RHO);

Ks = blkdiag(Ks1,Ks2,Ks3);
M = blkdiag(M1,M2,M3);
F = [F1; F2; F3];

% Contribution from bending strips
Ks(gs1_dofs,gs1_dofs) = Ks(gs1_dofs,gs1_dofs) +bK1(strip1_dofs(:),strip1_dofs(:));
Ks(gs2_dofs,gs2_dofs) = Ks(gs2_dofs,gs2_dofs) +bK2(strip2_dofs(:),strip2_dofs(:));

C = alpha_*M +beta_*Ks;
omega = 2*pi*freq;
Kd = Ks +1j*omega*C -omega*omega*M;

% Solution
u_init = zeros(length(Kd),1);
F(free_dofs) = F(free_dofs) -Kd(free_dofs, dr_dofs)*u_init(dr_dofs);
u_init(free_dofs) = Kd(free_dofs, free_dofs) \ F(free_dofs);
u_init(ms(:,2)) = u_init(ms(:,1)); % Master/slave displacements

% Objective Functions
aW_init = real(0.5*omega*omega*(u_init'*C*u_init));
velocity0 = -1j*omega*u_init;
v2_init = real(velocity0'*velocity0);

f1 = @multiPatchVesselObjectivesAndSensitivvies;
f2 = @multiPatchVesselObjectives;

s = whos;
initial_output = cell2struct({s.name}.',{s.name});
eval(structvars(initial_output,0).');
end