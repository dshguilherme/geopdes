function initial_output = firstStepPressureVessel(parameters, problem_data)
data_names = fieldnames(parameters);
for iopt = 1:numel(data_names)
    eval([data_names{iopt} '=parameters.(data_names{iopt});']);
end
clear method_data
method_data.degree = [degree degree];
method_data.regularity = [degree-1 degree-1];
method_data.nsub = nsub;
method_data.nquad = [degree+1 degree+1];

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
F = pressure*op_pn_v(sp1,msh1,problem_data.p(x{:})); % Pressure

%% get MS Relations
% master_boundary = sp.boundary(1).dofs;
% slave_boundary = sp.boundary(2).dofs;

% ms(:,1) = master_boundary(:);
% ms(:,2) = slave_boundary(:);
% degen_dofs = [sp.boundary(3).dofs sp.boundary(4).dofs];
% master = degen_dofs(1:sp.boundary(3).ndof_dir:end);
% slaves = setdiff(degen_dofs,master);
% 
% side_slaves = sp.boundary(2).dofs;
% side_master = sp.boundary(1).dofs;
% % 
% [~, ~, ib] = intersect(slaves,side_slaves);
% side_slaves(ib) = [];
% side_master(ib) = [];
% 
% ms = zeros(length(slaves),2);
% 
% ms(:,2) = slaves;
% ms(:,1) = sort(repmat(master',length(slaves)/6,1));
% ms = [ms; side_master(:), side_slaves(:)];
% ms = sort(ms);
% ms = zeros(size(ms));
% ms = [side_master(:), side_slaves(:)];

%% Generate Bending Strip
% 4 Bending strips, one in each 90 degree line of the cylinder
domain = geometry.nurbs;
% Find the strips of control points
% number = domain.number(1);
% strips_idx = 1+(number-1)/4:(number-1)/4:number;

% Generate the Bending Strips
% PP1 = [0 0 0];
% PP2 = [0 0 0];
% line = nrbline(PP1,PP2);
% [~, sz1, sz2] = size(domain.coefs);
% interval = linspace(0,1,sz2);
% strip = nrbkntins(line, interval(2:end-1));
% strip.number = [strip.number 3];
% strip.knots = {strip.knots, [0 0 0 1 1 1]};
% strip.order = [2 3];

% strip1 = strip;
% coefs = zeros(4,length(interval),3);
% coefs(:,:,1) = domain.coefs(:,strips_idx(1)-1,:);
% coefs(:,:,2) = domain.coefs(:,strips_idx(1),:);
% coefs(:,:,3) = domain.coefs(:,strips_idx(1)+1,:);
% strip1.coefs = coefs;
% 
% strip2 = strip;
% coefs(:,:,1) = domain.coefs(:,strips_idx(2)-1,:);
% coefs(:,:,2) = domain.coefs(:,strips_idx(2),:);
% coefs(:,:,3) = domain.coefs(:,strips_idx(2)+1,:);
% strip2.coefs = coefs;
% 
% strip3 = strip;
% coefs(:,:,1) = domain.coefs(:,2,:);
% coefs(:,:,2) = domain.coefs(:,1,:);
% coefs(:,:,3) = domain.coefs(:,end-1,:);
% strip3.coefs = coefs;
% 
% strip4 = strip;
% coefs(:,:,1) = domain.coefs(:,3,:);
% coefs(:,:,2) = domain.coefs(:,2,:);
% coefs(:,:,3) = domain.coefs(:,1,:);
% strip4.coefs = coefs;
% 
% strip5 = strip;
% coefs(:,:,1) = domain.coefs(:,end-2,:);
% coefs(:,:,2) = domain.coefs(:,end-1,:);
% coefs(:,:,3) = domain.coefs(:,end,:);
% strip5.coefs = coefs;

% [bK1, ssp1] = bending_strip_stiffness(strip1);
% [bK2, ssp2] = bending_strip_stiffness(strip2);
% [bK3, ssp3] = bending_strip_stiffness(strip3);
% [bK4, ssp4] = bending_strip_stiffness(strip4);
% [bK5, ssp5] = bending_strip_stiffness(strip5);

% bK = cell(2,1);
% bK{1} = bK1;
% bK{2} = bK2;
% bK{1} = bK3;
% bK{2} = bK4;
% bK{3} = bK5;

% DOF mapping from strip to global dof
% strip_dofs = cell(3,1);
% strip1_dofs = sp.boundary(1).dofs'+(number-1)/4;
% strip_dofs{1} = [strip1_dofs-1,strip1_dofs,strip1_dofs+1];
% 
% strip2_dofs = sp.boundary(1).dofs'+2*((number-1)/4);
% strip_dofs{2} = [strip2_dofs-1,strip2_dofs,strip2_dofs+1];
% % 
% strip3_dofs = sp.boundary(1).dofs'+1;
% strip_dofs{1} = [strip3_dofs strip3_dofs-1 sp.boundary(2).dofs'-1];

% strip4_dofs = sp.boundary(1).dofs'+4*(number-1)/4;
% strip_dofs{4} = [strip4_dofs-1,strip4_dofs,sp.boundary(1).dofs'+1];
% strip4_dofs = sp.boundary(1).dofs';
% strip_dofs{2} = [strip4_dofs+2 strip4_dofs+1 strip4_dofs];
% % 
% strip5_dofs = sp.boundary(2).dofs';
% strip_dofs{3} = [strip5_dofs-2 strip5_dofs-1 strip5_dofs];
%% Set No displacement conditions
% cpoints = reshape(geometry.nurbs.coefs,4,numel(geometry.nurbs.coefs)/4);
% Find all controlpoints with y below -A/2
% idx_y = find(cpoints(2,:) <= -A/2);

% Find all controlpoints B/4 away from the origin, with a 10% margin
% idx_x1 = find(cpoints(1,:) <= B/4+B/40);
% idx_x2 = find(cpoints(1,:) >= B/4-B/40);
% idx_x = intersect(idx_x1,idx_x2);
% 
% idx_x3 = find(cpoints(1,:) >= -B/4-B/40);
% idx_x4 = find(cpoints(1,:) <= -B/4+B/40);
% neg_idx_x = intersect(idx_x3,idx_x4);
% 
% idx_x = union(idx_x,neg_idx_x);

% Intersection of x and y points
% idx = intersect(idx_x,idx_y);

% dr_dofs = [idx, idx+sp.cumsum_ndof(2), idx+sp.cumsum_ndof(3)];
% free_dofs = setdiff(1:sp.ndof,dr_dofs);
% free_dofs = setdiff(free_dofs,ms(:,2));
dr_dofs = [];
for i=1:length(problem_data.drchlt_sides)
    bdofs = sp.boundary(problem_data.drchlt_sides(i)).dofs;
    dr_dofs = union(dr_dofs,bdofs);
end
dr_dofs = unique(dr_dofs);
free_dofs = setdiff(1:sp.ndof,dr_dofs);
%% Initialize MMA constants
% 
tmax = max_thickness;
tmin = min_thickness;
factor = 100*(thickness-tmin)/(tmax-tmin);
dfactor = (tmax-tmin)/100;
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
for i=1:length(strip_dofs)
    sd = strip_dofs{i}(:);
    Ks(sd,sd) = Ks(sd,sd)+bK{i};
end

% for i=1:length(ms)
%     Ks(ms(i,:),ms(i,:)) = 1e10*[1 -1;-1 1];
% end
    

C = alpha_*M +beta_*Ks;
omega = 2*pi*freq;
Kd = Ks +1j*omega*C -omega*omega*M;

% Solution 
u_init = zeros(length(Kd),1);
u_init(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs); % No need for F = F - K*u(dr_dofs) because u(dr_dofs) = 0
% u_init(ms(:,2)) = u_init(ms(:,1)); % Master/slave displacements

% Objective Functions
aW_init = real(0.5*omega*omega*(u_init'*C*u_init));
velocity0 = -1j*omega*u_init;
v2_init = real(velocity0'*velocity0);

f1 = @vesselEvalShellObjectivesAndSensitivities;
f2 = @vesselEvalShellObjectives;

s = whos;
initial_output = cell2struct({s.name}.',{s.name});
eval(structvars(initial_output,0).');
end