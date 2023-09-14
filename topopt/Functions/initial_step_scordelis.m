function initial_step_scordelis(parameters, problem_data)
data_names = fieldnames(parameters);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= parameters.(data_names{iopt});']);
end
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

%% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data,method_data);

%% Pre-calculate element volume, stiffness, mass, force
msh1 = msh_precompute(msh);
sp1 = sp_precompute(sp,msh1, 'value', true,'gradient', true, 'hessian', true);
for idim = 1:msh1.rdim
      x{idim} = reshape (msh1.geo_map(idim,:,:), msh1.nqn, msh1.nel);
end
% Stiffness Matrix
Ke = op_KL_shell_elements(sp1,sp1,msh1,problem_data.E_coeff, problem_data.nu_coeff);
% Mass Matrix
Me = op_u_v_elements(sp1, sp1, msh1, problem_data.density);
% Location matrix (connectivity)
lm = sp1.connectivity';

% Force vector
F = op_f_v_tp (sp, msh, problem_data.f);


%% Initialize MMA constants
m = 1; % Number of Restrictions
n = msh.nel; % Number of variables
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = vol_frac*eeen;
xold1   = xval;
xold2   = xval;
xmin    = 0*eeen;
xmax    = eeen;
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

%% get BC dofs
dr_dofs = [];
for iside = 1:numel(drchlt_sides)
  side = drchlt_sides(iside);
  if (~exist('drchlt_components','var'))
    components = 1:3;
  else
    components = drchlt_components{iside};
  end
  for icomp = components
    dr_dofs = union (dr_dofs, sp.boundary(side).dofs(sp.boundary(side).comp_dofs{icomp}));
  end
end

free_dofs = setdiff(1:sp.ndof, dr_dofs);


%% Filter
[h, Hs] = density_filter(rmin, nsub);
filter_options.type = philter;
filter_options.h = h;
filter_options.Hs = Hs;
filter_options.subshape = nsub;

%% Objective Function
Cs0 = 1;
aW0 = 1;
tmp = objective_function;
objective_function = "Initial";

f1 = @EvalShellObjectivesAndSensitivities;
f2 = @EvalShellObjectives;

%% Initial Functions
save('init_scordelis.mat')
[f0, ~] = f2(xval);
Cs0 = f0(1);
aW0 = f0(2);
objective_function = tmp;
save('init_scordelis.mat');
end