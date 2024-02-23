function initial_output = firstStepKLShell(parameters,problem_data)
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
sp1 = sp_precompute_param(sp,msh1, 'value', true,'gradient', true, 'hessian', true);
for idim = 1:msh1.rdim
      x{idim} = reshape (msh1.geo_map(idim,:,:), msh1.nqn, msh1.nel);
end
% Stiffness Matrix
[Bke, Ske] = op_KL_shells_elements(sp1,sp1,msh1,problem_data.E_coeff(x{:}), problem_data.nu_coeff(x{:}));
% Mass Matrix
Me = op_u_v_elements(sp1, sp1, msh1, RHO);
% Location matrix (connectivity)
lm = sp1.connectivity';

% Restriction
Ve = (msh1.element_size.^2)'; % Element area

% Force vector
F = op_f_v_tp (sp, msh, problem_data.f);
F = Fmag*F/sum(F);

tmax = max_thickness;
tmin = min_thickness;
factor = 100*(thickness-tmin)/(tmax -tmin);
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

%% get BC dofs
dr_dofs = [];
for iside = 1:numel(problem_data.drchlt_sides)
  side = problem_data.drchlt_sides(iside);
  components = problem_data.drchlt_components{iside};
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

%% Initial Solution
t = (tmax-tmin)*xval/100 +tmin;
t = apply_x_filter(filter_options, t);
[Ks, M] = shellMatricesFromElements(Bke, Ske, Me, lm, t, YOUNG, RHO);
if initial_guess == 1
    D = YOUNG*thickness^3/(1-POISSON*POISSON);
    const = sqrt(D/RHO/thickness)*pi*pi;
    n_eigs = ceil(sqrt(2*pi*freq*(1/const)));
    n_eigs = 2*numel(ndgrid(0:n_eigs+1,0:n_eigs+1));
    t = getInitialGuess(Ks,M,dr_dofs,n_eigs,freq,sp, geometry,nsub);
    t = reshape(t,size(xval));
    t = norm(xval)/norm(t)*t;
    xval = t;
    t = (tmax-tmin)*t/100 +tmin;
    [Ks, M] = shellMatricesFromElements(Bke, Ske, Me, lm, t, YOUNG, RHO);
end

u_init = cell(nfreq,2);
aW_init = cell(nfreq,1);
for i=1:nfreq
C = alpha_(i)*M +beta_(i)*Ks;
Kd = Ks +1j*omega(i)*C -omega(i)*omega(i)*M;
dr_values = zeros(length(dr_dofs),1);
u_init{i,1} = SolveDirichletSystem(Kd, F, dr_dofs, free_dofs, dr_values);
% us = SolveDirichletSystem(Ks, F, dr_dofs, free_dofs, dr_values);
% [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
% vals = diag(vals);
%% Objective Functions
aW_init{i} = real(0.5*omega(i)*omega(i)*(u_init{i,1}'*C*u_init{i,1}));
velocity0 = -1j*omega(i)*u_init{i,1};
% % gap0 = vals(1)-vals(2);
% % eig0 = -vals(1);
% pts = generatePoints(geometry);
% k = omega/320; % 320m/s = speed of sound in air k is the angular wavenumber
% R0 = zeros(length(pts));
% const = omega*omega*sum(Ve)*1.204/(4*pi*320); %1.204 -> density of air
% for i=1:length(pts)
%     P = pts(i,:);
%     dist = P-pts;
%     dist = vecnorm(dist')';
%     R0(i,:) = const*sin(k*dist)./(k*dist);
%     R0(i,i) = const*1;
% end
% R0 = sparse(R0);
% R0 = sparse(blkdiag(R0,R0,R0));
% V0 = real(velocity0'*R0*velocity0);
u_init{i,2} = real(velocity0'*velocity0);
end

f1 = @parallelEvalShellObjectivesAndSensitivities;
f2 = @parallelEvalShellObjectives;

s = whos;
initial_output = cell2struct({s.name}.',{s.name});
eval(structvars(initial_output,0).');

end