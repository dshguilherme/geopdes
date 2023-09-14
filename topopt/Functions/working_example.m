% Physical domain, defined as NURBS map
clear problem_data
radius = 25.0;
theta = 2.0 / 9.0 * pi;

nrb = nrbreverse (nrbcirc(25, [0 0 0], pi/2 - theta, pi/2 + theta));
srf = nrbextrude (nrbtform(nrb,vecrotx(pi/2)), [0 50 0]);
problem_data.geo_name = srf;

% Type of boundary conditions for each side of the domain
% Only homogeneous Dirichlet conditions have been implemented so far.
problem_data.drchlt_sides = [3 4];
problem_data.drchlt_components = {[1 3] [1 3]};

% Physical parameters
E = 4.32e8;
nu = 0.0;
thickness = 0.25;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) -90*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

% Discretization parameters
deg = 3;

clear method_data
method_data.degree     = [deg deg];
method_data.regularity = [deg-1 deg-1];
method_data.nsub       = [16 16];
method_data.nquad      = [deg+1 deg+1];

data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

if (any(degree <= 1) || any(regularity == 0))
  error ('The degree must be at least two, and the regularity at least C^1')
end

geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, ~, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);

% Construct space structure
sp_scalar = sp_nurbs (geometry.nurbs, msh);
scalar_spaces = repmat ({sp_scalar}, 1, msh.rdim);
space = sp_vector (scalar_spaces, msh);

K_tp = op_KL_shells_tp (space, space, msh, E_coeff, nu_coeff, thickness);

% Precompute mesh, space and x
K_ele = sparse(zeros(size(K)));
for i=1:msh.nel
        msh1 = msh_evaluate_element_list(msh,i);
        sp1 = sp_evaluate_element_list_param(space, msh1, 'value', true,'gradient', true, 'hessian', true);
    for idim = 1:msh1.rdim
          x{idim} = reshape (msh1.geo_map(idim,:,:), msh1.nqn, msh1.nel);
    end

% Calculate K element by element
K_ele = K_ele+ op_KL_shells(sp1, sp1, msh1, E_coeff(x{:}), nu_coeff(x{:}), 0.25);
end
% Frobenius norm of the difference between matrices
norm(K_tp - K_ele,'fro')
