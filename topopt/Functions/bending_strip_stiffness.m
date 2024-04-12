function [K, dofs] = bending_strip_stiffness(strip)
% Construct the msh structure
rule = msh_gauss_nodes([3 3]);
zeta = {unique(strip.knots{1}), unique(strip.knots{2})};
[qn, qw] = msh_set_quad_nodes(zeta, rule);
msh = msh_cartesian(zeta, qn, qw, geo_load(strip));

space_scalar = sp_nurbs(strip, msh);
scalar_spaces = repmat({space_scalar}, 1, msh.rdim);
space = sp_vector(scalar_spaces, msh);
clear space_scalar scalar_spaces;

K = op_bending_strip_tp(space, space, msh, @(x,y,z) 1e10*ones(size(x)), @(x,y,z) zeros(size(x)), 1);
dofs = space.boundary(1).dofs;
dofs = [dofs' dofs'+1 dofs'+2];
end