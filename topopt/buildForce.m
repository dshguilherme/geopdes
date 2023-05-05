function rhs = buildForce(sp, msh, problem_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end


rhs = op_f_v_tp(sp,msh,f);
% Apply Neumann boundary conditions
for iside = nmnn_sides
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) g(varargin{:},iside);
  dofs = sp.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) + op_f_v_tp (sp.boundary(iside), msh.boundary(iside), gside);
end

% Apply pressure conditions
for iside = press_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  sp_side  = sp_eval_boundary_side (sp, msh_side);

  x = cell (msh_side.rdim, 1);
  for idim = 1:msh_side.rdim
    x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
  end
  pval = reshape (p (x{:}, iside), msh_side.nqn, msh_side.nel);

  rhs(sp_side.dofs) = rhs(sp_side.dofs) - op_pn_v (sp_side, msh_side, pval);
end


end