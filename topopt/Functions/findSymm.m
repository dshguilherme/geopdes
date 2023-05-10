function symm_dofs = findSymm(sp, msh, problem_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end

% Apply symmetry conditions
symm_dofs = [];
for iside = symm_sides
  if (~strcmpi (sp.transform, 'grad-preserving'))
    error ('The symmetry condition is only implemented for spaces with grad-preserving transform')
  end
  msh_side = msh_eval_boundary_side (msh, iside);
  for idim = 1:msh.rdim
    normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
  end

  parallel_to_axes = false;
  for ind = 1:msh.rdim
    ind2 = setdiff (1:msh.rdim, ind);
    if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
      symm_dofs = union (symm_dofs, sp.boundary(iside).dofs(sp.boundary(iside).comp_dofs{ind}));
      parallel_to_axes = true;
      break
    end
  end
  if (~parallel_to_axes)
    error ('solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
  end

end
end