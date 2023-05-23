function [Me, dofs] = elementaryMass(element, sp1, sp2, msh, rho)
    iel = element;
    msh_col = msh_evaluate_element_list (msh, iel);
    sp1_col = sp_evaluate_element_list (sp1, msh_col, 'value', true, ...
                               'gradient', false, 'divergence', true);
    sp2_col = sp_evaluate_element_list (sp2, msh_col, 'value', true, ...
                               'gradient', false, 'divergence', true);
   for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
Me = op_u_v(sp1_col, sp2_col, msh_col, rho(x{:}));
nel_dof = sp1.nsh_max;
idx = sp1_col.connectivity;

Me = reshape(Me(idx, idx), nel_dof, nel_dof);
dofs = sp1_col.connectivity;
end