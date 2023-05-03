function [Ke, dofs] = elementaryStiffness(element, sp1, sp2, msh, lambda, mu)
    iel = element;
    msh_col = msh_evaluate_element_list (msh, iel);
    sp1_col = sp_evaluate_element_list (sp1, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);
    sp2_col = sp_evaluate_element_list (sp2, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);
   for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
Ke = op_su_ev (sp1_col, sp2_col, msh_col, lambda (x{:}), mu (x{:}));
nel_dof = sqrt(nnz(Ke));

Ke = reshape(nonzeros(Ke), nel_dof, nel_dof);
dofs = sp1_col.connectivity;
end