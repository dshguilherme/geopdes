function Ve = elementVolume(element, sp, msh)
    coeff = @(x,y) [ones(size(x)); zeros(size(x))];
    iel = element;
    msh_col = msh_evaluate_element_list (msh, iel);
    sp_col = sp_evaluate_element_list (sp, msh_col, 'value', true, ...
                               'gradient', false);
   for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
Fe = op_f_v (sp_col, msh_col, coeff(x{:}));
Ve = sum(Fe);
end