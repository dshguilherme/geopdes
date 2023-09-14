function ke = op_gradgradu_gradgradv_elements(spu, spv, msh, coeff)  
  der2u = reshape (spu.shape_function_hessians, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  der2v = reshape (spv.shape_function_hessians, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
  ndir = size (der2u, 2);
  ke = zeros(msh.nel,spu.nsh_max*spv.nsh_max);
  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      der2u_iel = reshape (der2u(:,:,:,:,iel), spu.ncomp*ndir, msh.nqn, 1, spu.nsh_max);
      der2v_iel = reshape (der2v(:,:,:,:,iel), spv.ncomp*ndir, msh.nqn, spv.nsh_max, 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
      
      jacdet_der2u = bsxfun (@times, jacdet_iel, der2u_iel);
      tmp1 = sum (bsxfun (@times, jacdet_der2u, der2v_iel), 1);
      elementary_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);

      ke(iel,:) = elementary_values(:)';

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradgradu_gradgradv_elements: singular map in element number %d', iel)
    end
  end

end