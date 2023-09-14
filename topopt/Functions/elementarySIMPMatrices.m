function [ke, me, dofs] = elementarySIMPMatrices(space, malha)
msh = msh_precompute(malha);
sp = sp_precompute(space,msh, 'value', true,'gradient', true, 'divergence', true);
spu = sp;
spv = sp;
for idim = 1:msh.rdim
      x{idim} = reshape (msh.geo_map(idim,:,:), msh.nqn, msh.nel);
end
lambda = @(x, y) ((0.3)/((1.3)*(0.4)) * ones (size (x)));
mu = @(x, y) (1/(2*(1.3)) * ones (size (x)));
coeff = @(x,y) ones(size(x));

shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

ndir = size (gradu, 2);

jacdet_weights = msh.jacdet .* msh.quad_weights;
jacdet_weights_mu = jacdet_weights .* mu(x{:});
jacdet_weights_lambda = jacdet_weights .* lambda(x{:});
jacdet_weights = msh.jacdet .* msh.quad_weights.*coeff(x{:});

ke = zeros(msh.nel,spu.nsh_max*spu.nsh_max);
me = ke;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      shpu_iel = reshape (shpu(:, :, :, iel), spu.ncomp, msh.nqn, 1, spu.nsh_max);
      shpv_iel = reshape (shpv(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max, 1);  
      
      gradu_iel = reshape (gradu(:,:,:,:,iel), spu.ncomp, ndir, msh.nqn, spu.nsh_max);
      epsu_iel = (gradu_iel + permute (gradu_iel, [2 1 3 4]))/2;
      epsu_iel = reshape (epsu_iel, [spu.ncomp*ndir, msh.nqn, 1, spu.nsh_max]);
      
      gradv_iel = reshape (gradv(:,:,:,:,iel), spv.ncomp, ndir, msh.nqn, spv.nsh_max);
      epsv_iel = (gradv_iel + permute (gradv_iel, [2 1 3 4]))/2;
      epsv_iel = reshape (epsv_iel, [spv.ncomp*ndir, msh.nqn, spv.nsh_max, 1]);

      divu_iel = reshape (spu.shape_function_divs(:,:,iel), [msh.nqn, 1, spu.nsh_max]);
      divu_iel = repmat (divu_iel, [1, spv.nsh_max, 1]);
      divv_iel = reshape (spv.shape_function_divs(:,:,iel), [msh.nqn, spv.nsh_max, 1]);
      divv_iel = repmat (divv_iel, [1, 1, spu.nsh_max]);

      jacdet_mu_iel = reshape (jacdet_weights_mu(:,iel), [1,msh.nqn,1,1]);
      jacdet_lambda_iel = reshape (jacdet_weights_lambda(:,iel), [msh.nqn,1,1]);
      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
        
      jacdet_epsu = bsxfun (@times, jacdet_mu_iel, epsu_iel);
      jacdet_shpu = bsxfun (@times, jacdet_iel, shpu_iel);
      tmp1 = sum (bsxfun (@times, jacdet_shpu, shpv_iel), 1);
      
      aux_val1 = 2 * sum (bsxfun (@times, jacdet_epsu, epsv_iel), 1);
      aux_val2 = bsxfun (@times, jacdet_lambda_iel, divu_iel .* divv_iel);

      stiffness_values = reshape (sum (aux_val1, 2), spv.nsh_max, spu.nsh_max) + ...
          reshape (sum(aux_val2, 1), spv.nsh_max, spu.nsh_max);
      mass_values = reshape (sum (tmp1, 2), spv.nsh_max, spu.nsh_max);
      
      ke(iel,:) = stiffness_values(:)';
      me(iel,:) = mass_values(:)';
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end
dofs = sp.connectivity;
end