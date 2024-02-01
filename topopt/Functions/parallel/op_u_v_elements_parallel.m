% OP_U_V: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i).
%
%   mat = op_u_v (spu, spv, msh, coeff);
%   [rows, cols, values] = op_u_v (spu, spv, msh, coeff);
%
% INPUT:
%
%  spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%  spv:   structure representing the space of test functions  (see sp_scalar/sp_evaluate_col)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%  coeff: reaction coefficient
%
% OUTPUT:
%
%  mat:    assembled mass matrix
%  rows:   row indices of the nonzero entries
%  cols:   column indices of the nonzero entries
%  values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function me = op_u_v_elements_parallel (spu, spv, msh, coeff)
  
  shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
  
  me = zeros(msh.nel,spu.nsh_max*spv.nsh_max);
  if (~isempty (msh.jacdet))
    jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  else
    jacdet_weights = [];
  end
  
    uncomp = spu.ncomp;
    nqn = msh.nqn;
    vncomp = spv.ncomp;
    unsh_max = spu.nsh_max;
    vnsh_max = spv.nsh_max;
  
  parfor iel = 1:msh.nel
      shpu_iel = reshape (shpu(:, :, :, iel), uncomp, nqn, 1, unsh_max);
      shpv_iel = reshape (shpv(:, :, :, iel), vncomp, nqn, vnsh_max, 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,nqn,1,1]);

      jacdet_shpu = bsxfun (@times, jacdet_iel, shpu_iel);
      tmp1 = sum (bsxfun (@times, jacdet_shpu, shpv_iel), 1);
      elementary_values = reshape (sum (tmp1, 2), vnsh_max, unsh_max);

      me(iel,:) = elementary_values(:)';
  end

end


%% COPY OF THE FIRST VERSION OF THE FUNCTION (MORE UNDERSTANDABLE)
% 
% function mat = op_u_v (spu, spv, msh, coeff)
%   
%   mat = spalloc(spv.ndof, spu.ndof, 1);
%   shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, spu.nsh_max, msh.nel);
%   shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
%   
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%       mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
%       for idof = 1:spv.nsh(iel)
%         ishp = reshape (shpv(:,:,idof,iel), spv.ncomp, []);
%         for jdof = 1:spu.nsh(iel)
%           jshp = reshape (shpu(:,:,jdof,iel), spu.ncomp, []);
% % The cycle on the quadrature points is vectorized
%           %for inode = 1:msh.nqn
%             mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
%               sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
%               sum(ishp .* jshp, 1).' .* coeff(:, iel));
%           %end  
%         end
%       end
%       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
%         mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) +  mat_loc;
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_u_v: singular map in element number %d', iel)
%     end
%   end
% 
% end
