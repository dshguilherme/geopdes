% OP_KL_SHELLS: assemble the Kirchhoff-Love stiffness matrix K.
%
%   mat = op_KL_shells (spu, spv, msh, E_coeff, nu_coeff, t_coeff);
%
% INPUT:
%
%  spu:   structure representing the space of trial functions (see sp_vector/sp_evaluate_col)
%  spv:   structure representing the space of test functions  (see sp_vector/sp_evaluate_col)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%  E_coeff:  coefficients for the Young's modulus
%  nu_coeff: coefficients for the Poisson's ratio
%  t_coeff:  thickness of the shell
%
% OUTPUT:
%
%  mat:    assembled stiffness matrix
% 
% Copyright (C) 2018, 2019 Pablo Antolin, Luca Coradello
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

function [Be, Se] = op_KL_shells_elements_parallel (ref_sp_u, ref_sp_v, msh, E_coeff, nu_coeff)

   wq = msh.jacdet .* msh.quad_weights;
   t_coeff = 1;
   [bending_stress, Kappa] = op_KL_bending_stress (ref_sp_u, ref_sp_v, msh, E_coeff .* wq, nu_coeff, t_coeff);
   [membrane_stress, Epsilon] = op_KL_membrane_stress (ref_sp_u, ref_sp_v, msh, E_coeff .* wq, nu_coeff, t_coeff);
   
   bending_stress = reshape(bending_stress, [], msh.nqn, ref_sp_u.nsh_max, msh.nel);
   membrane_stress = reshape(membrane_stress, [], msh.nqn, ref_sp_u.nsh_max, msh.nel);

   Kappa = reshape(Kappa, [], msh.nqn, ref_sp_v.nsh_max, msh.nel);
   Epsilon = reshape(Epsilon, [], msh.nqn, ref_sp_v.nsh_max, msh.nel);
  
   Be = V;
   Se = V;
   ndof_u = ref_sp_u.nsh_max;
   ndof_v = ref_sp_v.nsh_max;
   parfor iel = 1:msh.nel
       bending_stress_iel = reshape(bending_stress (:, :, :, iel), [], msh.nqn, 1, ndof_u);
       membrane_stress_iel = reshape(membrane_stress (:, :, :, iel), [], msh.nqn, 1, ndof_u);
       
       Kappa_iel = reshape(Kappa (:, :, :, iel), [], msh.nqn, ndof_v, 1);
       Epsilon_iel = reshape(Epsilon (:, :, :, iel), [], msh.nqn, ndof_v, 1);
       
       tmp1 = sum (bsxfun (@times, bending_stress_iel, Kappa_iel), 1);
       tmp2 = sum (bsxfun (@times, membrane_stress_iel, Epsilon_iel), 1);
       
       B = reshape (sum (tmp1, 2), ndof_v, ndof_u);
       S = reshape (sum (tmp2, 2), ndof_v, ndof_u);
       Be(iel,:) = B(:);
       Se(iel,:) = S(:);
   end
end

