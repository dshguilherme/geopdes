function K = Elasticity_WQ(msh, space, geometry, YOUNG, POISSON, xPhys)
for idim = 1:msh.ndim
    sp1d = space.scalar_spaces{1}.sp_univ(idim);
    Connectivity(idim).neighbors = cellfun (@(x) unique (sp1d.connectivity(:,x)).', sp1d.supp, 'UniformOutput', false);
    Connectivity(idim).num_neigh = cellfun (@numel, Connectivity(idim).neighbors);

    Quad_rules(idim) = quadrule_stiff_fast(sp1d);
 
    brk{idim} = [space.scalar_spaces{1}.knots{idim}(1), space.scalar_spaces{1}.knots{idim}(end)];
    qn{idim} = Quad_rules(idim).all_points';
end

% Generate a new GeoPDEs mesh with the new quadrature points
new_msh = msh_cartesian (brk, qn, [], geometry);
space_wq = space.constructor (new_msh);

% compute the function values and gradients at the quadrature points
for idim = 1:msh.ndim
    sp1d = space_wq.scalar_spaces{1}.sp_univ(idim);
    for ii = 1:sp1d.ndof
      BSval{idim,ii} = sp1d.shape_functions(Quad_rules(idim).ind_points{ii}, Connectivity(idim).neighbors{ii}).'; 
	  BSder{idim,ii} = sp1d.shape_function_gradients(Quad_rules(idim).ind_points{ii}, Connectivity(idim).neighbors{ii}).'; 
    end 	
end

% compute the coefficient at the quadrature points
x = cell (msh.ndim, 1);
[x{:}] = ndgrid(qn{:});

nsd = msh.ndim;
aux_size = cellfun (@numel, qn);
jac = msh.map_der(qn);

% Incidence Matrix and Stress-Strain Tensor
if nsd == 3
    E = zeros(6,3,3);
    E(:,:,1) = [1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 0 1 0];
    E(:,:,2) = [0 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 0; 1 0 0];
    E(:,:,3) =  [0 0 0; 0 0 0; 0 0 1; 0 1 0; 0 0 1; 0 0 0];

else
    E = zeros(3,2,2);
    E(:,:,1) = [1 0; 0 0; 0 1];
    E(:,:,2) = [0 0; 0 1; 1 0];
end

Y = YOUNG;
v = POISSON;
lambda = Y*v/(1+v)/(1-2*v);
mu = Y/2/(1+v);

total_size = nchoosek(nsd+2-1,2); % Size of tensor Dijkl in Voigt Notation
small_size = total_size-nsd; % Size of the shear stress matrix
D1 = lambda*ones(nsd) +2*mu*eye(nsd);
D2 = mu*eye(small_size);
D = blkdiag(D1,D2);

rho = SplineInterp(xPhys,x{:});
rho = reshape(rho,1,1,numel(rho));
SIMP_penalty = 1e-3+ (rho.^3).*(YOUNG - 1e-3);
D = bsxfun(@times,D,SIMP_penalty);

C = C_ijkl(E,jac,D,nsd,aux_size);
K11 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{1,1});
K12 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{1,2});
% K21 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{2,1});
K21 = K12';
K22 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{2,2});
if space.ncomp == 3
K13 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{1,3});
K23 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{2,3});
K33 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{3,3});
% K31 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{3,1});
K31 = K13';
% K32 = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,C{3,2});
K32 = K23';
K = [K11, K12, K13; K21, K22, K23; K31, K32, K33];
else
    K = [K11, K12; K21, K22];
end
end

function Stiff_matrix = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,aux_val)

d = ndims(aux_val)-2;
N_dof = space.scalar_spaces{1}.ndof;

nonzeros = prod(arrayfun(@(i)sum(Connectivity(i).num_neigh),1:d));
cols = zeros(1,nonzeros); rows = cols; values = cols;
ncounter = 0;

n_size = space.scalar_spaces{1}.ndof_dir;
indices = cell(1,d);
[indices{:}] = ind2sub(n_size, 1:N_dof);
indices = cell2mat(indices);  indices = reshape(indices,[N_dof d]);
points = cell(1,d); j_act = cell(1,d); len_j_act = zeros(1,d); n_index = zeros(1,d);
for ll = 1:d
    n_index(ll) = prod(n_size(1:ll-1));
end

for ii = 1:N_dof

	ind = indices(ii,:);
	for ll = 1:d
		points{ll} = Quad_rules(ll).ind_points{ind(ll)}; 
		j_act{ll} = Connectivity(ll).neighbors{ind(ll)}; 
		len_j_act(ll) = length(j_act{ll});
	end
	i_nonzeros = prod(len_j_act);
	
	for k1 = 1:d % derivative index on the test function
		for k2 = 1:d % derivative index on the trial function
			C = aux_val(points{:},k1,k2); % coefficient tensor
			for ll = d:-1:1 % sum_factorization loop
				if (k1 == k2 && ll == k1)
					Q = Quad_rules(ll).quad_weights_11{ind(ll)};
					B = BSder{ll,ind(ll)}(1:len_j_act(ll),:);
				elseif (k1 ~= k2 && ll == k1)
					Q = Quad_rules(ll).quad_weights_10{ind(ll)};
					B = BSval{ll,ind(ll)}(1:len_j_act(ll),:);				
				elseif (k1 ~= k2 && ll == k2)
					Q = Quad_rules(ll).quad_weights_01{ind(ll)};
					B = BSder{ll,ind(ll)}(1:len_j_act(ll),:);
				else
					Q = Quad_rules(ll).quad_weights_00{ind(ll)};
					B = BSval{ll,ind(ll)}(1:len_j_act(ll),:);
				end
				B = bsxfun(@times,Q,B);
				if (d == 3)
					C = tprod__(B,C,ll);
				elseif (d == 2)
					if (ll == 1)
						C = B*C;
					elseif (ll == 2)
						C = C*B';
					end
				end
			end
			values(ncounter+1:ncounter+i_nonzeros) = values(ncounter+1:ncounter+i_nonzeros) + C(:)';
		end
	end
	
	% row indices
	rows(ncounter+1:ncounter+i_nonzeros) = ii;
	
	% compute the column indices
	i_col = zeros(d,i_nonzeros);
	for ll = 1:d
		rep = len_j_act; rep(ll) = 1;
		perm = ones(1,d); perm(ll) = len_j_act(ll);
		ap = repmat(reshape(j_act{ll}',perm),rep);
		i_col(ll,:) = ap(:)';      
	end
	cols(ncounter+1:ncounter+i_nonzeros) = 1 + n_index*(i_col-1);
% 	ap = cell(1,d);
% 	[ap{:}] = ndgrid(j_act{:});
% 	for k = d:-1:2
% 		cols(ncounter+1:ncounter+i_nonzeros) = cols(ncounter+1:ncounter+i_nonzeros) + (ap{k}(:)' - 1)*(n_index(k));
% 	end
% 	cols(ncounter+1:ncounter+i_nonzeros) = cols(ncounter+1:ncounter+i_nonzeros) + ap{1}(:)';
    
    ncounter = ncounter + i_nonzeros;
   
end

Stiff_matrix = sparse (rows, cols, values, N_dof, N_dof); % assemble the matrix

end