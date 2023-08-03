function K = Elasticity_WQ_Bulk(msh, space, geometry, YOUNG, POISSON)
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

C = C_ijkl(E,jac,D,nsd,aux_size);
K = Stiff_fast_bulk(space,Quad_rules,Connectivity,BSval,BSder,C);
end

function Stiff_matrix = Stiff_fast_bulk(space,Quad_rules,Connectivity,BSval,BSder,aux_val)
dim = numel(space.scalar_spaces);
d = ndims(aux_val)-2;
N_dof = space.scalar_spaces{1}.ndof;

nonzeros = prod(arrayfun(@(i)sum(Connectivity(i).num_neigh),1:d));
cols = zeros(1,nonzeros); rows = cols; values = cols;
[values11, values21, values22] = deal(values);
if dim == 3
    [values31, values32, values33] = deal(values);
end
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
            C11 = aux_val{1,1};
            C21 = aux_val{2,1};
            C22 = aux_val{2,2};
            C11 = C11(points{:},k1,k2);
            C21 = C21(points{:},k1,k2);
            C22 = C22(points{:},k1,k2);
            if dim == 3
                C31 = aux_val{3,1};
                C32 = aux_val{3,2};
                C33 = aux_val{3,3};
                C31 = C31(points{:},k1,k2);
                C32 = C32(points{:},k1,k2);
                C33 = C33(points{:},k1,k2);
            end

% 			C = aux_val2(points{:},k1,k2); % coefficient tensor
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
                    C11 = tprod__(B,C11,ll);
                    C21 = tprod__(B,C21,ll);
                    C22 = tprod__(B,C22,ll);
                    if dim == 3
                        C31 = tprod__(B,C31,ll);
                        C32 = tprod__(B,C32,ll);
                        C33 = tprod__(B,C33,ll);
                    end
% 					C = tprod__(B,C,ll);
				elseif (d == 2)
					if (ll == 1)
                        C11 = B*C11;
                        C21 = B*C21;
                        C22 = B*C22;
                        if dim == 3
                            C31 = B*C31;
                            C32 = B*C32;
                            C33 = B*C33;
                        end
% 						C = B*C;
					elseif (ll == 2)
                        C11 = C11*B';
                        C21 = C21*B';
                        C22 = C22*B';
                        if dim == 3
                            C31 = C31*B';
                            C32 = C32*B';
                            C33 = C33*B';
                        end
% 						C = C*B';
					end
				end
            end
            values11(ncounter+1:ncounter+i_nonzeros) = values11(ncounter+1:ncounter+i_nonzeros) +C11(:)';
            values21(ncounter+1:ncounter+i_nonzeros) = values21(ncounter+1:ncounter+i_nonzeros) +C21(:)';
            values22(ncounter+1:ncounter+i_nonzeros) = values22(ncounter+1:ncounter+i_nonzeros) +C22(:)';
                if dim == 3
                    values31(ncounter+1:ncounter+i_nonzeros) = values31(ncounter+1:ncounter+i_nonzeros) +C31(:)';
                    values32(ncounter+1:ncounter+i_nonzeros) = values32(ncounter+1:ncounter+i_nonzeros) +C32(:)';
                    values33(ncounter+1:ncounter+i_nonzeros) = values33(ncounter+1:ncounter+i_nonzeros) +C33(:)';
                end
% 			values(ncounter+1:ncounter+i_nonzeros) = values(ncounter+1:ncounter+i_nonzeros) + C(:)';
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
k11 = sparse(rows, cols, values11, N_dof, N_dof);
k21 = sparse(rows, cols, values21, N_dof, N_dof);
k22 = sparse(rows, cols, values22, N_dof, N_dof);
if dim == 3
    k31 = sparse(rows, cols, values31, N_dof, N_dof);
    k32 = sparse(rows, cols, values32, N_dof, N_dof);
    k33 = sparse(rows, cols, values33, N_dof, N_dof);
end
% Stiff_matrix = sparse (rows, cols, values, N_dof, N_dof); % assemble the matrix
if dim == 3
    Stiff_matrix = [k11, k21', k31'; k21, k22, k32'; k31, k32, k33];
elseif dim == 2
    Stiff_matrix = [k11, k21'; k21, k22];
end
end