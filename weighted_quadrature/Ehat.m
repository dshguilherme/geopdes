function Eh_ik = Ehat(E,DF_minus_T,i,k)
[~, ~, sz3] = size(DF_minus_T);
[sz1, ~, sz2] = size(E);
Eh_ik = zeros(sz1, 1, sz3);
for idx=1:sz2
    Eh_ik = Eh_ik+E(:,i,idx).*DF_minus_T(idx,k,:);
end
end