function M = shellMassFromElements(Me, lm, thickness)
    ndof = max(max(lm));
    M = zeros(ndof);
    [sz1 sz2] = size(lm);
    for i=1:sz1
        t = thickness(i);
        m = t*Me(i,:);
        m = reshape(m,sz2,sz2);
        idx = lm(i,:);
        M(idx, idx) = M(idx,idx) +m;
    end
    M = sparse(M);
end