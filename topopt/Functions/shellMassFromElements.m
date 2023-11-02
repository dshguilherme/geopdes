function M = shellMassFromElements(Me, lm, x, thickness, RHO, modo)
    ndof = max(max(lm)); 
    M = zeros(ndof);
    [sz1 sz2] = size(lm);
    switch modo
        case "SIMP"
            idx = find(x > 0.1);
            nidx = setdiff(1:length(x),idx);
            density = x;
            density(idx) = 1e-3 +(RHO-1e-3)*density(idx);
            density(nidx) = 1e-3 +((RHO-1e-3))*(density(nidx).^9);
        case "Continuous"
            density = RHO*ones(size(thickness));
    end
    for i=1:sz1
        t = thickness(i);
        m = density(i)*t*Me(i,:);
        m = reshape(m,sz2,sz2);
        idx = lm(i,:);
        M(idx, idx) = M(idx,idx) +m;
    end
    M = sparse(M);
end