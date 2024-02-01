function M = shellMassFromElements_parallel(Me, lm, thickness, RHO)
    M = zeros(max(max(lm)));
    Me = RHO*thickness.*Me;
    [sz1, sz2] = size(lm);
    for i=1:sz1
        idx = lm(i,:);
        M(idx,idx) = M(idx,idx) +reshape(Me(i,:),sz2,sz2);
    end
    M = sparse(M);
end
