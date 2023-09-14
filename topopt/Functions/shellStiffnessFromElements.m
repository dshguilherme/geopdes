function K = shellStiffnessFromElements(Bke, Ske, lm, thickness)
    ndof = max(max(lm));
    K = zeros(ndof);
    [sz1 sz2] = size(lm);
    for i=1:sz1
        t = thickness(i);
        Ke = (t^3)*Bke(i,:) +t*Ske(i,:);
        Ke = reshape(Ke,sz2,sz2);
        idx = lm(i,:);
        K(idx, idx) = K(idx,idx) +Ke;
    end
    K = sparse(K);
end

