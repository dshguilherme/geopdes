function K = shellStiffnessFromElements(Bke, Ske, lm, x, thickness, YOUNG, modo)
    ndof = max(max(lm));
    K = zeros(ndof);
    [sz1 sz2] = size(lm);
    switch modo
        case "SIMP"
            E = 1e-6*YOUNG.*((1e-6 -1).*x +1).^-1;
        case "Continuous"
            E = YOUNG*ones(size(thickness));
    end
   for i=1:sz1
        t = thickness(i);
        Ke = (t^3)*Bke(i,:) +t*Ske(i,:);
        Ke = E(i)*Ke;
        Ke = reshape(Ke,sz2,sz2);
        idx = lm(i,:);
        K(idx, idx) = K(idx,idx) +Ke;
    end

    K = sparse(K);
end

