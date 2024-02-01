function [K, M] = shellStiffnessFromElements_parallel(Bke, Ske, lm, thickness, YOUNG,RHO)
    K = zeros(max(max(lm)));
    Ke = YOUNG*((thickness.^3).*Bke +thickness.*Ske);
    [sz1, sz2] = size(lm);
    for i=1:sz1
        idx = lm(i,:);
        K(idx,idx) = K(idx,idx) +reshape(Ke(i,:),sz2,sz2);
    end
    K = sparse(K);
 end