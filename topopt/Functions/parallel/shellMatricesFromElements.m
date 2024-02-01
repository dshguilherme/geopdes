function [K, M] = shellMatricesFromElements(Bke, Ske, Me, lm, thickness, YOUNG,RHO)
[sz1, sz2] = size(lm);
Ke = YOUNG*((thickness.^3).*Bke +thickness.*Ske);
Me = RHO*thickness.*Me;
idx = uint16(zeros(sz1*sz2*sz2,2));
kvals = zeros(sz1*sz2*sz2,1);
mvals = kvals;
for i=1:sz1
    [id1, id2] = ndgrid(lm(i,:),lm(i,:));
    idx(1+(i-1)*sz2*sz2:i*sz2*sz2,:) = uint16([id1(:), id2(:)]);
    kvals(1+(i-1)*sz2*sz2:i*sz2*sz2) = Ke(i,:);
    mvals(1+(i-1)*sz2*sz2:i*sz2*sz2) = Me(i,:);
end
K = sparse(idx(:,1),idx(:,2),kvals);
M = sparse(idx(:,1),idx(:,2),mvals);
end