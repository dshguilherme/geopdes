function de = fastStiffnessSentivities(ke,lm,lambda,u,xPhys, YOUNG)
[sz1, sz2] = size(ke);
aux_size = [sqrt(sz2), sqrt(sz2)];
SIMP_penalty = (3*xPhys(:).^2)*(YOUNG -1e-3);
de = arrayfun(@(i) (lambda(lm(:,i)).')*reshape(ke(i,:),aux_size)*u(lm(:,i)),1:sz1);
de = SIMP_penalty.*de.';
end