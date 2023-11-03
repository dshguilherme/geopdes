function df0dx = CalculateSensivities(ell, u, lm, dk)
nel = size(lm,1);
df0dx = zeros(nel,1);
for i=1:nel
    dofs = lm(i,:)';
    df0dx(i) = real(transpose((ell(dofs)))*squeeze(dk(i,:,:))*u(dofs));
end
end