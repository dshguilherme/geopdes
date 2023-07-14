function coeff = test(nurbs, x)
tt = zeros(length(x),prod(size(x{1})));
for i=1:length(x)
    tt(i,:) = x{i}(:);
end
p = nrbeval(nurbs,tt);
coeff = p(3,:);
coeff = coeff(:);
end