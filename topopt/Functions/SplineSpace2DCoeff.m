% function rho = SplineSpace2DCoeff(nurbs,x)
%     tt = zeros(length(x),prod(size(x{1})));
%     for i=1:length(x)
%         tt(i,:) = x{i}(:);
%     end
%     p = nrbeval(nurbs,tt);
%     rho = p(3,:);
%     rho = rho(:);
% end

function rho = SplineSpace2DCoeff(nurbs,x,y)
shape = size(x);
p = nrbeval(nurbs,[x(:), y(:)]);
rho = p(3,:)';
rho = reshape(rho,shape);
end

%  arrayfun(@(x,y) SplineSpace2DCoeff(nurbs,x,y),array(:,1),array(:,2),'UniformOutput',false)
