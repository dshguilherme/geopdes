function pts = generatePoints(geometry)
nurbs = geometry.nurbs;
knots = nurbs.knots;
u = cell(length(knots),1);
for i=1:length(knots)
    n = length(knots{i}) -nurbs.order(i);
    u{i} = unique(conv(knots{1},[0.5 0.5],'valid')); % Midpoints of knots to represent the fun max
end
[u{1:1:end}] = ndgrid(u{1:1:end});
nn = numel(u);
pts = reshape(cat(nn,u{:}),[],nn); % Make a combination of each u,v point;
pts = nrbeval(nurbs,pts);
pts = pts';
end