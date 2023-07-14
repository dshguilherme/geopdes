function rho = SplineInterp2D(xPhys,x,y)
[m,n] = size(xPhys);
[X,Y] = meshgrid(linspace(0,1,m),linspace(0,0.5,n));
shape = size(x);
rho = interp2(X,Y,xPhys',x(:),y(:),'spline');
rho = reshape(rho,shape);
end