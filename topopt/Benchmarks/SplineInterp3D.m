function rho = SplineInterp3D(xPhys,x,y,z)
[m,n,p] = size(xPhys);
[X,Y,Z] = meshgrid(linspace(0,1,m),linspace(0,0.5,n),linspace(0,0.1,p));
shape = size(x);
xPhys = permute(xPhys,[2 1 3]);
rho = interp3(X,Y,Z,xPhys,x(:),y(:),z(:),'spline');
rho = reshape(rho,shape);
end