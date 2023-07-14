function rho = rhoPhys2D(xPhys,x,y)
[m,n] = size(xPhys);
x = x*n;
y = 0.5*y*m;
rho = interp2(xPhys,x,y);
end