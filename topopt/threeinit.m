%  This is the file threeinit.m
%  in which some vectors for the "three bar
%  truss problem" defined in the file threemain.m
%  are initialized.
%
%  Written in May 1999 by
%  Krister Svanberg <krille@math.kth.se>
%  Department of Mathematics
%  SE-10044 Stockholm, Sweden.
%
m = 3;
n = 3;
xval  = ones(n,1);
xold1 = xval;
xold2 = xval;
xmin  = 0.001*ones(n,1);
xmax  = 2*ones(n,1);
low   = xmin;
upp   = xmax;
c = 1000*ones(m,1);
d = zeros(m,1);
a0 = 1;
a = zeros(m,1);
iter = 0;
[f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = three(xval);
outvector = [iter f0val fval' xval']'
maxite=2;

