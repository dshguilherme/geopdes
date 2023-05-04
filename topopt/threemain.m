%  This is the file threemain.m
%  which is used as a main program for
%  the following "three bar truss problem".
%
%  mininimize x1+x2+x3
%  subject to:
%      pi'*ui(x) <= 1, for i=1,2,3
%    0.001 <= xj <= 2, for j=1,2,3.
%
%  where
%  xj = volume of the j:th bar,
%  pi'*ui(x) = compliance for the i:th loadcase,
%
%  (Since the length of each bar is = 1,
%  xj is also the cross section area of the j:the bar.)
%  
%  Bar 1 connects the nodes 1 and 4.
%  Bar 2 connects the nodes 2 and 4.
%  Bar 3 connects the nodes 3 and 4.
%  The coordinates of node 1 are (-1,0).
%  The coordinates of node 2 are (-1/sqrt(2),-1/sqrt(2)).
%  The coordinates of node 3 are (0,-1).
%  The coordinates of node 4 are (0,0).
%  The nodes 1,2,3 are fixed, while the load vectors for
%  the different loadcases are applied at node 4.
%  The load vector p1 = (1,0)' (loadcase 1).
%  The load vector p2 = (1,1)' (loadcase 2).
%  The load vector p3 = (0,1)' (loadcase 3).
%
%  The displacement vector ui is obtained from the system
%  K(x)*ui = pi, i=1,2,3, where K(x) is the stiffness matrix
%  and pi is the load vector.
%  The stiffness matrix is given by K(x) = R*D(x)*R',
%  where D(x) is a diagonal matrix with diagonal elements x1,x2,x3.
%  The constant matrix R is given in the code below.
%  The derivatives of the functions fi(x) are then given by
%  dfi/dxj = -(rj'*ui)^2 ,
%  where rj is the j:th column of the matrix R.
%
%  The problem is written on the following form required by MMA.
%
%  minimize  x1 + x2 + x3 + z + 1000*(y1+y2+y3)
%  subject to the constraints:
%            pi'*ui(x) - 1 - yi <= 0, i=1,2,3
%                       0 <= xj <= 2, j=1,2,3
%                            yi >= 0, i=1,2,3
%                             z >= 0.
%
%
%  Written in May 1999 by
%  Krister Svanberg <krille@math.kth.se>
%  Department of Mathematics
%  SE-10044 Stockholm, Sweden.
%
itte = 0;
while itte < maxite
  iter = iter+1;
  itte = itte+1;

  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);

  xold2 = xold1;
  xold1 = xval;
  xval = xmma;

  [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = three(xval);
  outvector = [iter f0val fval' xval']'
end
