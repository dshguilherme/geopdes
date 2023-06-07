if outeriter < 0.5
    [f0val, df0dx, fval, dfdx] = cantilever2(xval, msh, sp, Ke, Me, F, Ve, ...
        lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0, h, Hs, eta);
    innerit = 0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval';
end

% Outer iterations:
kktnorm = kkttol+10;
outit = 0;
while kktnorm > kkttol && outit < maxoutit
    outit = outit+1;
    outeriter = outeriter+1;
    % low, upp, raa0 and raa
    [low, upp, raa0, raa] = asymp(outeriter, n, xval, xold1, xold2, xmin, ...
        xmax, low, upp, raa0, raa, raa0eps, raaeps, df0dx, dfdx);
    % Filter densities
    
    % GCMMA subprolem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
    gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
               raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
    xval = reshape(xval, nsub);
    xval = conv2(xval,h,'same')./Hs;
    xval = xval(:);
    % Calculate function values w/o gradients
    [f0valnew, fvalnew, ~, ~] = cantilever1(xmma, msh, sp, Ke, Me, F, Ve, ...
        lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0);
    
    % Check if approximation is conservative
    [conserv] = concheck(m, epsimin, f0app, f0valnew, fapp, fvalnew);
    
    % While hte approximations are non-conservative (conserv = 0),
    % Repeated inner iterations are made:
    innerit=0;
    if conserv == 0
        while conserv == 0 & innerit <= 15
            innerit = innerit+1;
            % New values on the parameters raa0 and raa are calculated:
            [raa0, raa] = ...
                raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                f0app, fapp, raa0, raa, raa0eps, raaeps, epsimin);
            % The GCMMA subproblem is solved with these new raa0 and raa:
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
            gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
               raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
           % Calculate function values w/o gradients of objective and
           % constraint functions at point xmma
           xval = reshape(xval, nsub);
           xval = conv2(xval,h,'same')./Hs;
           xval = xval(:);
           [f0valnew,fvalnew, ~, ~] = cantilever1(xmma, msh, sp, Ke, Me, F, Ve, ...
               lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0);
           % Check conservative
           [conserv] = concheck(m, epsimin,f0app,f0valnew,fapp,fvalnew);
        end
    end
    % No more inner iters. Update vectors:
    xold2 = xold1;
    xold1 = xval;
    xval = xmma;
    % Calculate function values and gradients at xval
    xval = reshape(xval, nsub);
    xval = conv2(xval,h,'same')./Hs;
    xval = xval(:);
    [f0val, df0dx, fval, dfdx] = cantilever2(xval, msh, sp, Ke, Me, F, Ve, ...
        lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0, h, Hs, eta);
    % Calculate the residual vector of the KKT conditions:
    [residu, kktnorm, residumax] = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
  outvector1 = [outeriter innerit f0val fval'];
  outvector2 = xval';
end
