function [xval, fobj, fres, x_history] = GCMMA(f1,f2,init_mat, filter_options)
    % Initialize
    eval_objective_and_constraints = f1;
    eval_objective = f2;
    load(init_mat);
    
    
    % Initial values for the functions
    [f0val, df0dx, fval, dfdx] = eval_objective_and_constraints(xval);
    change = 1;
    x_history = zeros(numel(xval),maxoutit);
    fobj = zeros(maxoutit,1);
    fres = zeros(maxoutit,m);
    % Start outer iterations
    kktnorm = kkttol+10;
    outit = 0;
    while (kktnorm > kkttol) & (outit < maxoutit)  &  (change > change_min)
        outit = outit+1;
        outeriter = outeriter+1;
    %     Calculate low, upp, raa0 and raa
      [low,upp,raa0,raa] = ...
      asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
            raa0,raa,raa0eps,raaeps,df0dx,dfdx);

    % Solve GCMMA subproblem at xval
      [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
      gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
               raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
    
    % Calculate new function values at xmma
    [f0valnew, fvalnew] = eval_objective(xmma);

    % Check if it is conservative
      [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);

    % While the approximations are not conservative, make inner iterations
    innerit = 0;
    if conserv == 0
        while conserv == 0 & innerit <= 15
            innerit = innerit+1;
    %         New values for raa0 and raa
          [raa0,raa] = ...
          raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                    f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);
    %         GCMMA is solved with new raa0 and raa
          [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
          gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
                   raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
               
    %        Calculate the new fval and check if it is conservative
          [f0valnew, fvalnew] = eval_objective(xmma);
          [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
        end
    end
    % End of inner iterations. Update values
%     change = max(abs(xmma(:)-xval(:)))/norm(xval);
%     change = 100*abs((f0val-f0valnew)/f0val);
    change = max(abs(xmma(:)-xval(:)))/parameters.thickness;
    xold2 = xold1;
    xold1 = xval;
    xval = xmma;
    
    % Filter new value
    xPhys = apply_x_filter(filter_options, xmma);
    % Calculate function and restrictions
    [f0val, df0dx, fval, dfdx] = eval_objective_and_constraints(xval);
    
    % Residual vector of the KKT conditions:
      [residu,kktnorm,residumax] = ...
      kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
               xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
           
    % Output Stuff
        x_history(:,outit) = xval; x_plot = reshape(xPhys,nsub); fobj(outit) = f0val; fres(outit,:) = fval; 
        fprintf(' Iteration: %3i | Objective: %1.2e | Mass: %4.1f kg | Restrictions: %2.2e | %2.2e | Change: %1.2e\n', ...
      outit, f0val, RHO*sum(Ve.*xPhys), fval(1), fval(2), change);
    grafo = nrbplot(geometry.nurbs,nsub); colormap(jet); grafo.CData = x_plot; colorbar; drawnow; %imagesc(rot90(x_plot)); colorbar; caxis([0 100]); axis equal;axis off;drawnow;
    end
    fobj = fobj(1:outit); fres = fres(1:outit,:); x_history = x_history(:,1:outit);
end