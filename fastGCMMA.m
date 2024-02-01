function [xval, fobj, fres, x_history] = fastGCMMA(io)
    % Initialize
    eval_objective_and_constraints = io.f1;
    eval_objective = io.f2;
    
    % Load the GCMMA parameters
    xval = io.xval;
    xold1 = io.xold1;
    xold2 = io.xold2;
    xmin = io.xmin;
    xmax = io.xmax;
    low = io.low;
    upp = io.upp;
    raa0 = io.raa0;
    raa = io.raa;
    
    % Initial values for the functions
    [f0val, df0dx, fval, dfdx] = eval_objective_and_constraints(io, xval);
    change = 1;
    x_history = zeros(numel(xval),io.maxoutit);
    fobj = zeros(io.maxoutit,1);
    fres = zeros(io.maxoutit,io.m);
    % Start outer iterations
    kktnorm = io.kkttol+10;
    outit = 0;
    outeriter = io.outeriter;
    while (kktnorm > io.kkttol) & (outit < io.maxoutit)  &  (change > io.change_min)
        outit = outit+1;
        outeriter = outeriter+1;
    %     Calculate low, upp, raa0 and raa
      [low,upp,raa0,raa] = ...
      asymp(outeriter,io.n,xval,xold1,xold2,io.xmin,io.xmax,low,upp, ...
            raa0,raa,io.raa0eps,io.raaeps,df0dx,dfdx);

    % Solve GCMMA subproblem at xval
      [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
      gcmmasub(io.m,io.n,outeriter,io.epsimin,xval,io.xmin,io.xmax,low,upp, ...
               raa0,raa,f0val,df0dx,fval,dfdx,io.a0,io.a,io.c,io.d);
    
    % Calculate new function values at xmma
    [f0valnew, fvalnew] = eval_objective(io, xmma);

    % Check if it is conservative
      [conserv] = concheck(io.m,io.epsimin,f0app,f0valnew,fapp,fvalnew);

    % While the approximations are not conservative, make inner iterations
    innerit = 0;
    if conserv == 0
        while conserv == 0 & innerit <= 15
            innerit = innerit+1;
    %         New values for raa0 and raa
          [raa0,raa] = ...
          raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                    f0app,fapp,raa0,raa,io.raa0eps,io.raaeps,io.epsimin);
    %         GCMMA is solved with new raa0 and raa
          [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
          gcmmasub(io.m,io.n,outeriter,io.epsimin,xval,xmin,xmax,low,upp, ...
                   raa0,raa,f0val,df0dx,fval,dfdx,io.a0,io.a,io.c,io.d);
               
    %        Calculate the new fval and check if it is conservative
          [f0valnew, fvalnew] = eval_objective(io, xmma);
          [conserv] = concheck(io.m,io.epsimin,f0app,f0valnew,fapp,fvalnew);
        end
    end
    % End of inner iterations. Update values
    change = max(abs(xmma(:)-xval(:)));
%     change = 100*abs((f0val-f0valnew)/f0val);
%     change = max(abs(xmma(:)-xval(:)))/parameters.thickness;
    xold2 = xold1;
    xold1 = xval;
    xval = xmma;
    
    % Filter new value
    xPhys = apply_x_filter(io.filter_options, xmma);
    xPhys = io.dfactor*xPhys +io.tmin;
    % Calculate function and restrictions
    [f0val, df0dx, fval, dfdx] = eval_objective_and_constraints(io, xval);
    
    % Residual vector of the KKT conditions:
      [residu,kktnorm,residumax] = ...
      kktcheck(io.m,io.n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
               xmin,xmax,df0dx,fval,dfdx,io.a0,io.a,io.c,io.d);
           
    % Output Stuff
        x_history(:,outit) = xval; x_plot = reshape(xPhys,io.nsub); fobj(outit) = f0val; fres(outit,:) = fval; 
        fprintf(' Iteration: %3i | Objective: %3.1f | Mass: %4.1f kg | Restrictions: %3.1f | %3.1f | Change: %1.2e\n', ...
      outit, f0val, io.RHO*sum(io.Ve.*xPhys), fval(1), fval(2), change);
    grafo = nrbplot(io.geometry.nurbs,io.nsub); view(0,90); colormap(jet); grafo.CData = x_plot; colorbar; caxis([io.tmin io.tmax]); drawnow; %imagesc(rot90(x_plot)); colorbar; caxis([0 100]); axis equal;axis off;drawnow;
    end
    fobj = fobj(1:outit); fres = fres(1:outit,:); x_history = x_history(:,1:outit);
end