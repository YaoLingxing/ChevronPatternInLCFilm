function [y,reskeep,errs1,iters]=...
  onestep(m,t0,y0,dt,n,tc,kmax,gtol,etol,k0,A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  one time step advancement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iprob

t=dt*tc+t0;  % set up the Gaussian points in [t0,t1].
t1=dt+t0;    % t1 is the end time.

    yp0=newRHS(m,y0,t0);
    ypnew=predictor(yp0,t,m,n);
    ynew=dt*A*ypnew+ones(n+1,1)*y0; 
    rhs=newRHS(m,ynew,t)-ypnew;     % initial RHS
    yp=ypnew(n+1,:);
    y=dt*B*ypnew+y0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    errest=max(max(abs(rhs))); % estimate error (residue).
    lerrest=errest;     % used to check if residual decays or not.
    lcount=0;           % used to check if the number of nondcay error iterations.
    reskeep=[];         % gmres error.
    iters=[];           % gmres iterations.
    errs1=errest;       % right hand side error, the residual.
    count1=0;           % the number of GMRES corrections.
    error=[];           % each gmres error.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GMRES loop.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while count1<kmax & errest>etol
      count1=count1+1;

      dfdy=JacFCN(m,ynew,t);  % set the jacobian

      % Step 1: Compute initial RHS
       epsilon=update(rhs,t,dfdy,n,m);

      % Step 2: Solve the linear equation
      delta=zeros(m*(n+1),1); b=epsilon(:);
      [delta, error, total_iters] = gmres(delta, b, 'atv', [gtol, k0],t,dt,A,dfdy,n,m);
      delta=reshape(delta,n+1,m);

      % Step 3: Check error.
      ypnew=ypnew+delta;
      ynew=dt*A*ypnew+ones(n+1,1)*y0;

      yp=ypnew(n+1,:);
      y=dt*B*ypnew+y0;
      %Compute residual
      rhs=newRHS(m,ynew,t)-ypnew;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      errest=max(max(abs(rhs)));   % raw error.
      reskeep=[reskeep error];
      errs1=[errs1,errest];
      iters(count1)=total_iters;

      if errest>=0.5*lerrest
        if lcount==4, break, else lcount=lcount+1; end;
      else
        lerrest=errest;
      end;
    end


return
