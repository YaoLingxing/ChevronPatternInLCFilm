function [ysol, ysolfin,res,indres,errrhs,inderr,iter]=...
  ALCPSDC(iprob,m,t0,tfinal,y0,h0,n,kmax,gtol,etol,k0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INTENT (out)
%   ysolfin (row vector) solution y at final time
%   res (row vector) residues outputed by gmres
%   indres (row vector) index for res : res(indres(i-1)+1:indres(i)), 
%          residues outputs from gmres' during the i_th time step
%        
%   err (row vector) error after each gmres
%   inderr (row vector) index for err: err(inderr(i-1)+1:inderr(i)),
%          errors outputs from gmres' during the i_th time step
%   iter (row vector) iteration numbers outputed by each gmres.
% 
% INTENT (in) 
%   Most are scalars.
%   m size of the problem
%   t0 initial time
%   tfinal fintal time
%   y0 (row vector) initial value for y
%   h0 step size
%   n number of grid points used for each time step
%   kmax maximal number that the gmres is done
%   gtol tolerance for gmres
%   etol tolerance for err/res
%   k0 gmresk selecting parameter
%
% Last change: Jingfang Huang, 03/10/2005.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initiallization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A,B,tc]=tgauss(n);  %set up the integration matrix, and construct tc.

tnow=t0; dt=h0; ynow=y0;
ysol = []; % Now set up the initial values for iteration.
count=0; res=[]; errrhs=[];  iter=[]; %monitor values.
indres=[];inderr=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   the main marching scheme. No adaptive steps yet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while tnow<tfinal,
  count=count+1;

  if dt>tfinal-tnow,
    dt=tfinal-tnow;  % find the right time step. This is important for last step.
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % march one-step. main code.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [ynow,reskeep,errs1,iters]=...
    onestep(m,tnow,ynow,dt,n,tc,kmax,gtol,etol,k0,A,B);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % output the Hamiltonian system constant and the solution for the
  % geodesic Flows problem.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   if (iprob==7) then
%       figure(10);
%       plot(ynow(3),ynow(4),'*');
%       hold on;
%       figure(11)
%       plot(ynow(1),ynow(2),'*');
%       hold on;  %plot the solution at different times.
%   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % record errors for each iteration.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ysol = [ysol, ynow'];
  res=[res, reskeep];   %GMRES error
  errrhs=[errrhs,errs1]; % right hand side error.
  iter=[iter,iters];       % number of GMRES iterations.
  indres(count)=length(res); % index for res.
  inderr(count)=length(errrhs);  %index for right hand side error.
  tnow=tnow+dt;
end
count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output the final solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ysolfin=ynow;

return
