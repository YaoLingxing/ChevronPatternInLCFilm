function Ax = atv(x,t,dt,A,dfdy,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Matrix vector product
%
% Input:
%   x: input vector.
%   t: points where  solutions are located
%   dt: step size.
%   A: integration matrix.
%   dfdy: the Jacobian matrix.
%   n: number of Gaussian points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    x=reshape(x,n+1,m);
    dum=dt*A*x;

    for i=1:n+1
      dum(i,:)=dum(i,:)*dfdy(:,:,i); 
    end

    rhs=x-dum;
    delta=update(rhs,t,dfdy,n,m);  % Apply preconditioner.

    Ax=delta(:);

return
