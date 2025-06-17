function delta=update(rhs,t,dfdy,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with the input rhs, returns
% (1-h*tilde{A}*alpha)^(-1)*rhs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E=eye(m);
    delta=zeros(n+1,m);
    sdelta=zeros(1,m);

    for i=1:n,
      dum=dfdy(:,:,i+1);
      delta(i+1,:)=(rhs(i+1,:)+sdelta*dum)/(E-(t(i+1)-t(i))*dum);
      sdelta=sdelta+(t(i+1)-t(i))*delta(i+1,:);
    end
return
