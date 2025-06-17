function dfdy=JacFCN(m,y,t)

% m is dim of the problem
% t is the time to evaluate the jacobian
% y(1:m,1:length(t)) stores y values
% dfdy(j,i,k) is df(i)/dy(j) at time t(k)

global iprob mypar;
if (iprob==7) % alpha1 is not 0  nn=length(t);
  nn=length(t);
  eps=1e-6;
  gp0=1/gp(0,mypar);
  coe=mypar.ac*mypar.Er*0.25*gp0;
  for j=1:nn
      atemp=zeros(m,m);
      p1=y(j,1);
      p2=y(j,2);
      %y3=y(j,3);
      %y4=y(j,4);
      gp0=1/gp(p1,mypar);
      atemp(1,1)=0; atemp(1,2)=1;
      %atemp(2,1)=2*p2*p2*cos(2*p1)+coe*(2*cos(2*p1)*(mypar.gamma1-mypar.gamma2*cos(2*p1))+sin(2*p1)*mypar.gamma2*2*sin(2*p1) );
      atemp(2,1)=coe*(-mypar.gamma1*2*cos(2*p1)/mypar.gamma2+2*cos(4*p1) )*gp0;
      atemp(2,2)=0;
      atemp = atemp*mypar.T;
      dfdy(:,:,j)=atemp';
  end
elseif (iprob==8) % the active LCP, angle only.
  nn=length(t);
  eps=1e-6;
  %gp0=1/gp(0,mypar);
  coe=mypar.ac*mypar.gamma2*mypar.Er*0.25;
  for j=1:nn
      atemp=zeros(m,m);
      p1=y(j,1); %phi
      p2=y(j,2); %phi'
      %y3=y(j,3);
      %y4=y(j,4);
      atemp(1,1)=0; atemp(1,2)=1;
      atemp(2,1)=coe*(-mypar.gamma1*2*cos(2*p1)/mypar.gamma2+2*cos(4*p1) );
      atemp(2,2)=0;
      %atemp = atemp*mypar.T;
      dfdy(:,:,j)=atemp';
  end
elseif (iprob==9) % the active LCP, angle and 
  nn=length(t);
  eps=1e-6;
  gp0=1/gp(0,mypar);
  coe=mypar.ac*mypar.Er*0.25*gp0;
  for j=1:nn
      atemp=zeros(m,m);

      y1=y(j,1);
      y2=y(j,2);
      y3=y(j,3);
      y4=y(j,4);
      fsval=fsp(y3,mypar);
      atemp(1,1)=0; atemp(1,2)=1; atemp(1,3)=0; atemp(1,4)=0;
      atemp(2,1)=2*y2*y2*cos(2*y1)+(coe/y3/y3)*(2*cos(2*y1)*(mypar.gamma1-mypar.gamma2*cos(2*y1))+sin(2*y1)*mypar.gamma2*2*sin(2*y1) );
      atemp(2,2)=2*y2*sin(2*y1)-2*y4/y3;
      atemp(2,3)=2*y4*y2/y3/y3-(coe*2/y3/y3/y3)*sin(2*y1)*(mypar.gamma1-mypar.gamma2*cos(2*y1));
      atemp(2,4)=-2*y2/y3;
      atemp(3,1)=0;
      atemp(3,2)=0;
      atemp(3,3)=0;
      atemp(3,4)=1;
      atemp(4,1)=-mypar.beta3*mypar.Er*mypar.ac*gp0*sin(2*y1)*cos(2*y1);
      atemp(4,2)=2*y3*y2;
      atemp(4,3)=mypar.Er*mypar.nu*0.5*fsval+y2*y2;
      atemp(4,4)=0;
      dfdy(:,:,j)=atemp';
  end
else
    stop
end

return
