function f=newRHS(m,y,t)
% m is size of the problem
% t is the time to evaluate the jacobian, could be a scalar or a colomn vector
% y(1:m,1:length(t)) is y values
% f(k,i) is f(i) at time t(k)
global iprob mypar;

if (iprob==7) % alpha1 is not 0
  l=length(t);
  if (m~=2)
      stop
  end
  gp0=1/gp(0,mypar);
  
  f=zeros(l,m);
  eps=1.0e-6;
  coe=0.25*mypar.ac*mypar.Er*gp0;
  for i=1:l
    p1=y(i,1);
    p2=y(i,2);
    f(i,1)=p2;
    f(i,2)=mypar.gamma2*mypar.ac*mypar.Er*(-sin(2*p1)*mypar.gamma1/mypar.gamma2 + 0.5*sin(4*p1) );
    
  end
elseif (iprob==8) %the angle only ActiveLCP
  l=length(t);
  if (m~=2)
      stop
  end
  gp0=1/gp(0,mypar);
  
  f=zeros(l,m);
  eps=1.0e-6;
  coe=0.25*mypar.ac*mypar.Er*gp0;
  for i=1:l
    p1=y(i,1);
    p2=y(i,2);
    f(i,1)=p2;
    f(i,2)=mypar.gamma2*mypar.ac*mypar.Er*(-sin(2*p1)*mypar.gamma1/mypar.gamma2 + 0.5*sin(4*p1) );
    
  end
elseif (iprob==9) %the angle & s
  l=length(t);
  if (m~=4)
      stop
  end
  gp0=1/gp(0,mypar);
  
  f=zeros(l,m);
  eps=1.0e-6;
  coe = 0.25*mypar.ac*mypar.Er*gp0;
  
  for i=1:l
    y1=y(i,1);
    y2=y(i,2);
    y3=y(i,3);
    y4=y(i,4);
    fsval=fs(y3,mypar);
    f(i,1)=y2;
    f(i,2)=sin(2*y1)*(y2*y2+coe*(mypar.gamma1-mypar.gamma2*cos(2*y1)))-2*y2*y4/y3;
    f(i,3)=y4;
    f(i,4)=mypar.nu*mypar.Er*0.5*fsval+y3*y2*y2-0.25*mypar.beta3*mypar.Er*mypar.ac*gp0*sin(2*y1)*sin(2*y1);
  end
else
    stop
end

return
