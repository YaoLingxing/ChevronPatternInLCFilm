%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiallization. Note: You need to add your table path.
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long e; clc; %close all;
clf;

addpath 'GaussNode' % add path for node file

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  iprob: global variable for angle only or angle+order parameter
%
%- problem 8: Active LCP single phi ode
%- problem 9: Active LCP ODEs for phi and s
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iprob mypar;

mypar = MALCPparams(0); % parameter set for active LCP


gp0=gp(0,mypar);

omega0 = sqrt(mypar.ac*mypar.Er/gp0);
%mypar.T = 2*pi/omega0;

iprob=8;        % the index of the problem

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch iprob
  case 7 % ALCP alpha1 is not 0
      m=2;
      t0=0;
      tfinal=10;
      yinit(1) = 1.07;
      yinit(2) = 0;
      yanalfinal(  1) =  0;      
      yanalfinal(  2) =  0;          
  case 8 % ALCP phi ode system
      m=2;
      t0=0;
      tfinal=10.;
      yinit(1) = 1.07;
      yinit(2) = 0;
      yanalfinal(  1) =  0;      
      yanalfinal(  2) =  0;      
  case 9
	  m=4;
	  t0=0;
	  tfinal=4;
	  yinit(1)=pi/4; %angle
	  yinit(2)=.0; % angle'
	  yinit(3)=.5; %s 
	  yinit(4)=0; %s'
      yanalfinal(  1) =  0;      
      yanalfinal(  2) =  0;      
      yanalfinal(  3) =  0;      
      yanalfinal(  4) =  0;      
  otherwise
      stop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameter settings. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kmax=100;            %  maximum outer GMRES iterations
gtol=1e-2;           %  Tolerance for GMRES call(original 1e-11)

etol=1e-14;          %  Error tolerance for SDC step
nsteps=2500;
for nk = 1:1
  n=16; % fixed number of nodes used now
  for k =1:1         % not used

    %h0=tfinal*0.5^(k+1);  %  optional for time step size
    %h0=0.01;
    h0=(tfinal-t0)/nsteps;
    k0=max(n+4,81);

%%%%%%%%%%%%%%%%%%%%%%Main Subroutine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [ysol,ysolfin,res,indres,errrhs,inderr,iter]=...
        ALCPSDC(iprob,m,t0,tfinal,yinit,h0,n,kmax,gtol,etol,k0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end
end

%%%%%%%%%%%%%%%%%%%%The End***********************************************
figure(1)
%plot solution to angle
subplot(4,1,1)
t=[0:h0:tfinal];
plot(t,ysol(1,:),'.')
ylabel('$\phi(x)$','interpreter','latex')
xlabel(strcat('$x$, with',' $E_r=$',num2str(mypar.Er),', $\mathcal{A}$=',num2str(mypar.ac)),'interpreter','latex');

%%% data processing
subplot(4,1,2)

rlen=length(ysol(1,:));
istart=floor(nsteps/exp(1)+347);
phi=ysol(1,istart:rlen);
tt = t(1,istart:rlen);

gp0=gp(0,mypar);
w0=sqrt(0.25*mypar.ac*mypar.Er*(mypar.gamma1-mypar.gamma2));


[ppks,plocs]=findpeaks(phi);
nprd = mean(diff(tt(plocs)))

Perddifference = (diff(tt(plocs))-nprd)/nprd

findpeaks(ysol(1,:))
%tstr=strcat('FFT, period= ',num2str(Perd(1)), ' ');

title(strcat('use max value for lane= ',num2str(nprd)), 'interpreter','latex')
xlabel('index n')
ylabel('$\phi$', 'interpreter','latex')
% 

% plot velocity?
subplot(4,1,3)

for i =1:length(ysol(1,:))
    %vrs(1,i)=-sin(2*phi(i))/gp(phi(i),mypar);
    vrs(1,i)=-sin(2*ysol(1,i));
end
%vrs=-sin(2*phi)./gp(0,mypar);
velo=h0*cumsum(vrs);
plot(t,velo)
%xlabel('t')
ylabel('$w(x)$', 'interpreter','latex'); 
xlabel(strcat('$x$, $\alpha_1=$',num2str(mypar.alpha1), ....
    ', $\alpha_2=$',num2str(mypar.alpha2),',$\alpha_3=$',num2str(mypar.alpha3),...
    ', $\alpha_4=$',num2str(mypar.alpha4),',$\alpha_5=$',num2str(mypar.alpha5),...
    ', $\alpha_6=$',num2str(mypar.alpha6)),'interpreter','latex');
title('w(x)=$\int \sin(2\phi)/g(2\phi)dx$','interpreter','latex');

subplot(4,1,4)
ratio = mypar.gamma1/mypar.gamma2; 

xangle=[-pi:0.01:pi];
curve=-ratio*sin(2*xangle)+0.5*sin(4*xangle);
plot(xangle,curve)
xlabel('$x$', 'interpreter','latex')
ylabel('$-\gamma_1/\gamma_2 sin(2\phi)+0.5sin(4\phi)$', 'interpreter','latex');
title(strcat('$\gamma_1/\gamma_2$ ', num2str(ratio)), 'interpreter','latex');
