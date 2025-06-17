% parameters used in the model
function par=MALCPparams(isel)
%
    %global alpha1  alpha2  alpha3 alpha4 alpha5 alpha6

    %global Er activeC velV char_l char_l1 char_l2 beta3; 

    %global k1 k2 gamma1 gamma2;

    par.Re=1;     % Reynolds number 
    par.Er=1;     % Erickson number 
    par.ac=4/10/1.8365;%.388*10^4/par.Er;     % active number 
    %alpha = par.ac;
    alpha = 4;
    %par.ac=alpha;
    par.alpha1 = 0;
par.alpha2 = -1.4221*alpha;
par.alpha3 = -0.4144*alpha;
par.alpha4 = 4/3*10^-3;
par.alpha5 = 1 - par.alpha4 + 1.4221*alpha;
par.alpha6 = par.alpha5 - 1.8365*alpha;
par.gamma1 = par.alpha3 - par.alpha2
par.gamma2 = par.alpha6 - par.alpha5

%keep |gamma2|/4*A *Er=1
par.ac=-4/par.Er/par.gamma2; 
-par.gamma2 * par.ac * par.Er /4

par.gamma1/par.gamma2
%     par.alpha1= 0;    % viscosity
%     par.alpha2= -1.5;       % viscosity
%     par.alpha3= -0.5;       % viscosity
%     par.alpha4= 2;    % viscosity
%     par.alpha5= 2;    % viscosity
%     par.alpha6= 0;    % viscosity
%     par.gamma1= par.alpha3-par.alpha2;    % viscosity 
%     par.gamma2= par.alpha6-par.alpha5;   % viscosity 
    par.beta3=.45;
    par.nu=.5;
    %par.gamma1= 1.;                       % viscosity 
    %par.gamma2=-2.;                      % viscosity 
    par.phi0=0.35;
    par.L=.023191;
    %par.L = .005*pi;
    par.s0=1.0;
    par.T=1;

    par.threeD=false; % need to be false, 3D system (with v) is singular

%save(['March19Data_ActivityVSlanespace_noscale_phi1p0758.data'],'alpha','-ASCII', '-append')
%save(['March21Data_ActivityVSlanespace_noscale_phi1p0758.data'],'alpha','-ASCII', '-append')
