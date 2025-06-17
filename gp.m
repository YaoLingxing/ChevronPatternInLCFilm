function out=gp(phi,par)
alp1 = par.alpha1;
alp2 = par.alpha2;
alp3 = par.alpha3;
alp4 = par.alpha4; 
alp5 = par.alpha5;
alp6 = par.alpha6;


% new version;
out = 0.5*alp1*(sin(2*phi).*sin(2*phi)+sin(2*pi).*(1-cos(2*phi))) + ...
    0.5*(alp2+alp5).*(1-cos(2*phi))+0.5*(alp6-alp3).*(1+cos(2*phi))+alp4;
