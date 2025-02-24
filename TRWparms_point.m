function [cgx,cgy,dkdt,dldt] =TRWparms_point(k,l,sig,Ei) 
% function [cgx,cgy,dkdt,dldt] =TRWparms_point(k,l,sig,Ei) 
% at loc=(x,y) for wave defined by (sig,k,l), calculates parms for TRWpath
%   called by TRWpath, which uses parms to do Runge-Kutta integration
%       assumes (k,l sig, Ei) are all input from TRWpath, 
%   and the input loc is (x,y)   
%
%   input variable is Ei, used to calc [cgx,cgy,dkdt,dldt] 
%
%   includes beta, uses constant NB

% environment parms Ei= [h hx hy hxx hxy hyy NB f beta]; %are INPUT
% USAGE WILL BE
h=     Ei(1); hx= Ei(2);   hy= Ei(3);
hxx=   Ei(4); hxy=Ei(5);  hyy= Ei(6); 
NB=    Ei(7); f=  Ei(8); beta= Ei(9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kbsqr =(k*k + l*l + k*beta/sig);    % K with beta     - uses input sig

if(Kbsqr <= 0)
    'in TRWparms Kb^2 is < 0'
    cgx=999;
    cgy=999;
    dkdt=999;
    dldt=999;
    return % abort simulation
end %if
Kb = sqrt(Kbsqr);
kb2s=k+ beta/(2*sig) ;         % dKb/dk term, with beta - uses input sig 
KXgh = (hy*k - hx*l);             % "K cross grad h"
mu= NB*Kb/f;                     % mu has Kb
zeta= mu*h;                      % argument in tanh and sech^2
Tz = tanh(zeta);                  % the tanh term
TTT=(1-Tz*Tz)/Tz;

%%%%%%%%%%%%%%%get the answers  %%%%%%%%%%%%%%%%
cgx =sign(f)*sig *(( hy/KXgh) -((kb2s)/Kb/Kb)*(1 +zeta*TTT));
cgy =sign(f)*sig *((-hx/KXgh) -(    l/Kb/Kb) *(1 +zeta*TTT));

dkdt = -sign(f)*sig *(((k*hxy -l*hxx)/KXgh) -mu*TTT*hx); 
dldt = -sign(f)*sig *(((k*hyy -l*hxy)/KXgh) -mu*TTT*hy); 

%  keyboard  %debugging 

 return

 % for debugging only
% calculate fraction due to tanh term zeta*TTT - only affects (cgx, cgy)
% zeta*TTT; , compared to (tiny?) beta term (1)
% calculate ratios of (variable tanh part)/(const tanh part)

% ratiodkdt = -mub*TTT*hx/((k*hxy -l*hxx)/KXgh);
% ratiodldt = -mub*TTT*hy/((k*hyy -l*hxy)/KXgh);
% ratiocgx = -((kb2s)/Kb/Kb)*(1 +zeta*TTT)/( hy/KXgh);
% ratiocgy = -(  l/Kb/Kb  ) *(1 +zeta*TTT)/(-hx/KXgh);
% [ratiodkdt ratiodldt ratiocgx ratiocgy zeta*TTT ]
% 


