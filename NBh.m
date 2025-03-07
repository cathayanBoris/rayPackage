function NB = NBh(h)
% function NB = NBh(h)
% make buoyancy freq NB depend in exponental downward decay with depth h
%       NB =NN*exp(-h/HH)
% the two constants NN and HH are fitted by two value-depth pairs 
% hi,Ni,TB_hr =(1500, 8e-4, 2.2) &(3000, 2.4e-4, 7.3); Aryan and Boris profiles
% test point from (hi,Ni,TB_hr) = (2500, 2.7e-4, 6.5);
%   HH= -(h2-h1)/ln(N2/N1);  %  NN= N1*exp(h1/HH); % fits the two constants  
%  NB= NN*exp(-h/HH);

h = h-500;

h1=1500; N1=8e-4;
h2=3000; N2=2.4e-4;
HH= -(h2-h1)/log(N2/N1);
NN= N1*exp(h1/HH);

% NB= NN*exp(-[1500:500:3000]/HH)     %test values look good
NB =NN*exp(-h/HH); 
end

% NB = 1.4E-3;
% NB = 2.65E-4;
% test values:
%   htest=500:500:4000;  NBh(htest)*1000   % units 1000*(1/s)
% 1.7852    1.1950    0.8000    0.5355    0.3585    0.2400    0.1607    0.1076

