function Ei = Env_parms(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg,f,beta)
% function Ei = Env_parms(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg,NB,f,beta);
%
%    calculates Ei(1:9) at position (x,y)
%        position units for (x,y) are (m);
% Ei =[h,hx,hy,hxx,hxy,hyy,NB,f,beta];
% Ei(1)=h;   Ei(2)=hx;  Ei(3)=hy;
% Ei(4)=hxx; Ei(5)=hxy; Ei(6)=hyy;
% Ei(7)=NB;   Ei(8)=f;   Ei(9)=beta;
% USAGE WILL BE
% h=    Ei(1); hx= Ei(2); hy= Ei(3);
% hxx=  Ei(4); hxy=Ei(5); hyy=Ei(6);
% NB=    Ei(7); f=  Ei(8); beta= Ei(9);

% 2D interpolated local environmental parameters
[hq,hxq,hyq,hxxq,hxyq,hyyq]=bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg);


Ei(1)=hq;      Ei(2)=hxq;    Ei(3)=hyq;
Ei(4)=hxxq;    Ei(5)=hxyq;   Ei(6)=hyyq;
Ei(7)=NBh(hq);                  % function NBh(h) finds new NB at depth hq
Ei(8)=f;      Ei(9)=beta;
% % Ei(7)=Nfb(1);  Ei(8)=Nfb(2); Ei(9)=Nfb(3);
%Ei(7) =  1.4E-3; % want constant N
return


% %  Tbuoy= 1.0; %(buoyancy period,hours) %T=1 hr would be NB=0.001745
% NB=2*pi/(Tbuoy*3600);  % units (1/s) % could later make depth-dependent N(h)
% f=sw_f(lat);
% beta=2.0E-11;
% NB=Nfb(1); f=Nfb(2); beta=Nfb(3);
% Tbuoy=Nfb(4);  lat=Nfb(5);




