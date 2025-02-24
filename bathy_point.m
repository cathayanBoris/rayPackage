function [hq,hxq,hyq,hxxq,hxyq,hyyq] =bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg)
% function
% [hq,hxq,hyq,hxxq,hxyq,hyyq]=bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg);
% use % Vq = interp2(X,Y,V,Xq,Yq) interpolates to find Vq at one loc (xq,yq)
%  last letter "q" is for a single point
%  last letter "g" is whole grid
% inputs: 
% Rslope=25E3; A=2500; dmax=4000; ifplots=1;
% Bath=[Rslope,A,dmax,Rflat,Redge,dx]; typically:
% %                      Bath= [25E3,2500,4000,40E3,15E3,1E3];

%%%%%%%%%%%%%%%%%%% output for whole x,y,grid   %%%%%%%%%%%%%%%%%%%%
%[xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg]=bathy_sim(Bath,ifplots); %whole x,y grid
%%%%%%%%%%%%%%%%%%%%

hq=interp2(xg,yg,hg,xq,yq,'spline');  % agree exactly at grid points

hxq=interp2(xg,yg,hxg,xq,yq,'spline');
hyq=interp2(xg,yg,hyg,xq,yq,'spline');

hxxq=interp2(xg,yg,hxxg,xq,yq,'spline');
hxyq=interp2(xg,yg,hxyg,xq,yq,'spline');
hyyq=interp2(xg,yg,hyyg,xq,yq,'spline');