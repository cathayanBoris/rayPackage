function [Bgxy, Bgrid] =bathy_sim(Bath, ifplots)
% function [xg,yg,hg, hxg, hyg, hxxg, hxyg, hyyg] =bathy_sim(Bath,ifplots);
% units all SI for Rslope, A, dmax, (x,y),h, ...
%
% xg and yg vectors define the grid 
% all "h...g" variable outputs are matrices that define the whole field 
% for mapping h(x,y) and its derivatives, which all
% get calculated in this m-file by applying gradient.m to the gridded h 

% % clear, close, Bath= [25E3,2500,4000,40E3,15E3,1E3]; ifplots=1;
% Bath=[Rslope,A,dmax,Rflat,Redge,dx]; typical #'s above

% dx=grid spacing for a contour plot (x,y,h)
% This bathymetry will be circularly symmetric, composed of
% (1) a flat region in the middle radius Rflat, depth dmax
% (2) then a paraboloid ring, rising -A*r^2, radius r=(0, Rslope)
% (3) lastly a flat outer ring from radius (Rflat+Rslope) out to Redge. 
% NOTE all three R's are integer multiples of dx
% Rslope is an input variable, typically Rslope=25E3;
Rslope= Bath(1);
A=      Bath(2);
dmax=   Bath(3);
% Rflat= 40E3;
% Redge= 15E3;
Rflat= Bath(4);
Redge= Bath(5);
dx=    Bath(6);  dy=dx;
Rtot =Rflat+Rslope+Redge;

l=Rtot/dx; % integer # grid points along full radius
L=2*l+1;                   % integer # grid points along diameter

xg=dx*(-l:l);    
yg=dy*(-l:l)';      %transposed, dx and dy have SAME magnitude value
% axis([XMIN XMAX YMIN YMAX]) %XMIN=-Rtot; Xmax=Rtot; YMIN=-Rtot;YMAX=Rtot;

if( L~=length(xg) && L~=length(yg))   %write a warning if you get here
    'something in bathy_sim is screwed up, stop'
    return
end 

%%%%%
hg=ones(L, L);   %preallocate all outputs
hxg=hg; hyg=hg; hxxg=hg; hxyg=hg; hyyg=hg;   %preallocate all outputs
%%%%% 
r2=ones(L,L);
r2=xg.*xg + yg.*yg;  %size(r2) % radius squared to point (x,y)

% keyboard 

% the parabolic ring  h=dmax - A*(r-Rflat)^2
   
Rflat2=Rflat*Rflat;
hg(find(r2<Rflat2))=dmax;  % deep depth of central flat region

Rinner_edge2=(Rflat+Rslope)*(Rflat+Rslope);
dmin=dmax-A;
hg(find(r2>=Rinner_edge2))= dmin;  % shallow depth of outer ring beyond slope

hg(find((r2>=Rflat2)&(r2<Rinner_edge2)))...
    =dmax-(A/Rslope^2)*(sqrt(r2(find((r2>=Rflat2)&(r2<Rinner_edge2))))-Rflat).^2; % paraboloidal slope

[hxg,hyg] =gradient(hg); hxg= hxg/dx; hyg=hyg/dy;  % equal grid spacing dy=dx

[hxxg,hxyg] =gradient(hxg); hxxg=hxxg/dx; hxyg=hxyg/dy;
[hyxg,hyyg] =gradient(hyg); hyxg=hyxg/dx; hyyg=hyyg/dy;  %unused hyx is=hxy
% now we have computed all the output variables

Bgxy= [xg, yg'];
Bgrid= [hg, hxg, hyg, hxxg, hxyg, hyyg];  % =bathy_sim(Bath,ifplots);

% keyboard

if(ifplots==1)

%     figure(1) %always want this plot - bring outside this funcion
%     CI=250; 
%     v=[500:CI:4000];
%     contour(xg,-yg,hg,v),
%      axis square
%      hold on
%     [C,lab]=contour(xg,-yg,hg,[3000,1500],'k');  clabel(C,lab)
%     title('paraboloidal ring bathymetry h(x,y) (m)')
%     xlabel('x (m)'); ylabel('y (m)');
%     text(65e3, -85e3, ['CI ', num2str(CI,4),'(m)'])

    figure(2)
    % [hx,hy] =gradient(h); hx= hx/dx; hy=hy/dx;  % equal grid spacing dy=dx
    alpha = sqrt(hxg.*hxg + hyg.*hyg);           % size(alpha) %just testing - okay 
    CI=0.04; v=-0.40:CI:0.40;
    subplot(2,2,3), contour(xg,-yg,hxg,v)
       axis square,  title('hx (dimensionless)')
    subplot(2,2,2), contour(xg,-yg,hyg,v)
       axis square,  title('hy (dimensionless)')
    subplot(2,2,1), [C,lab]=contour(xg,-yg, alpha,v);  
       axis square, 
       clabel(C,lab,[0, 0.12])
       title('alpha=sqrt(hx^2+hy^2)')
       text(45e3, -83e3, ['CI= ', num2str(CI,3)])

  % [max(max(hxx)),min(min(hxx)),max(max(hxy)),min(min(hxy)),max(max(hyy)),min(min(hyy))];
  % has isolated extreme values at edges where 2nd derivatives are not
  % continuous;  don't worry about it - path should not approach these edges

    % % just checked that gradient operator is symmetric
    % diffhyx = hyx-hxy; % expect = zero; numerical vals <~+/- 1e-19; 
    % [ max(max(diffhyx)), min(min(diffhyx)) ]  % YES!
    %    figure(4),  contour(x,-y,diffhyx), axis square  %pretty pattern!

    v=1.0E-5*[-1.0:0.05:+1.0];
    figure(3) 
    subplot(3,1,1), contour(xg,-yg,hxxg,v)
       axis square, title('hxx (1/m)')
    subplot(3,1,2), contour(xg,-yg,hxyg,v)
       axis square, title('hxy (1/m)')
    subplot(3,1,3), contour(xg,-yg, hyyg,v)
       axis square, title('hyy (1/m)')

end %ifplots

return


% expecting that Rslope =25E3 usually, and A=2500; 
% Rflat= 40E3;
% Redge= 15E3;

