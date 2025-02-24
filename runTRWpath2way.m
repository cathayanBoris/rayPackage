function runTRWpath2way(X0km,Y0km,T0day,LAM0km,ndaysrun,dt_hr,sgn_dt,...
    ifcalcbathy,Bath) 
% function runTRWpath2way(X0km,Y0km,T0day,LAM0km,ndaysrun,dt_hr,sgn_dt,...
%    ifcalcbathy,Bath)
% automatically runs and reverses the path (2 way)
% these input variables are the most common ones to enter in different runs
%     
%  X0km,Y0km,T0day,LAM0km are launch position, and TRW period and wavelength 
%  example starting locations at mid-quartile for DRW's simulated bathy
%       xq=+62E3;  yq= +62E3; % Q1   
%       xq=-64.5E3;  yq= +64.5E3; % Q2, near h=1000m shallowest edge
%       xq= -64.5E3*sqrt(2); yq=0; %is at equivalent depth along western bdy
%       xq=-62E3;  yq= -62E3; % Q3
%       xq=+62E3;  yq= -62E3; % Q4 
%
%  ndaysrun,dt_hr,sgn_dt are time-interval, abs val of time-step, 
%      and initial direction, e.g., sgn_dt=-1 runs reverse path first
%                                or sgn_dt=+1 runs forward path first 
%
% hard-wired environmental variables are h-dependent NB(h), f(lat), beta
%
% bathymetry can be calculated from Bath when ifcalcbathy=1, and saved as 
% GLOBAL VARIABLES:  Bgxy= [xg, yg'];  
%                    Bgrid= [hg, hxg, hyg, hxxg, hxyg, hyyg];
% subsequently set ifcalcbathy=0, and the GLOBAL VARIABLES will be used 

% clear 
close all
xq= X0km*1000;    yq=Y0km*1000;        %convert km to m
T0 =T0day*86400;  LAM0 =LAM0km *1000;  %convert day to sec, km to m

% sgn_dt =-1; now is an input % +1 forward-time path & -1 backward-time path


global xg yg hg hxg hyg hxxg hxyg hyyg % declare global bathy grid vars

%CAUTION RE VARIABLE NAMES! last-letter "g" are grid variables length(xg) X length(yg),
%               for grid size(xg)= 1 X length(xg), size(yg)= length(yg) X 1
%  last-letter "q" are single-point variables 1 X 1
%  last-letter "p" are path variables, 1 X nsteps+1 along the ray-path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% you may rarely want to change these path-length variables %%%%
ndays=ndaysrun;  % input var =number of days to run the set of TRW paths
dt_hr= sgn_dt*dt_hr;   %prelim runs say < 3 hrs and >0.5 hrs is pretty good
             %and 0.05 hr (3 minutes) runs SLOW, makes little improvement
     %number of steps one-way transit
nsteps=floor(abs(ndays*24/dt_hr)); % integer abs to also work if dt_hr is -ve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ifcalcbathy=0; %is an input now %set =1 first time as global variable; thereafter set =0
if(ifcalcbathy) % then must sppecify simulated Bath parameters
    %   Bath=[Rslope,A,dmax,Rflat,Redge,dx]; typical #'s listed as follows
 Bath= [34E3,3000,3800,60E3,2E3,0.50E3]; %units (m) % 0.5 km bathy grid,integer multiples of dx
    %  Bath= [33E3,3000,3800,60E3,2E3,1E3];  % this 1km grid is what I suspected 
    %                may be too coarse??

%%%%%%%%%%%%%%%%% whole gridded fields of bathymetry & gradients %%%%%%%%%
    ifplots=1;  % first pass, calc bathy and gradients
    [Bgxy, Bgrid] =bathy_sim(Bath,ifplots);
    lx=length(Bgxy); xg=Bgxy(1:lx/2); yg=Bgxy(lx/2+1:lx)';
    [nr,~]=size(Bgrid);
    hg =Bgrid(:,1:nr); hxg=Bgrid(:,nr+1:2*nr); hyg=Bgrid(:,2*nr+1:3*nr);
    hxxg=Bgrid(:,3*nr+1:4*nr); hxyg=Bgrid(:,4*nr+1:5*nr); hyyg=Bgrid(:,5*nr+1:6*nr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this grid of bathymetric variables is used in Env_parms below
else
    % these are saved GLOBAL VARIABLES after calculating in 1st pass
    Bgxy= [xg, yg'];  % define here if you don't recalculate bathy_sim
    Bgrid= [hg, hxg, hyg, hxxg, hxyg, hyyg];
end %ifcalcbathy  %otherwise use the already-calculated global bathy vars

    figure(1) %always want this plot 
    CI=250; 
    v=500:CI:4000;
    contour(xg,-yg,hg,v),
     axis square
     hold on
    [C,lab]=contour(xg,-yg,hg,[3000,1500],'k');  clabel(C,lab)
    title('TRW path on ring of paraboloidal bathymetry h(x,y) (m)')
    xlabel('x (m)'); ylabel('y (m)');
    text(50e3, -90e3, ['depth CI ', num2str(CI,4),'(m)'])


% initially need starting depth h=hq at initial point (xq,yq) to calc NB(h)
%[hq,hxq,hyq,hxxq,hxyq,hyyq] =bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg)
% [hq,~,~,~,~,~] =bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg); 
% NB=NBh(hq); % NOW IS NOT USED

%%%%%%%%%%%%%%%% plot Tmin = 1/( N*alpha ) ************
   figure(10)   % simplified case T=1/N*alpha
   % Tmin= 2*pi/sigmax/86400; with sigmax=N*alpha 
   con=2*pi/86400; % convert period from seconds to days
    CId=20;        % (days)
    v=20:CId:200 ;
    contour(xg,-yg,con./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)),v), 
       axis square
       hold on
%     [C,lab]=contour(xg,-yg,con./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)),[2,5,15,45],'k');  clabel(C,lab)
     [C,lab]=contour(xg,-yg,con./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)),[2,7,20,60],'k');  clabel(C,lab)
    title('Tmin=1/N*alpha (d) for tanh=1, beta=0')
    xlabel('x (m)'); ylabel('y (m)');
%     text(63e3, -87e3, ['CId ', num2str(CId,3),'(days)'])
    text(40e3, -91e3, ['contour labels (days)'])

    % fig 11 is a more accurate map of Tmin, accounting for the tanh term
%    figure(11) % 
%   % Tmin= 2*pi/sigmax/86400; with sigmax=N*alpha/tanh(zeta);
%   % zeta=(K*h*N/f); small K "illustrate" = Kill=2*pi/160e3;
%   Lill=160e3;     % choose longish wavelength to illustrate
%   Kill=2*pi/Lill; % illustrate for a longish wave, smallish K, smallish Tmin
%   con=2*pi/86400;
%   lat=27.2;  % latitude (decimal degrees)
%   omega=7.292e-5;   %% Earth's rotation rate(s-1)   A.E.Gill p.597
%   f= 2*omega*sind(lat); % same as f=sw_f(lat);
%     CId=20;        % (days)
%     v=20:CId:200 ;
%        Tanhillg= tanh(Kill.*hg.*NBh(hg)/f);
%     contour(xg,-yg,con*Tanhillg./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)),v), 
%        axis square
%        hold on
%     [C,lab]=contour(xg,-yg,con*Tanhillg./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)),[2,15,40],'k');  clabel(C,lab);
%         title('Tmin=tanh/(N*alpha) (d),  tanh<1, const LAM')
%         xlabel('x (m)'); ylabel('y (m)');
%         text(50e3, -88e3, ['Lill= ', num2str(Lill/1e3,3), '(km)'] )

    % return  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ifckkl=0;      % don't confuse ifckkl with ifbathy

lat=27.2;  % latitude (decimal degrees), Peter Hamilton used 27.2
omega=7.292e-5;   %% Earth's rotation rate(s-1)   A.E.Gill p.597
f= 2*omega*sind(lat);        % exactly same formula as f=sw_f(lat);
beta=2.0E-11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HERE'S A WAY TO COMPARE THE REVERSE PATH AND FORWARD PATH %%%%%
%%%% GOOD IDEA TO FIRST RUN WITH -VE dt %%%%%
npathrun=1;  % just do both ways in one pass; eliminate the old for-loop
    Startvars =[xq, yq, T0day,LAM0km, ndaysrun,dt_hr,sgn_dt, npathrun];
% %%%%%  THEN KEEP (sig, k,l, Ei) without reinitializing AND 
%  JUST REVERSE +dt DIRECTION INSIDE TRWpath_sig1_2way %%%%%
% % use Endvars as next Startvars, but reverse sign of dt
% % Endvars is [xq, yq, T0day,LAM0km, ndaysrun,dt_hr,sgn_dt, npathrun];
    if(npathrun>1)
        xq=Endvars(1); yq=Endvars(2); T0day =Endvars(3); LAM0km=Endvars(4);
        ndaysrun=Endvars(5); dt_hr=Endvars(6);
        sgn_dt=-Endvars(7);   %I found a simpler way to reverse the path
        % SO THESE ARE JUST FOR RECORD-KEEPING; save reverse track startup
        Startvars =[xq, yq, T0day, LAM0km,ndaysrun,dt_hr,sgn_dt,npathrun]; 
    end %if
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODIFIED TRWpath_sig1_2way to use existing (k,l,sig) ON RETURN PATH 
%        and NOT re-initialize

Ei = Env_parms(xq,yq, xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg,f,beta);
% 9 vars of output packs into Ei

% PVAR=TRWpath_sig1_2way(T0 ,LAM0 ,xq,yq, Ei, dt_hr,nsteps,Bath); %TRWpath does the hard work
PVAR =TRWpath_sig1_2way(T0 , LAM0 , xq,yq, Ei, dt_hr,nsteps, Bgxy, Bgrid);

% output PVAR=[xp,yp,kp,lp,K0p,LAMkmp,sigp,Tdayp,cgxp,cgyp,hp,hxp,hyp]; %along path

 [m,~] = size(PVAR);
 % Endvars is [xq, yq, T0day,LAM0km, ndaysrun,dt_hr,sgn_dt, npathrun];
Endvars=[PVAR(m,1),PVAR(m,2),PVAR(m,8),PVAR(m,6),ndaysrun,dt_hr,sgn_dt,npathrun]; 
