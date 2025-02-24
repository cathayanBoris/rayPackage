function runTRWpath2ways(lonInput,latInput,initialPeriodInDays,initialWavelengthInKM,numberOfDaysToRun, ...
    timestepInHr,forwardBackward,needBathy,smoothingParameterInKM)
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
% xq= X0deg*1000;    yq=Y0deg*1000;        %convert km to m
InitialPeriodInS =initialPeriodInDays*86400;  initialWavelengthInM =initialWavelengthInKM *1000;  %convert day to sec, km to m
smoothingParameterInM = smoothingParameterInKM * 1000;
% sgn_dt =-1; now is an input % +1 forward-time path & -1 backward-time path


global xg yg hg hxg hyg hxxg hxyg hyyg % declare global bathy grid vars
global meandx meandy originLon originLat omega
%CAUTION RE VARIABLE NAMES! last-letter "g" are grid variables length(xg) X length(yg),
%               for grid size(xg)= 1 X length(xg), size(yg)= length(yg) X 1
%  last-letter "q" are single-point variables 1 X 1
%  last-letter "p" are path variables, 1 X nsteps+1 along the ray-path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timestepInHr= forwardBackward*timestepInHr;   %prelim runs say < 3 hrs and >0.5 hrs is pretty good
%and 0.05 hr (3 minutes) runs SLOW, makes little improvement
%number of steps one-way transit
nsteps=floor(abs(numberOfDaysToRun*24/timestepInHr))+1; % the first entry is reserved for t=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ifcalcbathy=0; %is an input now %set =1 first time as global variable; thereafter set =0
% if(ifcalcbathy) % then must sppecify simulated Bath parameters
%     %   Bath=[Rslope,A,dmax,Rflat,Redge,dx]; typical #'s listed as follows
%  Bath= [34E3,3000,3800,60E3,2E3,0.50E3]; %units (m) % 0.5 km bathy grid,integer multiples of dx
%     %  Bath= [33E3,3000,3800,60E3,2E3,1E3];  % this 1km grid is what I suspected
%     %                may be too coarse??
%
% %%%%%%%%%%%%%%%%% whole gridded fields of bathymetry & gradients %%%%%%%%%
%     ifplots=1;  % first pass, calc bathy and gradients
%     [Bgxy, Bgrid] =bathy_sim(Bath,ifplots);
%     lx=length(Bgxy); xg=Bgxy(1:lx/2); yg=Bgxy(lx/2+1:lx)';
%     [nr,~]=size(Bgrid);
%     hg =Bgrid(:,1:nr); hxg=Bgrid(:,nr+1:2*nr); hyg=Bgrid(:,2*nr+1:3*nr);
%     hxxg=Bgrid(:,3*nr+1:4*nr); hxyg=Bgrid(:,4*nr+1:5*nr); hyyg=Bgrid(:,5*nr+1:6*nr);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % this grid of bathymetric variables is used in Env_parms below
% else
%     % these are saved GLOBAL VARIABLES after calculating in 1st pass
%     Bgxy= [xg, yg'];  % define here if you don't recalculate bathy_sim
%     Bgrid= [hg, hxg, hyg, hxxg, hxyg, hyyg];
% end %ifcalcbathy  %otherwise use the already-calculated global bathy vars
if needBathy == 1
    [xtopo,ytopo,ztopoS] = getMeSmoothed(smoothingParameterInM,-92,-88,latInput-2,latInput+2);
    %
    % [ztopo] = getMePlateTopo(2800,0.005,xtopo,ytopo,ztopo);

    hg = -ztopoS;
    [hxg,hyg,~,dx,dy] = getMeGradZ(xtopo,ytopo,hg);
    [hxxg,hxyg,~,~,~] = getMeGradZ(xtopo,ytopo,hxg);
    [~,hyyg,~,~,~] = getMeGradZ(xtopo,ytopo,hyg);
    originLon = nanmean(xtopo);
    originLat = nanmean(ytopo);
    meandx = nanmean(dx); % NOTE: dx is not uniform
    meandy = nanmean(dy);
    xg = (xtopo - originLon) * meandx * 60;
    yg = (ytopo - originLat) * meandy * 60;

end

% starting point lonlat -> xy
xq = (lonInput - originLon) * meandx * 60;
yq = (latInput - originLat) * meandy * 60;

% Live Map
figure(100) 
CI = 100;
v=0:CI:4000;
contour(xg,yg,hg',v,"ShowText","off"),
clim([0 4000])
%axis square
hold on
[C,lab]=contour(xg,yg,hg',[3000:-500:1500],'k');  clabel(C,lab)
title(['TRW path ',num2str(initialPeriodInDays),'d',num2str(initialWavelengthInKM),'km'])
xlabel('x (m)'); ylabel('y (m)');
text(0, -2e5, ['depth CI ', num2str(CI,4),'(m)'])

% initially need starting depth h=hq at initial point (xq,yq) to calc NB(h)
%[hq,hxq,hyq,hxxq,hxyq,hyyq] =bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg)
% [hq,~,~,~,~,~] =bathy_point(xq,yq,xg,yg,hg,hxg,hyg,hxxg,hxyg,hyyg);
% NB=NBh(hq); % NOW IS NOT USED

%%%%%%%%%%%%%%%% plot Tmin = 1/( N*alpha ) ************
% figure(10)   % simplified case T=1/N*alpha
% % Tmin= 2*pi/sigmax/86400; with sigmax=N*alpha
% con=2*pi/86400; % convert period from seconds to days
% CId=20;        % (days)
% v=20:CId:200 ;
% contour(xg,yg,(con./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)))',v),
% %axis square
% hold on
% %     [C,lab]=contour(xg,-yg,con./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)),[2,5,15,45],'k');  clabel(C,lab)
% [C,lab]=contour(xg,yg,(con./(NBh(hg).*sqrt(hxg.*hxg+hyg.*hyg)))',[2,7,20,60],'k');  clabel(C,lab)
% title('Tmin=1/N*alpha (d) for tanh=1, beta=0')
% xlabel('x (m)'); ylabel('y (m)');
% %     text(63e3, -87e3, ['CId ', num2str(CId,3),'(days)'])
% text(-20e4, 15e4, 'contour labels (days)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lat=27.2;  % latitude (decimal degrees), Peter Hamilton used 27.2
omega=7.292e-5;   %% Earth's rotation rate(s-1)   A.E.Gill p.597
f = 2*omega*sind(latInput);        % exactly same formula as f=sw_f(lat);
% beta=2.0E-11;
beta = 2*omega*(sind(latInput+0.01) - sind(latInput-0.01))/(0.02*pi*6371000/180);
% beta=0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% HERE'S A WAY TO COMPARE THE REVERSE PATH AND FORWARD PATH %%%%%
% %%%% GOOD IDEA TO FIRST RUN WITH -VE dt %%%%%
% npathrun=1;  % just do both ways in one pass; eliminate the old for-loop
%     Startvars =[xq, yq, T0day,LAM0km, ndaysrun,dt_hr,sgn_dt, npathrun];
% % %%%%%  THEN KEEP (sig, k,l, Ei) without reinitializing AND
% %  JUST REVERSE +dt DIRECTION INSIDE TRWpath_sig1_2way %%%%%
% % % use Endvars as next Startvars, but reverse sign of dt
% % % Endvars is [xq, yq, T0day,LAM0km, ndaysrun,dt_hr,sgn_dt, npathrun];
%     if(npathrun>1)
%         xq=Endvars(1); yq=Endvars(2); T0day =Endvars(3); LAM0km=Endvars(4);
%         ndaysrun=Endvars(5); dt_hr=Endvars(6);
%         sgn_dt=-Endvars(7);   %I found a simpler way to reverse the path
%         % SO THESE ARE JUST FOR RECORD-KEEPING; save reverse track startup
%         Startvars =[xq, yq, T0day, LAM0km,ndaysrun,dt_hr,sgn_dt,npathrun];
%     end %if
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODIFIED TRWpath_sig1_2way to use existing (k,l,sig) ON RETURN PATH
%        and NOT re-initialize

% 9 vars of output packs into Ei: Environmental parameters
Ei = Env_parms(xq,yq,xg,yg,hg',hxg',hyg',hxxg',hxyg',hyyg',f,beta);


lxg = length(xg);
lyg = length(yg);
% turn xg, yg into xg-yg matrix
xg = repmat(xg,1,lyg); % dimension lxg x lyg
yg = repmat(yg,1,lxg); % dimension lyg x lxg

Bgxy= [xg, yg'];  % dimension: lxg x 2lyg
Bgrid= [hg', hxg', hyg', hxxg', hxyg', hyyg']; % dimension: lyg x 6lxg

xg = xg(:,1);
yg = yg(:,1);

% PVAR=TRWpath_sig1_2way(T0 ,LAM0 ,xq,yq, Ei, dt_hr,nsteps,Bath); %TRWpath does the hard work
PVAR =TRWpath_sig1_2ways(InitialPeriodInS , initialWavelengthInM , xq,yq, Ei, timestepInHr,nsteps, Bgxy, Bgrid, smoothingParameterInM);
if PVAR == 999
    'The execution was a failure.'
    return
end
% output PVAR=[xp,yp,kp,lp,K0p,LAMkmp,sigp,Tdayp,cgxp,cgyp,hp,hxp,hyp]; %along path

%  [m,~] = size(PVAR);
%  Endvars is [xq, yq, T0day,LAM0km, ndaysrun,dt_hr,sgn_dt, npathrun];
% Endvars=[PVAR(m,1),PVAR(m,2),PVAR(m,8),PVAR(m,6),ndaysrun,dt_hr,sgn_dt,npathrun];
