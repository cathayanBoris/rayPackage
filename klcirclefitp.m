function [k,l, sigfit] = klcirclefitp(sig, K0, Ei, iUPDN)
%  function [k,l,sigfit] = klcirclefitp(sig, K0, Ei, iUPDN)
%
% "p" appended for path, requiring fine-comb initial (k,l); NO ifcheck 
%
% calculate TRW wavevector components (k,l) 
% from input freq (sig, 1/s) and wavenumber magnitude (K0, 1/m)
% fit in (k,l) plane along 1D angle on K0 circle
% can optimize for best direct sig fit and tanh fit

% for ifcheck=1, plots rhs vs lhs crossings after best k,l chosen

%    input Ei= Environ_parms
h =  Ei(1);    hx=Ei(2);   hy=Ei(3);  
% hxx =Ei(4);  hxy=Ei(5);  hyy=Ei(6); 
beta=Ei(9);    f=Ei(8);    NB=Ei(7);  

% if( ifcheck~=1) % DRW deleted ifcheck input variable
%     ifcheck=0; %default no printed check that (k,l) gives the correct sig
% end

%%%%%%%%%
% B   = [hx, hy];  % bathymetric grad(h) downslope vector
Bang=atan2(hy,hx);  % math angle of grad(h) vector
% search in (k,l) plane along K0 circle, as a function of angle
% Nsearch=400000;   % e.g., 1000 would be approx 0.1 degree steps; 
%     the T0day and sig1 is offset ~100ppm 
% Nsearch=16E5;   % sig1 has 10ppm diff
% Nsearch-32E5; % used in dispersion curves 
%        and less than 1m different end-location after 2.5d at 1 hr steps
Nsearch=264E5; % runs fast enough and matches end-path of 128E5
% Nsearch=128E5; %  end-point agrees to 1m
% Nsearch=160E5;   % try 8X less fine; runs fast enough on RW's MacBook; 
% and path end-point again agrees to 1m  
% Nsearch=128E6;  % OVERKILL for initial (k,l) in TRWpath
ang=(pi/(2*Nsearch))*(0:(Nsearch)); %divide a quadrant into Nsearch angles

% iUPDN = +ve or -ve determines usual UPSLOPE, DNSLOPE phase prop
%             and respectively DN and UP group vel propagation cg
angadd=pi/2;
if(iUPDN>0)   % UPSLOPE phase propagation
    bang=(Bang-angadd)-ang;   %angle to right of vector BXz
elseif(iUPDN<0)         % normal case (-1) is DNSLOPE phase propagation
    %angadd=0;
    bang=Bang-ang;      % angle in quadrant to right of angle B=grad(h)
else   % iUPDN is zero, so do both
%     bang=Bang-ang-ang; % 2*ang spans half-plane of "topographic westward" 
    bang= [ (Bang-ang), ((Bang-angadd)-ang)] ;
end 
% bang is the search-set of angles in the proper quadrants

% calculate cos and sin just Nsearch times for one input K0
% % already know Lbang=Nsearch+1 ;   
Lbang= length(bang); % Lbang is =Nsearch or for IUPDN=0, is =2*Nsearch
cb=cos(bang);  
sb=sin(bang); 
kt=K0*cb;   % test set of (k,l) vectors along K0 circle
lt=K0*sb;
%%%%%%%%%%%%%%%%%%%%%%%%%

% % re-express sig = NB*KXgh/(Kb*tanh(zeta)); 
% % in terms of k,l   which are on K0 circle

lhs= ones(1, Lbang);   %initialize sizes
rhs= ones(1, Lbang); 

%%%%% sig minimization loop to prep LHS and RHS for optimization
for ja = 1:Lbang  %fine comb search for best kt,lt ang fit

   kta = kt(ja);  % wavevector x-component, includes +/- vals if needed
   lta=  lt(ja);  % wavevector l component, includes +/- vals if needed
   % work with one value at a time
   Kbsq =K0*K0 +beta*kta/sig; % K0 & sig= input vars
   if Kbsq<=0          % skip imaginary calculations 
       rhs(ja)= NaN;  % want to ignore, would be imaginary number
   else       % Kb will be real, +ve
       Kb = sqrt(Kbsq);  % is used separately AND contained in mub
%        mub= Kb*(NB/abs(f)); % scalar, real +ve in N and S hemisphere
       zeta=h*Kb*NB/abs(f);
       KXgh  =  (hy*kta - hx*lta);  
          %LHS, RHS
       lhs(ja) = sig; %sig minimization, don't place outside ja loop
%        rhs(ja) = NB*KXgh/Kb/tanh(mub*h); 
       rhs(ja) = NB*KXgh/Kb/tanh(zeta);

   end %if Kbsq

end % ja loop that fills variables for lhs vs rhs optimization

    [~, iwant] = min(abs(rhs-lhs));  % min ignores NaNs
    k = kt(iwant);      % best k,l, either could be +ve or -ve
    l = lt(iwant);
%     sigfit=lhs(iwant);  %this is the exact original sig value
    sigfit=rhs(iwant);  %this is the exact fitted sigma value for k,l, Ei
%      not used: Tcheck=2*pi/86400/rhs(iwant);  % single value, not a vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if ifcheck==1     %SLOW      % ifcheck==1 makes a plot of crossings
%         rad2deg=180/pi;
%         figure(3); hold on;
%         plot(bang*rad2deg, lhs, 'k.', bang*rad2deg, rhs, 'r.', ...
%             bang(iwant)*rad2deg, lhs(iwant), 'ko', ...
%             bang(iwant)*rad2deg, rhs(iwant) ,'r*')
%         title('fitted crossings lhs and rhs in klcirclefitp')
%         bangmin=rad2deg*min(bang); bangmax=rad2deg*max(bang);
%         rhmin=min(rhs); rhmax=max(rhs);     %range of sig & lhs vals
%         axis([bangmin bangmax rhmin rhmax]); 
%         xlabel('angle'); 
%         ylabel('lhs=sig (blk), rhs=NB*KXgh/(Kb*tanh(mub*h) (red), (1/s)');
%         hold on
%     end % if ifcheck==1
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          

% BXz = [hy, -hx]; % 'B cross z' to the right of B
%                  % BXz is "topographic westward in N hemisphere


