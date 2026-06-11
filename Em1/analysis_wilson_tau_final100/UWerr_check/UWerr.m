function [value,dvalue,ddvalue,tauint,dtauint,Qval, Wopt] = ...
UWerr(Data,Stau,Nrep,Name,Quantity,varargin)
%
% call:
% [value,dvalue,ddvalue,tauint,dtauint,Qval] = ...
% UWerr(Data,Stau,Nrep,Name,Quantity,P1,P2,....)
%
% autocorrelation-analysis of MC time-series
% following the Gamma-method in
% ``Monte Carlo errors with less errors''
% by Ulli Wolff, hep-lat/0306017
% 
%----------------------------------------------------------------------------
%  Ulli Wolff,   Nov. 2006 Version V6
%  new: errobars on normalized autocorrelation (Luescher formula, see below)
%       subleading terms in error of tauint changed again (see v4 of paper)
%       protection against error^2 < 0 estimates 
% Changes by Tomasz Korzec, 2016:
%       - replacement of mlock/munlock mechanism
%       - additional windowing conditions:
%           2) Sum to at least t=MIN_W_OPT
%           3) tau_int(Wopt) and tau_int(2*Wopt) have to be compatible
%       - negative Stau: -Stau is used as Wopt
%       - compatibility with octave
%----------------------------------------------------------------------------
% input arguments beyond Data maybe omitted or set to []
% in the input list; they then get default value indicated in [D=..]
%
%
% Data     -- (N x Nalpha) matrix of measured (equilibrium!) data 
%              N = total number of measurements
%              Nalpha = number of (primary) observables
%              if your data are in a different format, consider 
%              the MATLAB commands `reshape' and `permute'
%
% Stau     -- guess for the ratio S of tau/tauint [D=1.5]
%             if set = 0, absence of autocorrelations is *assumed*
%             if negative: |Stau| is used as summation window
%
% Nrep     -- vector [N1 N2 N3 ...] specifying a breakup of the N rows
%             of Data into replica of length N1,N2,.. (N=sum(Nrep)!)  [D=[N]]
%             The replica distribution is histogrammed and a Q-value 
%             (=probability to find this much or more scatter) is given 
%             for R >= 2 replica; if you have one history you may set Nrep to 
%             artificially split it into >=2 replica to see their distribution                
%
% Name     -- if string: name of observable for titles of generated plots
%             if not string: all plots are supressed [D='NoName']
%
% Quantity -- either:
%             scalar function handle (@functionname) for the derived 
%             observable F; it has to operate on a row-vector of length 
%             Nalpha as first argument; optional parameters P1,P2,... are 
%             passed on to this function as 2nd, 3rd, .... argument
%          -- or:
%             integer between 1 and Nalpha to analyze primary observable [D=1]
%
%----------------------------------------------------------------------------
% output::
%  value   -- estimate for F
%  dvalue  -- its error (incl. autocorrelation effects)
%  ddvalue -- statistical error of the error
%  tauint  -- integrated autocorrelation time
%  dtauint -- statistical error of tauint
%  Qval    -- Q-value of the replica distribution if R >= 2
%             (goodness of fit to a constant)
%----------------------------------------------------------------------------
%

% hard-wired parameters
MEAN_BIAS_MIN_REP = 2;   % minimal number of replica that has to be present
                         % in order for bias cancelation of the mean to be
                         % applied
ERROR_BIAS_MIN_REP = 2;  % minimal number of replica that has to be present
                         % in order for bias cancelation of the error to be
                         % applied
DERIVATIVE_H       = 1e-5;  % h in numerical derivatives is std*DERIVATIVE_H
MIN_W_OPT          = 5;  % Wopt cannot be smaller than this

% analyze input arguments:

[N,Nalpha]=size(Data); % number of measurements and observables

if nargin < 5 || isempty(Quantity)  Quantity=1;    end
if nargin < 4 || isempty(Name),     Name='NoName'; end
if nargin < 3 || isempty(Nrep),     Nrep=N;        end
if nargin < 2 || isempty(Stau),     Stau=1.5;      end

if any(Nrep ~= round(Nrep)) || any(Nrep < 1) || sum(Nrep) ~= N
  error('inconsistent N,Nrep')
end

if isnumeric(Quantity)
   if round(Quantity) ~= Quantity || Quantity < 1 || Quantity > Nalpha
       error('illegal numeric value of Quantity')
   end
   primary = 1; % logical variable
   primaryindex=Quantity;
else
   primary = 0;
end

Nrep=Nrep(:);   % make column
R=length(Nrep); % number of Replica

%----------------------------------------------------------------------------

% means of primaries:

abb=mean(Data,1);    %    total mean, primary obs.
abr=zeros(R,Nalpha); % replicum-mean, primary obs.
i0=1;
for r=1:R
  i1=i0-1+Nrep(r);
  abr(r,:)=mean(Data(i0:i1,:),1);
  i0=i0+Nrep(r);
end

% means of derived (or primary, depending on Quantity):

if primary
  Fbb=abb(primaryindex);
  Fbr=abr(:,primaryindex);	
else
  Fbb=feval(Quantity,abb,varargin{:}); % total mean, derived obs.
  Fbr=zeros(R,1);                      % replica-means, derived obs.
  for i=1:R
    Fbr(i)=feval(Quantity,abr(i,:),varargin{:});
  end
end
Fb=sum(Fbr.*Nrep)/N;  % weighed mean of replica means

%----------------------------------------------------------------------------

% form the gradient of f and project fluctuations:

if primary
  delpro = Data(:,primaryindex)-abb(primaryindex);
else
  fgrad=zeros(Nalpha,1);
  h=std(Data,1)/sqrt(N)*DERIVATIVE_H; % spread for num. derivative
  ainc=abb;
  for alpha=1:Nalpha
    if h(alpha) == 0     % Data for this observable do not fluctuate
      fgrad(alpha)=0;
    else
      ainc(alpha) =abb(alpha)+h(alpha);
      fgrad(alpha)=feval(Quantity,ainc,varargin{:});
      ainc(alpha) =abb(alpha)-h(alpha);
      fgrad(alpha)=fgrad(alpha)-feval(Quantity,ainc,varargin{:});
      ainc(alpha) =abb(alpha);
      fgrad(alpha)=fgrad(alpha)/(2*h(alpha));
    end
  end

% projected deviations: 
  delpro = Data*fgrad-abb*fgrad;
end

%----------------------------------------------------------------------------
% This is returned if the error analysis fails in the following:
value   = Fbb;
dvalue  = 0;
ddvalue = 0;
tauint  = 0.5;
dtauint = 0;
Qval=[];
  
%----------------------------------------------------------------------------  
% Computation of Gamma, automatic windowing
% note: W=1,2,3,.. in Matlab-vectorindex <-> W=0,1,2,.. in the paper!

% values for W=0:
GammaFbb(1) = mean(delpro.^2);
% sick case:
if GammaFbb(1) == 0
  %beep
  fprintf('\n WARNING: no fluctuations \n')
  return
end

% compute Gamma(t), find the optimal W

if Stau == 0                    % no autocorrelations assumed
  Wopt=0;
  tmax=0; 
  flag=false;
elseif Stau<0                   % using -Stau as Wopt
  Wopt=-Stau;
  tmax=min(min(Nrep),2*Wopt);
  if tmax<Wopt
     fprintf('\nWARNING: sequence too short to reach Wopt=%d, using Wopt=%d\n',Wopt,tmax-1);
     Wopt=tmax-1;
  end
  flag=false;
else
  tmax = floor(min(Nrep)/2);    % do not take t larger than this
  flag=true; 
  Gint=0;
  if(tmax < MIN_W_OPT+1)
     fprintf('\nERROR: sequence too short to estimate errors reliably\n')
     return
  end
end
t=1;
t2=1;
while t <= tmax,
  % compute rho and tauint_tmp up to 2*t
  for tt=t2:min(2*t,tmax)
     GammaFbb(tt+1) = 0;
     % sum over replica:
     i0=1;
     for r=1:R
       i1=i0-1+Nrep(r);
       GammaFbb(tt+1) = GammaFbb(tt+1) + sum(delpro(i0:i1-tt).*delpro(i0+tt:i1));
       i0=i0+Nrep(r);
     end
     GammaFbb(tt+1) = GammaFbb(tt+1)/(N-R*tt);
  end
  t2=min(2*t,tmax)+1; % last computed value of GammaFbb
  
  % temporary values for tauint
  rho=GammaFbb/GammaFbb(1);
  tauint_tmp=cumsum(rho)-0.5;

  if flag
    Gint=Gint+GammaFbb(t+1)/GammaFbb(1);
    if Gint <= 0, tauW=eps; else tauW = Stau/(log((Gint+1)/Gint)); end
    gW = exp(-t/tauW)-tauW/sqrt(t*N);
    if gW < 0                % Windowing condition 1 satisfied
      Wopt=t; 
      if t>= MIN_W_OPT       % Windowing condition 2 satisfied
         dtauint_tmp = (tauint_tmp(Wopt+1)*2*sqrt((Wopt-tauint_tmp(Wopt+1)+0.5)/N));
         if abs(tauint_tmp(t) - tauint_tmp(min(end,2*t))) < dtauint_tmp % Windowing condition 3 satisfied
            tmax=min(tmax,2*t); 
            flag = false; % Gamma up to tmax
         end
      end
    end
  end
  t=t+1;
end   % while-loop 

if flag,
  %beep
  fprintf('\n WARNING: windowing condition failed up to W = %i \n',tmax)
  Wopt=tmax;
end

CFbbopt=GammaFbb(1) + 2*sum(GammaFbb(2:Wopt+1));   % first estimate
if CFbbopt <= 0
  beep
  fprintf('\n WARNING: Gamma pathological: estimated error^2 < 0 \n')
  return
end

if R>=ERROR_BIAS_MIN_REP
   GammaFbb=GammaFbb+CFbbopt/N;                    % bias in Gamma corrected
end
CFbbopt=GammaFbb(1) + 2*sum(GammaFbb(2:Wopt+1));   % refined estimate
sigmaF=sqrt(CFbbopt/N);                            % error of F
rho=GammaFbb/GammaFbb(1);                          % normalized autocorr.
tauintFbb=cumsum(rho)-0.5;

%----------------------------------------------------------------------------

% bias cancellation for the mean value

if R >= MEAN_BIAS_MIN_REP
  bF = (Fb-Fbb)/(R-1);
  Fbb=Fbb - bF;
  if abs(bF) > sigmaF/4
    beep
    fprintf('\n WARNING: a %.1f sigma bias of the mean has been cancelled \n',bF/sigmaF)
  end
  Fbr = Fbr - bF*N./Nrep;
  Fb  = Fb  - bF*R;
end

%----------------------------------------------------------------------------

% answers to be returned:

value   = Fbb;
dvalue  = sigmaF;
ddvalue = dvalue*sqrt((Wopt+0.5)/N);
tauint  = tauintFbb(Wopt+1);
% changed subleading terms (prev. versions!) here:
dtauint = tauint*2*sqrt((Wopt-tauint+0.5)/N);

% Q-value for replica distribution if R >= 2:
 
if R >= 2
  chisq=sum((Fbr-Fb).^2.*Nrep)/CFbbopt;
  Qval=1-gammainc(chisq/2,(R-1)/2);
else
  Qval=[];
end

%----------------------------------------------------------------------------

% Plotting:

% make plots, if demanded:

if ~ ischar(Name)   % no plots
   return 
end

% find old UWerr figures
f1=[]; f2=[]; f3=[];
figHandles = findobj('Type','figure');
for l=1:length(figHandles)
   str = get(figHandles(l),'Name');
   if strcmp(str,'UWerr 1')
      f1 = figHandles(l);
   end
   if strcmp(str,'UWerr 2')
      f2 = figHandles(l);
   end
   if strcmp(str,'UWerr 3')
      f3 = figHandles(l);
   end
end

scrsz=get(0,'screensize'); % to place plots on the screen

% autocorrelation:

if Stau ~= 0

  if isempty(f1), 
    f1=figure;
    set(f1,'Name','UWerr 1');
    set(f1,'Position',[0.01*scrsz(3) 0.45*scrsz(4) 0.45*scrsz(3) 0.45*scrsz(4)]);
    set(f1,'DefaultAxesFontSize',16);
  else 
    figure(f1)
  end; 
  clf;

% plot of GammaF(t)/Gamma(0)=rho(t)
% contruct errors acc. to hep-lat/0409106 eq. (E.11)
% pad zeros to simplify summation:
  rho(tmax+2:2*tmax+Wopt+1)=0;

  for t=0:tmax
    k=[max(1,t-Wopt):t+Wopt];                          % summation range    
    drho(t+1)=sum((rho(k+t+1)+rho(abs(k-t)+1)-2*rho(t+1)*rho(k+1)).^2);
    drho(t+1)=sqrt(drho(t+1)/N);
  end
  subplot(2,1,1)
  errorbar(0:tmax,rho(1:tmax+1),drho);
  v=axis;v=[0 tmax min(v(3),-0.1) 1];axis(v);hold on;
  plot([0 tmax],[0 0],'--');                           % zero line
  plot([Wopt Wopt],v(3:4),'r-');                       % Bar at Wopt
  title(['normalized autocorrelation of ' Name])
  ylabel('\rho')

% plot of tauintF

  subplot(2,1,2)
% tauint vs. W:
  errorbar(0:tmax,tauintFbb,tauintFbb.*sqrt(([0:tmax])/N)*2); 
  v=axis;v(1:2)=[0 tmax];axis(v);hold on;
  plot([Wopt Wopt],v(3:4),'r-');               % Bar at Wopt
  plot([0 tmax],[tauint tauint],'r--')         % hor. line at estimate
  title(['\tau_{int} with statistical errors of ' Name]);
  ylabel('\tau_{int}')
  xlabel('W')
  drawnow;
else
  if ~isempty(f1)
     close(f1);
  end
end

% histogram of data, if primary:

if primary
  if isempty(f2)
    f2=figure; 
    set(f2,'Name','UWerr 2');
    set(f2,'Position',[0.51*scrsz(3) 0.52*scrsz(4) 0.45*scrsz(3) 0.35*scrsz(4)]);
    set(f2,'DefaultAxesFontSize',16);
  else 
    figure(f2) 
  end
  clf
  hist(Data(:,primaryindex),20);
  title(['Distribution of all data for ' Name]);
  drawnow;
else
  if ~isempty(f2)
     close(f2);
  end
end

% histogram of replica values if R >= 2:

if R >= 2
  p=(Fbr-Fb)./(sigmaF*sqrt(N./Nrep-1));  
  if isempty(f3), 
     f3=figure;
     set(f3,'Name','UWerr 3');
  else
      figure(f3);
  end;
  clf;
  set(f3,'Position',[0.51*scrsz(3) 0.03*scrsz(4) 0.45*scrsz(3) 0.35*scrsz(4)]);
  set(f3,'DefaultAxesFontSize',16);
  hist(p);
  title(['replica distribution (mean=0,var=1) for ' Name ' >> Q = ' num2str(Qval,'%.2g') ' <<'])
  drawnow;
else
  if ~isempty(f3)
     close(f3);
  end
end

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
