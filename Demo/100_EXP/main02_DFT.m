








% MILPE principle - 2
% 2 0 2 6 . 0 5 
% Reconstruct step-function with cos(t)

clear all; close all;






% ctrl
te = 2*pi;                  % simulation duration
dt = 1e-5;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots

% time-vector
t  = [0:m-1]*dt;

% input time-history(s)
x1  = cos(t* 1);              
x2  = cos(t* 3);       
x3  = cos(t* 5);              
x4  = cos(t* 7);              
x5  = cos(t* 9);              
x6  = cos(t*11);              

% output time-history(s) - choose one from examples or try adding your own harmonics 
% y1  = 2*cos(t) - 0.5*cos(3*t) + 0.25*cos(7*t);  % random harmonic function 1
y1  = sign(cos(t));                             % step function

% input-subspace
X  =  [ x1; x2; x3; x4; x5; x6 ];

% output-subspace
Y  =  [ y1 ];

% MILPE
UyUxP = MILPE(X,Y,0)

% sav time-history
TH0 = [t;X;Y];






% MILPE - Projection
y1_MILPE = UyUxP * X;

   



% verification - compare with DFT result
[amp,phs,f] =  DFT(y1,dt);                      % exec DFT for y1, frequency might dependent on circumstances
y1_DFT      =  DFT_RECON(amp,phs,f,dt,m,12);    % DFT reconstruction 





% time history 
figure(11)
plot(TH0(1,:),TH0(  2,:),'c'); hold on; % OG - x1
plot(TH0(1,:),TH0(  3,:),'c'); hold on; % OG - x2
plot(TH0(1,:),TH0(  4,:),'c'); hold on; % OG - x3
plot(TH0(1,:),TH0(  5,:),'c'); hold on; % OG - x4
plot(TH0(1,:),TH0(  6,:),'c'); hold on; % OG - x5
plot(TH0(1,:),TH0(  7,:),'c'); hold on; % OG - x6
plot(TH0(1,:),TH0(end,:),'k'); hold on; % OG - y1
plot(t,y1_MILPE,'k--');        % MILPE - y1
plot(t,y1_DFT,'m--');            % DFT - y1
xlabel('t');
ylabel('vars');

aaf
