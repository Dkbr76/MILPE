








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

% MILPE Projection
y1_MILPE = UyUxP * X;






% LSQ
a1 = L2(X,Y)

% LSQ Projection
y1_LSQ = a1 * X;
   





% DFT
[amp,phs,f] =  DFT(y1,dt);                      % exec DFT for y1, frequency might dependent on circumstances
y1_DFT      =  DFT_RECON(amp,phs,f,dt,m,12);    % DFT reconstruction 





% time history 
figure(11)
plot(t,y1,      'k'  ); hold on;   % y1 - OG
plot(t,y1_MILPE,'k--');            % y1 - MILPE 
plot(t,y1_LSQ,  'r'  );            % y1 - LSQ
plot(t,y1_DFT,  'y--');            % y1 - DFT 
xlabel('t');
ylabel('vars');
legend("OG","MILPE","LSQ","DFT");

aaf
