








% MILPE principle - 1
% 2 0 2 6 . 0 5 
% Reconstruct cos(t) with sin(t): 
%   - same as triagonal identity and probably related to time-delay problem?

clear all; close all;






% ctrl
te = 1e+1;                  % simulation duration
dt = 1e-2;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots
t  = [0:m-1]*dt;            % time-vector

% input time-history(s)
x1 = sin(t);                % <- baseline
x2 = sin(t + 1e-8);        % variant

% output time-history(s)
y1 = cos(t);

% input-subspace
X  = [x1;x2];

% output-subspace
Y  = [y1];

% MILPE
UyUxP = MILPE(X,Y,0)

% MILPE projection
y1_MILPE = UyUxP * X;






% LSQ projection
a1      =  L2(X,Y)
y1_LSQ  =  a1 * X;
   





% time history - x-position
figure(11)
plot(t,y1,'k'); hold on;   % OG - y
plot(t,y1_MILPE,'y--');    % MILPE - y
plot(t,y1_LSQ,'r--');      % LSQ - y
xlabel('t');
ylabel('vars');
legend("OG","MILPE","LSQ")

aaf




















