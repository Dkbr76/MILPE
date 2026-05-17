








% MILPE principle - 1
% 2 0 2 6 . 0 5 
% Reconstruct cos(t) with sin(t): 
%   - same as triagonal identity and probably related to time-delay problem?

clear all; close all; format long;






% ctrl
te = 1e+1;                  % simulation duration
dt = 1e-2;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots
t  = [0:m-1]*dt;            % time-vector

% input time-history(s)
x1 = sin(t);                % <- baseline
x2 = sin(t + 1e-12);        % variant
x3 = sin(t + 1e-10);        % variant
x4 = sin(t + 1e-8);        % variant
x5 = sin(t + 0.1);        % variant

% output time-history(s)
y1 = cos(t);

% input-subspace
X  = [x1;x2;x3;x4;x5];

% output-subspace
Y  = [y1];

% MILPE
[UyUxP,Ux,Uy,U,S] = MILPE2(X,Y,0);

% MILPE projection
y1_MILPE = UyUxP * X;

% LSQ projection
a1      =  L2(X,Y);
y1_LSQ  =  a1 * X;

% verification
U
S
Ux
UxP=pinv(Ux)
Uy
UyUxP
a1

% plot
figure(1);  
plot(t,y1,'g','LineWidth',2); hold on;   % OG - y
plot(t,y1_MILPE,'r--','LineWidth',2);    % MILPE - y
plot(t,y1_LSQ,'b--','LineWidth',2);      % LSQ - y
xlabel('t'); ylabel('vars');
legend("OG","MILPE","LSQ");
aaf;


















