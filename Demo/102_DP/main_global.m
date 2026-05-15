








% MILPE demo - Data-driven Double-pendulum prediction 
% Mar.2026 https://www.mdpi.com/3762868

clear all; close all;





%% OG Double-pendulum simulation 1
% purpose: data aquisition for the eigenvector extraction
% reference: Shinbrot et al. Chaos in a double pendulum. Am. J. Phys. 1992, 60, 491–499
disp("OG DP Sim #1");
te    =  1e+0;                   % [s] simulation duration
dt    =  1e-5;                  % [s] time-step
ic    =  [10 10 0 0];           % initial condition [theta1, theta2, theta1_dot, theta2_dot]
co    =  [1 1 1 1 9.8];         % DP parameters [L1, L2, m1, m2, g]
TH0   =  DP(te,dt,ic,co);       % Run DP func and obtain time-history


%% OG Double-pendulum simulation 2
% purpose: simulate DP longer time for the comparison against MILPE prediction
disp("OG DP Sim #2");
te    =  20.0;                  % [s] simulation duration
dt    =  1e-5;                  % [s] time-step
ic    =  [70 50 0 0];           % initial condition [theta1, theta2, theta1_dot, theta2_dot]
co    =  [1 1 1 1 9.8];         % DP parameters [L1, L2, m1, m2, g]
TH1   =  DP(te,dt,ic,co);       % Run DP func and obtain time-history



%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation
disp("Running MILPE Algorithm");

% var alloc from ref time-history used for eigenvector extraction
t1   = TH0(2,:);
t2   = TH0(3,:);
t1d  = TH0(4,:);
t2d  = TH0(5,:);
t1dd = TH0(6,:);
t2dd = TH0(7,:);

% input subspace X1
DENOM   = 3-cos(2*t1-2*t2);
X11     = sin(t1)                           ./DENOM;
X12     = sin(t1-2*t2)                      ./DENOM;
X13     = sin(t1-t2).*t2d.^2                ./DENOM;
X14     = sin(t1-t2).*t1d.^2.*cos(t1-t2)    ./DENOM;
X1      = [X11; X12; X13; X14]; 

% input subspace X2
X21     = sin(t1-t2).*t1d.^2                ./DENOM;
X22     = sin(t1-t2).*cos(t1)               ./DENOM;
X23     = sin(t1-t2).*t2d.^2.*cos(t1-t2)    ./DENOM;
X2      = [X21; X22; X23]; 

% output subspace Y1
Y1 = [t1dd];

% output subspace Y2
Y2 = [t2dd];

% MILPE
U1yU1xP = MILPE(X1,Y1,0); % MILPE low-rank governing equation (Uy*Ux+)
U2yU2xP = MILPE(X2,Y2,0); % MILPE low-rank governing equation (Uy*Ux+)






% LSQ
a1 = L2(X1,Y1);
a2 = L2(X2,Y2);






%% MILPE prediction
% purpose: Predict Double-pendulum with approximated governing equation

% controller
te    =  20;                    % [s] simulation duration
dt    =  1e-5;                  % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots

% factors
d2r   =  pi/180.0;              % multiplication factor for deg to rad conversion
r2d   =  180.0/pi;              % multiplication factor for rad to deg conversion  

% initial condition (put degree value before d2r)
t1   =  70*d2r;                 % [rad]   theta1, angle of the 1st pendulum
t2   =  50*d2r;                 % [rad]   theta2, angle of the 2nd pendulum     
t1d  =   0*d2r;                 % [rad/s] theta1_dot, angular velocity of the 1st pendulum 
t2d  =   0*d2r;                 % [rad/s] theta2_dot, angular velocity of the 2nd pendulum

% time-loop
for it=1:m

    % echo
    if ( it == 1 ) cnt = 20; end % counter
    if ( mod(it,0.2*m) < 1 ) disp("MILPE......"+cnt+"%"); cnt = cnt+20; end

    % input snapshot X1
    DENOM   = 3-cos(2*t1-2*t2);
    X11     = sin(t1)                        /DENOM;
    X12     = sin(t1-2*t2)                   /DENOM;
    X13     = sin(t1-t2)*t2d^2               /DENOM;
    X14     = sin(t1-t2)*t1d^2*cos(t1-t2)    /DENOM;
    X1      = [X11 X12 X13 X14];

    % input snapshot X2
    X21     = sin(t1-t2)*t1d^2               /DENOM;
    X22     = sin(t1-t2)*cos(t1)             /DENOM;
    X23     = sin(t1-t2)*t2d^2*cos(t1-t2)    /DENOM;
    X2      = [X21 X22 X23];

    % output snapshot Y1
    Y1 = U1yU1xP * X1';

    % output snapshot Y2
    Y2 = U2yU2xP * X2';

    % var alloc
    t1dd = Y1; 
    t2dd = Y2;
    
    % time-advance
    t1d = t1d + t1dd*dt;
    t2d = t2d + t2dd*dt;
    t1  = t1  +  t1d*dt;
    t2  = t2  +  t2d*dt;

    % store snapshot info
    SNP = [it*dt t1 t2 t1d t2d t1dd t2dd]';

    % init R0 at it==1
    if ( it == 1 ) TH2 = zeros(size(SNP,1),m); end

    % save as time-history
    TH2(:,it) = SNP;

end





%% Verification 1 - governing equation 

U1yU1xP
U2yU2xP

a1
a2




%% Verification 2 - time history comparison

% time history - trajectory
figure(1)
plot(TH1(2,:)*r2d,TH1(3,:)*r2d,'k'); hold on; % OG
plot(TH2(2,:)*r2d,TH2(3,:)*r2d,'r');          % MILPE
title("trajectory"); axis equal;
xlabel('\theta_1(deg)');
ylabel('\theta_2(deg)');
legend("OG","MILPE");

% time history - theta1
figure(11)
plot(TH1(1,:),TH1(2,:)*r2d,'k'); hold on; % OG
plot(TH2(1,:),TH2(2,:)*r2d,'r');          % MILPE
title("time history of \theta_1");
xlabel('t(s)');
ylabel('\theta_1(deg)');
legend("OG","MILPE");

% time history - theta2
figure(12)
plot(TH1(1,:),TH1(3,:)*r2d,'k'); hold on; % OG
plot(TH2(1,:),TH2(3,:)*r2d,'r');          % MILPE
title("time history of \theta_2");
xlabel('t(s)');
ylabel('\theta_2(deg)');
legend("OG","MILPE");

aaf

