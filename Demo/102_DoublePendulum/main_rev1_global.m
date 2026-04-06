








% MILPE demo - Data-driven Double-pendulum prediction 
% Mar.2026 https://www.mdpi.com/3762868

clear all; close all;






%% OG Double-pendulum simulation 1
% purpose: data aquisition for the eigenvector extraction
% reference: Shinbrot et al. Chaos in a double pendulum. Am. J. Phys. 1992, 60, 491–499

% controller
te    =  3.0;                   % [s] simulation duration
dt    =  0.00001;               % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots

% factors
d2r   =  pi/180.0;              % multiplication factor for deg to rad conversion
r2d   =  180.0/pi;              % multiplication factor for rad to deg conversion  

% parameters 
L1   =  1.0;                    % [m]   length of 1st arm
L2   =  1.0;                    % [m]   length of 2nd arm
m1   =  1;                      % [kg]  mass of 1st pendulum
m2   =  1;                      % [kg]  mass of 2nd pendulum
g    =  9.8;                    % [m/s2] gravitational acceleration

% initial condition (put degree value before d2r)
t1   =  10*d2r;                 % [rad]   theta1, angle of the 1st pendulum
t2   =  10*d2r;                 % [rad]   theta2, angle of the 2nd pendulum     
t1d  =   0*d2r;                 % [rad/s] theta1_dot, angular velocity of the 1st pendulum 
t2d  =   0*d2r;                 % [rad/s] theta2_dot, angular velocity of the 2nd pendulum

% time-loop
for it=1:m

    % echo
    if ( mod(it,10000) == 0 ) disp(it);  end

    % Double-Pendulum 
    a1  =  -g*(2*m1+m2)*sin(t1);
    b1  =  -m2*g*sin(t1-2*t2);
    c1  =  -2*sin(t1-t2)*m2*(t2d^2*L2);
    d1  =  -2*sin(t1-t2)*m2*(t1d^2*L1*cos(t1-t2));
    e1  =  L1*(2*m1+m2-m2*cos(2*t1-2*t2));

    a2  =  2*sin(t1-t2)*t1d^2*L1*(m1+m1);
    b2  =  2*sin(t1-t2)*g*(m1+m2)*cos(t1);
    c2  =  2*sin(t1-t2)*t2d^2*L2*m2*cos(t1-t2);
    e2  =  L2*(2*m1+m2-m2*cos(2*t1-2*t2));
    
    t1dd  = (a1+b1+c1+d1)/e1;
    t2dd  = (a2+b2+c2)/e2;

    % Time marching (1st Euler)
    t1d = t1d + t1dd*dt;
    t2d = t2d + t2dd*dt;
    t1  = t1  + t1d *dt;
    t2  = t2  + t2d *dt;
    
    % store snapshot info
    R = [it*dt t1 t2 t1d t2d t1dd t2dd]';

    % init R0 at it==1
    if ( it == 1 ) R0 = zeros(size(R,1),m); end

    % save as time-history
    R0(:,it) = R;

end





%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation

% var alloc (committed overwriting var names)
t1   = R0(2,:);
t2   = R0(3,:);
t1d  = R0(4,:);
t2d  = R0(5,:);
t1dd = R0(6,:);
t2dd = R0(7,:);

% input subspace X1
DENOM   = 3-cos(2*t1-2*t2);
X11     = sin(t1)                           ./DENOM;
X12     = sin(t1-2*t2)                      ./DENOM;
X13     = sin(t1-t2).*t2d.^2                ./DENOM;
X14     = sin(t1-t2).*t1d.^2.*cos(t1-t2)    ./DENOM;

X1 = [X11; X12; X13; X14]; 

% input subspace X2
DENOM   = 3-cos(2*t1-2*t2);
X21     = sin(t1-t2).*t1d.^2                ./DENOM;
X22     = sin(t1-t2).*cos(t1)               ./DENOM;
X23     = sin(t1-t2).*t2d.^2.*cos(t1-t2)    ./DENOM;

X2 = [X21; X22; X23]; 

% output subspace Y1
Y1 = [t1dd];

% output subspace Y2
Y2 = [t2dd];

% unified solution space Z1
Z1 = [X1;Y1];

% unified solution space Z2
Z2 = [X2;Y2];

% SVD for eigenvector extraction
[U1, S1, V1] = svd(Z1,'econ');
[U2, S2, V2] = svd(Z2,'econ');

% governing equations
nX1      =  size(X1,1);                % number of variables in input subspace 
U1x      =  U1(1:nX1     , 1:nX1);     % projection matrix on input subspace   
U1y      =  U1(nX1+1:end , 1:nX1);     % projection matrix on output subspace  
U1yU1xP  =  U1y*pinv(U1x);             % MILPE low-rank governing equation (Uy*Ux+)

nX2      =  size(X2,1);                
U2x      =  U2(1:nX2     , 1:nX2);     
U2y      =  U2(nX2+1:end , 1:nX2);     
U2yU2xP  =  U2y*pinv(U2x);             





%% OG Double-pendulum simulation 2
% purpose: Simulate OG Double-pendulum longer time for the comparison against MILPE prediction

% controller
te    =  20;                    % [s] simulation duration
dt    =  0.00001;               % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots

% factors
d2r   =  pi/180.0;              % multiplication factor for deg to rad conversion
r2d   =  180.0/pi;              % multiplication factor for rad to deg conversion  

% parameters 
L1   =  1.0;                    % [m]   length of 1st arm
L2   =  1.0;                    % [m]   length of 2nd arm
m1   =  1;                      % [kg]  mass of 1st pendulum
m2   =  1;                      % [kg]  mass of 2nd pendulum
g    =  9.8;                    % [m/s2] gravitational acceleration

% initial condition (put degree value before d2r)
t1   =  70*d2r;                 % [rad]   theta1, angle of the 1st pendulum
t2   =  50*d2r;                 % [rad]   theta2, angle of the 2nd pendulum     
t1d  =   0*d2r;                 % [rad/s] theta1_dot, angular velocity of the 1st pendulum 
t2d  =   0*d2r;                 % [rad/s] theta2_dot, angular velocity of the 2nd pendulum

% time-loop
for it=1:m

    % echo
    if ( mod(it,10000) == 0 ) disp(it);  end

    % Double-Pendulum 
    a1  =  -g*(2*m1+m2)*sin(t1);
    b1  =  -m2*g*sin(t1-2*t2);
    c1  =  -2*sin(t1-t2)*m2*(t2d^2*L2);
    d1  =  -2*sin(t1-t2)*m2*(t1d^2*L1*cos(t1-t2));
    e1  =  L1*(2*m1+m2-m2*cos(2*t1-2*t2));

    a2  =  2*sin(t1-t2)*t1d^2*L1*(m1+m1);
    b2  =  2*sin(t1-t2)*g*(m1+m2)*cos(t1);
    c2  =  2*sin(t1-t2)*t2d^2*L2*m2*cos(t1-t2);
    e2  =  L2*(2*m1+m2-m2*cos(2*t1-2*t2));
    
    t1dd  = (a1+b1+c1+d1)/e1;
    t2dd  = (a2+b2+c2)/e2;

    % Time marching (1st Euler)
    t1d = t1d + t1dd*dt;
    t2d = t2d + t2dd*dt;
    t1  = t1  + t1d *dt;
    t2  = t2  + t2d *dt;
    
    % store snapshot info
    R = [it*dt t1 t2 t1d t2d t1dd t2dd]';

    % init R0 at it==1
    if ( it == 1 ) R1 = zeros(size(R,1),m); end

    % save as time-history
    R1(:,it) = R;

end





%% MILPE prediction
% purpose: Predict Double-pendulum with approximated governing equation

% controller
te    =  20;                    % [s] simulation duration
dt    =  0.00001;               % [s] time-step
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
    if ( mod(it,10000) == 0 ) disp(it);  end

    % input snapshot X1
    DENOM   = 3-cos(2*t1-2*t2);
    X11     = sin(t1)                        /DENOM;
    X12     = sin(t1-2*t2)                   /DENOM;
    X13     = sin(t1-t2)*t2d^2               /DENOM;
    X14     = sin(t1-t2)*t1d^2*cos(t1-t2)    /DENOM;

    X1 = [X11 X12 X13 X14];

    % input snapshot X2
    DENOM   = 3-cos(2*t1-2*t2);
    X21     = sin(t1-t2)*t1d^2               /DENOM;
    X22     = sin(t1-t2)*cos(t1)             /DENOM;
    X23     = sin(t1-t2)*t2d^2*cos(t1-t2)    /DENOM;

    X2 = [X21 X22 X23];

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
    R = [it*dt t1 t2 t1d t2d t1dd t2dd]';

    % init R0 at it==1
    if ( it == 1 ) R2 = zeros(size(R,1),m); end

    % save as time-history
    R2(:,it) = R;

end





%% Verification 1 - governing equation 

U1yU1xP
U2yU2xP





%% Verification 2 - time history comparison

% time history - trajectory
figure(1)
plot(R1(2,:),R1(3,:),'k'); hold on; % OG
plot(R2(2,:),R2(3,:),'r');          % MILPE
xlabel('t1');
ylabel('t2');
legend("OG","MILPE");

% time history - theta1
figure(11)
plot(R1(1,:),R1(2,:),'k'); hold on; % OG
plot(R2(1,:),R2(2,:),'r');          % MILPE
xlabel('t');
ylabel('t1');
legend("OG","MILPE");

% time history - theta2
figure(12)
plot(R1(1,:),R1(3,:),'k'); hold on; % OG
plot(R2(1,:),R2(3,:),'r');          % MILPE
xlabel('t');
ylabel('t2');
legend("OG","MILPE");



