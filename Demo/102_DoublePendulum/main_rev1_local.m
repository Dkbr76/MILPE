








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
t1   =  30*d2r;                 % [rad]   theta1, angle of the 1st pendulum
t2   =  30*d2r;                 % [rad]   theta2, angle of the 2nd pendulum     
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
    R = [it*dt t1 t2 t1dd t2dd]';

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
t1dd = R0(4,:);
t2dd = R0(5,:);

% input subspace
X = [t1; t2; t1.^3; t2.^3]; 

% output subspace
Y = [t1dd;t2dd];

% unified solution space (Z)
Z = [X;Y];

% SVD for eigenvector extraction
[U, S, V] = svd(Z,'econ');

% governing equation
nX       =  size(X,1);              % number of variables in input subspace 
Ux       =  U(1:nX     , 1:nX);     % projection matrix on input subspace   
Uy       =  U(nX+1:end , 1:nX);     % projection matrix on output subspace  
UyUxP    =  Uy*pinv(Ux);            % MILPE low-rank governing equation (Uy*Ux+)





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
t1   =  30*d2r;                 % [rad]   theta1, angle of the 1st pendulum
t2   =  30*d2r;                 % [rad]   theta2, angle of the 2nd pendulum     
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
    R = [it*dt t1 t2 t1dd t2dd]';

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
t1   =  30*d2r;                 % [rad]   theta1, angle of the 1st pendulum
t2   =  30*d2r;                 % [rad]   theta2, angle of the 2nd pendulum     
t1d  =   0*d2r;                 % [rad/s] theta1_dot, angular velocity of the 1st pendulum 
t2d  =   0*d2r;                 % [rad/s] theta2_dot, angular velocity of the 2nd pendulum

% time-loop
for it=1:m

    % echo
    if ( mod(it,10000) == 0 ) disp(it);  end

    % input snapshot X
    X = [t1 t2 t1^3 t2^3];

    % output snapshot Y
    Y = UyUxP * X';

    % var alloc
    t1dd = Y(1); 
    t2dd = Y(2);
    
    % time-advance
    t1d = t1d + t1dd*dt;
    t2d = t2d + t2dd*dt;
    t1  = t1  +  t1d*dt;
    t2  = t2  +  t2d*dt;

    % store snapshot info
    R = [it*dt t1 t2 t1dd t2dd]';

    % init R0 at it==1
    if ( it == 1 ) R2 = zeros(size(R,1),m); end

    % save as time-history
    R2(:,it) = R;

end





%% Verification 1 - governing equation 

UyUxP





%% Verification 2 - time history comparison

% time history - trajectory
figure(1)
plot(R1(2,:)*r2d,R1(3,:)*r2d,'k'); hold on; % OG
plot(R2(2,:)*r2d,R2(3,:)*r2d,'r');          % MILPE
title("trajectory"); axis equal;
xlabel('\theta_1(deg)');
ylabel('\theta_2(deg)');
legend("OG","MILPE");

% time history - theta1
figure(11)
plot(R1(1,:),R1(2,:)*r2d,'k'); hold on; % OG
plot(R2(1,:),R2(2,:)*r2d,'r');          % MILPE
title("time history of \theta_1");
xlabel('t(s)');
ylabel('\theta_1(deg)');
legend("OG","MILPE");

% time history - theta2
figure(12)
plot(R1(1,:),R1(3,:)*r2d,'k'); hold on; % OG
plot(R2(1,:),R2(3,:)*r2d,'r');          % MILPE
title("time history of \theta_2");
xlabel('t(s)');
ylabel('\theta_2(deg)');
legend("OG","MILPE");



