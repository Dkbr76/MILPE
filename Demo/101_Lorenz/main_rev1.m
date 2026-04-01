








% MILPE demo - Lorenz 1963 Prediction 
% Apr.2026

clear all; close all;






%% OG Lorenz simulation 
% purpose: data aquisition for the eigenvector extraction


% controller 
te    =    1.0;                   % [s] simulation duration
dt    =    1e-6;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition
x     =    -8; 
y     =     7; 
z     =    27;

% parameters 
L1    =    10;                    % sigma
L2    =   8/3;                    % beta
L3    =    28;                    % rho

% time-loop
for it=1:m

    % echo
    if ( mod(it,1000) == 0 ) disp(it);  end

    % Lorenz 
    u = L1*(y-x);
    v = x*(L3-z)-y;
    w = x*y-L2*z;

    % time-advance
    x = x + u*dt;
    y = y + v*dt;
    z = z + w*dt;

    % store snapshot info
    R = [it*dt x y z u v w];

    % init R0 at it==1
    if ( it == 1 ) R0 = zeros(size(R,2),m); end

    % save as time-history
    R0(:,it) = R';

end





%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation


% var alloc (committed overwriting var names)
x = R0(2,:);
y = R0(3,:);
z = R0(4,:);
u = R0(5,:);
v = R0(6,:);
w = R0(7,:);

% input subspace
X = [x; y; z; x.*z; x.*y]; 

% output subspace
Y = [u; v; w];

% unified solution space (Z)
Z = [X;Y];

% SVD for eigenvector extraction
[U, S, V] = svd(Z,'econ');

% governing equation
nX       =  size(X,1);              % number of variables in input subspace 
Ux       =  U(1:nX     , 1:nX);     % projection matrix on input subspace   
Uy       =  U(nX+1:end , 1:nX);     % projection matrix on output subspace  
UyUxP    =  Uy*pinv(Ux);            % MILPE low-rank governing equation (Uy*Ux+)





%% OG Lorenz simulation 
% purpose: Simulate Lorenz longer time for the comparison against MILPE prediction


% controller 
te    =    10;                    % [s] simulation duration
dt    =    1e-6;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition
x     =    -8; 
y     =     7; 
z     =    27;

% parameters 
L1    =    10;                    % sigma
L2    =   8/3;                    % beta
L3    =    28;                    % rho

% time-loop
for it=1:m

    % echo
    if ( mod(it,10000) == 0 ) disp(it);  end

    % Lorenz 
    u = L1*(y-x);
    v = x*(L3-z)-y;
    w = x*y-L2*z;

    % time-advance
    x = x + u*dt;
    y = y + v*dt;
    z = z + w*dt;

    % store snapshot info
    R = [it*dt x y z u v w];

    % init R0 at it==1
    if ( it == 1 ) R1 = zeros(size(R,2),m); end

    % save as time-history
    R1(:,it) = R';

end





%% MILPE prediction
% purpose: Predict Lorenz with approximated governing equation

% controller 
te    =    10;                    % [s] simulation duration
dt    =    1e-6;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition
x     =    -8; 
y     =     7; 
z     =    27;

% time-loop
for it=1:m

    % echo
    if ( mod(it,10000) == 0 ) disp(it);  end

    % input snapshot X
    X = [x y z x*z x*y];

    % output snapshot Y
    Y = UyUxP * X';
    
    % var alloc
    u = Y(1); 
    v = Y(2);
    w = Y(3);

    % time-advance
    x = x + u*dt;
    y = y + v*dt;
    z = z + w*dt;

    % store snapshot info
    R = [it*dt x y z u v w];

    % init R0 at it==1
    if ( it == 1 ) R2 = zeros(size(R,2),m); end

    % save as time-history
    R2(:,it) = R';

end





%% Verification 1 - governing equation and NRMSE


% OG governing equation
A = [ -10  10   0   0   0; 
       28  -1   0  -1   0;
        0   0 -8/3  0   1 ]

% MILPE approximated governing equation
UyUxP

% governing equation difference
D = UyUxP - A

% NRMSE
NUMER = sqrt( sum(sum(D.^2,2)) / ( size(D,1)*size(D,2) ) ); %numerator
DENOM = sqrt( sum(sum(A.^2,2)) / ( size(A,1)*size(A,2) ) ); %denomenator
NRMSE = NUMER/DENOM





%% Verification 2 - time history comparison


% time history - trajectory
figure(1)
plot3(R1(2,:),R1(3,:),R1(4,:),'k'); hold on; % OG
plot3(R2(2,:),R2(3,:),R2(4,:),'r');          % MILPE
xlabel('x');
ylabel('y');
zlabel('z');

% time history - x-position
figure(11)
plot(R1(1,:),R1(2,:),'k'); hold on; % OG
plot(R2(1,:),R2(2,:),'r');          % MILPE
xlabel('t');
ylabel('x');

% time history - y-position
figure(12)
plot(R1(1,:),R1(3,:),'k'); hold on; % OG
plot(R2(1,:),R2(3,:),'r');          % MILPE
xlabel('t');
ylabel('y');

% time history - z-position
figure(13)
plot(R1(1,:),R1(4,:),'k'); hold on; % OG
plot(R2(1,:),R2(4,:),'r');          % MILPE
xlabel('t');
ylabel('z');



