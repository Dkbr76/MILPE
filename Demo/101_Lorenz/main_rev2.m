








% MILPE demo - Data-driven Lorenz 1963 Prediction 
% Mar.2026 https://www.mdpi.com/3762868

clear all; close all;





%% OG Lorenz simulation 1
% purpose: data aquisition for the eigenvector extraction

% controller 
te    =    1.0;                   % [s] simulation duration
dt    =    1e-6;                  % [s] time-step
ic    =    [-8 7 27];             % initial condition
co    =    [10 8/3 28];           % LRZ coeff

% Lorenz function
TH0    =    LRZ(te,dt,ic,co);






%% OG Lorenz simulation 2
% purpose: Simulate Lorenz longer time for the comparison against MILPE prediction

% controller 
te    =    10.0;                  % [s] simulation duration
dt    =    1e-6;                  % [s] time-step
ic    =    [-8 7 27];             % initial condition
co    =    [10 8/3 28];           % LRZ coeff

% Lorenz function
TH1    =    LRZ(te,dt,ic,co);






%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation

% var alloc (committed overwriting var names)
x = TH0(2,:);
y = TH0(3,:);
z = TH0(4,:);
u = TH0(5,:);
v = TH0(6,:);
w = TH0(7,:);

% input subspace
X = [x; y; z; x.*z; x.*y]; 

% output subspace
Y = [u; v; w];

% MILPE
UyUxP = MILPE(X,Y,0); % MILPE low-rank governing equation (Uy*Ux+)






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
    if ( it == 1 ) cnt = 20; end % counter
    if ( mod(it,0.2*m) < 1 ) disp("MILPE......"+cnt+"%"); cnt = cnt+20; end

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
    SNP = [it*dt x y z u v w]';

    % init R0 at it==1
    if ( it == 1 ) TH2 = zeros(size(SNP,1),m); end

    % save as time-history
    TH2(:,it) = SNP;

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
NRMSE1 = NRMSE(A,UyUxP)






%% Verification 2 - time history comparison

% time history - trajectory
figure(1)
plot3(TH1(2,:),TH1(3,:),TH1(4,:),'k'); hold on; % OG
plot3(TH2(2,:),TH2(3,:),TH2(4,:),'r');          % MILPE
xlabel('x');
ylabel('y');
zlabel('z');

% time history - x-position
figure(11)
plot(TH1(1,:),TH1(2,:),'k'); hold on; % OG
plot(TH2(1,:),TH2(2,:),'r');          % MILPE
xlabel('t');
ylabel('x');

% time history - y-position
figure(12)
plot(TH1(1,:),TH1(3,:),'k'); hold on; % OG
plot(TH2(1,:),TH2(3,:),'r');          % MILPE
xlabel('t');
ylabel('y');

% time history - z-position
figure(13)
plot(TH1(1,:),TH1(4,:),'k'); hold on; % OG
plot(TH2(1,:),TH2(4,:),'r');          % MILPE
xlabel('t');
ylabel('z');



