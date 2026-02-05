
%%
%
%
%
%
%
%
%                     
%           _____ ______   ___  ___       ________  _______        
%          |\   _ \  _   \|\  \|\  \     |\   __  \|\  ___ \     
%          \ \  \\\__\ \  \ \  \ \  \    \ \  \|\  \ \   __/|    
%           \ \  \\|__| \  \ \  \ \  \    \ \   ____\ \  \_|/__  
%            \ \  \    \ \  \ \  \ \  \____\ \  \___|\ \  \_|\ \ 
%             \ \__\    \ \__\ \__\ \_______\ \__\    \ \_______\
%              \|__|     \|__|\|__|\|_______|\|__|     \|_______|
%                                                      
%            MILPE 
%            Apr.2024
%
%
%            *Author info:    
%            Dkbr(Dkbr76@gmail.com)
%
%            *MILPE Preprint:
%            https://www.researchgate.net/profile/Dk-Br   
%           
%            *Description:
%            MILPE is an algorithm that effectively reconstructs the
%            original governing equation of data. The purpose is to seek 
%            a multivariate regression model (or multi-input--multi-output 
%            function) which has extrapolation capability, ultimately. 
%            In order to have a sufficient regression capability, 
%            the Hilbert space where the eigenvectors are extracted from 
%            should be multidimensional linear.
%
%            *Keywords:
%            #MILPE #SVD #LinearRepresentation 
%            #SystemIdentification #ROM #Koopman
%
%
%
%
%
%


%%

%ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø%
%                                                                     %
%                                                                     %
%              MILPE - Double-Pendulum Prediction (Demo)              %
%                                                                     %
%                                                                     %
%ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø%


%% Preprocess

% Controller (default)
te    =  3.0;                   % [s] simulation duration
dt    =  0.000001;              % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots
AFLAG =  1;                     % Global/Local approach flag    
                                % (1: Global, 2: Local)

% Factors
d2r   =  pi/180.0;
r2d   =  180.0/pi;


%% Double-pendulum simulation (orig) for eigenvector extraction

% Initial Condition
t1   =  10*d2r;
t2   =  10*d2r;
t1d  =  0;
t2d  =  0;

% Parameter for Double-Pendulum system
L1   =  1.0; % [m]
L2   =  1.0; % [m]
m1   =  1; % [kg]
m2   =  1; % [kg]
g    =  9.8; % [m/s2]

% Time-loop
for it=1:m

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

    % Time advancing
    t1d = t1d + t1dd*dt;
    t2d = t2d + t2dd*dt;
    t1  = t1  + t1d *dt;
    t2  = t2  + t2d *dt;
    
    % Input snapshots (Global approach)
    if ( AFLAG == 1 )
        denom = 3-cos(2*t1-2*t2);
        X11 = sin(t1)                       / denom;
        X12 = sin(t1-2*t2)                  / denom;
        X13 = sin(t1-t2)*t2d^2              / denom;
        X14 = sin(t1-t2)*t1d^2*cos(t1-t2)   / denom;
        X21 = sin(t1-t2)*t1d^2              / denom;
        X22 = sin(t1-t2)*cos(t1)            / denom;
        X23 = sin(t1-t2)*t2d^2*cos(t1-t2)   / denom;
    
        X1  = [X11 X12 X13 X14];
        X2  = [X21 X22 X23];
    end

    % Input snapshots (Local approach)
    if ( AFLAG == 2 )
        X11 = t1;
        X12 = t2;
        X13 = t1.^3;
        X14 = t2.^3;
    
        X1  = [X11 X12 X13 X14];
        X2  = [X11 X12 X13 X14];
    end

    % Output snapshots (NOTE: Capital Y is used as var name, used lowercase in preprint)
    Y1 = [t1dd];
    Y2 = [t2dd];

    % Initializing unified solution space (Z) at iter #1
    if ( it == 1 ) 
        Z1 = zeros( size(X1,2)+size(Y1,2) , m );
        Z2 = zeros( size(X2,2)+size(Y2,2) , m );
    end

    % Storing snapshots to unified solution space
    Z1(:,it) = [X1 Y1]'; % store snapshot as column vector
    Z2(:,it) = [X2 Y2]'; 

end


%% MILPE

% SVD for eigenvector extraction
[U1, S1, V1] = svd(Z1,'econ');
[U2, S2, V2] = svd(Z2,'econ');

% MILPE - 1
nX1       =  size(X1,2);               % number of variables in input subspace 
U1x       =  U1(1:nX1     , 1:nX1);    % projection matrix on input subspace   
U1y       =  U1(nX1+1:end , 1:nX1);    % projection matrix on output subspace  
U1yU1xP   =  U1y*pinv(U1x);            % MILPE low-rank governing equation     

% MILPE - 2
nX2       =  size(X2,2);             
U2x       =  U2(1:nX2     , 1:nX2);  
U2y       =  U2(nX2+1:end , 1:nX2);  
U2yU2xP   =  U2y*pinv(U2x);          


%% Verification 1 - Governing equation

U1yU1xP
U2yU2xP















