
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
%                 MILPE - Lorenz 1963 Prediction (Demo)               %
%                                                                     %
%                                                                     %
%ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø%


%% Preprocess

% Controller (default)
te    =  0.001;                 % [s] simulation duration
dt    =  0.000001;              % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots
tvec  =  [0:m-1]*dt;            % time vector


%% Lorenz 1963 simulation (orig) for eigenvector extraction

% Initial Condition
x   = -8;
y   =  7;
z   =  27;

% Parameter for Lorenz system
L1  =  10;                      % sigma
L2  =  8/3;                     % beta
L3  =  28;                      % rho

% Time-loop
for it=1:m

    % Lorenz 
    u = L1*(y-x);
    v = x*(L3-z)-y;
    w = x*y-L2*z;

    % Time advancing
    x = x + u*dt;
    y = y + v*dt;
    z = z + w*dt;

    % Input snapshot (NOTE: Capital x is used as var name, used lowercase in preprint)
    X = [x y z x*z x*y]; 

    % Output snapshot (NOTE: Capital Y is used as var name, used lowercase in preprint)
    Y = [u v w];

    % Initializing unified solution space (Z)
    if ( it == 1 ) 
        Z = zeros( size(X,2)+size(Y,2) , m );
    end

    % Storing snapshots to unified solution space
    Z(:,it) = [X Y]'; % store snapshot as column vector

end


%% MILPE 

% SVD for eigenvector extraction
[U, S, V] = svd(Z,'econ');

% MILPE 
nX       =  size(X,2);              % number of variables in input subspace 
Ux       =  U(1:nX     , 1:nX);     % projection matrix on input subspace   
Uy       =  U(nX+1:end , 1:nX);     % projection matrix on output subspace  
UyUxP    =  Uy*pinv(Ux);            % MILPE low-rank governing equation (Uy*Ux+)


%% Verification 1 - Governing equation
UyUxP





