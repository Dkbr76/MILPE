








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











