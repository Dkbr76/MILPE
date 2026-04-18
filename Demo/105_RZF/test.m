








% MILPE demo - Data-driven Riemann Zeta Function
% Apr.2026

clear all; close all;

FOGS = 0;       % switch for OG sim
FMLP = 1;       % switch for MILPE sim
SRID = 10;      % id-sav for OG sim
LRID = 10;      % id-load  
LPID = 100;      % id-post


%% OG RZF simulation 
% purpose: (1) data aquisition  
%          (2) eigenvector extraction

if( FOGS == 1 )

% controller
te    =    3;                     % [s] simulation duration
dt    =    1e-6;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots
t     =    [0:m-1]*dt;

% Riemann Zeta Function
s     =    0.5 + i*t;
z     =    zeta(s);
x     =    real(z); 
y     =    imag(z);

% derivatives
for i=1:m-1
    u(i) = (x(i+1)-x(i))/dt;
    v(i) = (y(i+1)-y(i))/dt;
end
u(m) = u(m-1);
v(m) = v(m-1);

% sav time-history
TH1 = [t;x;y;u;v];

% sav init cond (for later use)
x0    =    x(1);
y0    =    y(1);

% sav mat
save(SRID+".mat","TH1");

end % switch FOGS



%% MILPE pre
% purpose: extract eigenvectors and construct approximated governing equation

if ( FMLP == 1 )

% load
load(LRID+".mat");

% var alloc
x  = TH1(2,:);
y  = TH1(3,:);
u  = TH1(4,:);
v  = TH1(5,:);

% input subspace X1
X1 = [x; y; x.*y; x.*x.*y; x.*x.*x];

% input subspace X2
X2 = [x; y; x.*y; x.*x.*y; x.*x.*x];

% output subspace Y1
Y1 = [u];

% output subspace Y2
Y2 = [v];

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
U1yU1xP  =  U1y*pinv(U1x)              % MILPE low-rank governing equation (Uy*Ux+)

nX2      =  size(X2,1);                
U2x      =  U2(1:nX2     , 1:nX2);     
U2y      =  U2(nX2+1:end , 1:nX2);     
U2yU2xP  =  U2y*pinv(U2x)     

%% MILPE prediction
% purpose: Predict RZF with approximated governing equation

% controller
te    =    15;                    % [s] simulation duration
dt    =    1e-4;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots
t     =    [0:m-1]*dt;

% init cond
x     =    TH1(2,1);
y     =    TH1(3,1);

% time-loop
for it=1:m

    % echo
    if ( mod(it,100000) == 0 ) disp(it);  end

    % input snapshot X1
    X1 = [x y x*y x*x*y x^3];
    
    % input snapshot X2
    X2 = [x y x*y x*x*y x^3];

    % output snapshot Y1
    u = U1yU1xP * X1';

    % output snapshot Y2
    v = U2yU2xP * X2';

    % time-advance
    x = x + u*dt;
    y = y + v*dt;

    % store snapshot info
    R = [it*dt x y u v]';

    % init R0 at it==1
    if ( it == 1 ) TH2 = zeros(size(R,1),m); end

    % save as time-history
    TH2(:,it) = R;

end % time-loop


%% Verification 2 - time history comparison

% load TH1
load(LPID+".mat");

% time history - trajectory
figure(1)
plot(TH1(2,:),TH1(3,:),'k'); hold on; % OG
plot(TH2(2,:),TH2(3,:),'r');          % MILPE

figure(2); 
subplot(2,2,1); plot(TH1(1,:),TH1(2,:),'k'); hold on; plot(TH2(1,:),TH2(2,:),'r--'); title("x"); legend("OG","MILPE"); xlim([0,6]);
subplot(2,2,2); plot(TH1(1,:),TH1(3,:),'k'); hold on; plot(TH2(1,:),TH2(3,:),'r--'); title("y"); legend("OG","MILPE"); xlim([0,6]);
subplot(2,2,3); plot(TH1(1,:),TH1(4,:),'k'); hold on; plot(TH2(1,:),TH2(4,:),'r--'); title("u"); legend("OG","MILPE"); xlim([0,6]);
subplot(2,2,4); plot(TH1(1,:),TH1(5,:),'k'); hold on; plot(TH2(1,:),TH2(5,:),'r--'); title("v"); legend("OG","MILPE"); xlim([0,6]);


end % FMLP

AAF(2,3,3)


