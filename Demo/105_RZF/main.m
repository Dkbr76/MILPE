








% MILPE demo - Data-driven Riemann Zeta Function
% Apr.2026

clear all; close all;

% switches
FOGS = 0;       % switch for OG sim
FMLP = 1;       % switch for MILPE sim

% ids
SRID = 104;     % id-sav for OG sim
LRID = 104;     % id-load for MILPE sim 
LPID = 104;     % id-post for OG and MILPE comparison


%% OG RZF simulation 
% purpose: (1) data aquisition  
%          (2) eigenvector extraction

if( FOGS == 1 )

% controller
te    =    10;                       % [s] simulation duration
dt    =    1e-4;                    % [s] time-step
m     =    floor((te)/dt+1);        % number of snapshots
t     =    [0:m-1]*dt;

% Riemann Zeta Function
s     =    0.5 + i*(t+12);          % truncated t<12
z     =    zeros(1,m);
x     =    zeros(1,m);
y     =    zeros(1,m);
for it=1:m % running m times of zeta per each s is faster than putting m size s at once
    if ( mod(it,100) == 0 ) disp(it); end % echo
    z(it) =  zeta(s(it));
    x(it) =  real(z(it)); 
    y(it) =  imag(z(it));
end

% sav time-history
TH1 = [t;x;y];

% sav mat
save(SRID+".mat","TH1");

end % switch FOGS



%% MILPE pre
% purpose: extract eigenvectors and construct approximated governing equation

if ( FMLP == 1 )

% load TH1
load(LRID+".mat");

% var alloc
x    = TH1(2,:);
y    = TH1(3,:);

% #snap and dt
m    = size(x,2);
dt   = TH1(1,2)-TH1(1,1);

% 1st derivatives
for i=2:m
    u(i) = (x(i)-x(i-1))/dt;
    v(i) = (y(i)-y(i-1))/dt;
end
u(1) = u(2);
v(1) = v(2);

% 2nd derivatives
for i=2:m
    udot(i) = (u(i)-u(i-1))/dt;
    vdot(i) = (v(i)-v(i-1))/dt;
end
udot(1) = udot(2);
vdot(1) = vdot(2);

% attach back to ref TH1
TH1 = [TH1;u;v;udot;vdot];

% input subspace X (candidates)
X1 = [u; v; u.*v; u.^2.*v; u.*v.^2; u.^3; v.^3];
X2 = [u; v; u.*v; u.^2.*v; u.*v.^2; u.^3; v.^3];

% X1 = [u;v;u.*v];
% X2 = [u;v;u.*v];

% output subspace Y
Y1 = [udot];
Y2 = [vdot];

% MILPE
U1yU1xP = MILPE(X1,Y1,0);
U2yU2xP = MILPE(X2,Y2,0);



%% MILPE prediction
% purpose: Predict RZF with approximated governing equation

% controller
te    =    10.0;                    % [s] simulation duration
dt    =    1e-4;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots
t     =    [0:m-1]*dt;

% init cond
x     =    TH1(2,5);
y     =    TH1(3,5);
u     =    TH1(4,5);
v     =    TH1(5,5);
% udot  =    TH1(6,5);
% vdot  =    TH1(7,5);

% time-loop
for it=1:m

    % echo
    if ( mod(it,100000) == 0 ) disp(it);  end

    % input snapshot X (candidates)
    X1 = [u v u*v u^2*v u*v^2 u^3 v^3];
    X2 = [u v u*v u^2*v u*v^2 u^3 v^3];
%     X1 = [u v u*v];
%     X2 = [u v u*v];

    % output snapshot Y
    udot = U1yU1xP * X1';
    vdot = U2yU2xP * X2';
    
    % time-advance
    u = u + udot*dt;
    v = v + vdot*dt;

    % time-advance
    x = x + u*dt;
    y = y + v*dt;

    % store snapshot info
    R = [it*dt x y u v udot vdot]';

    % init R0 at it==1
    if ( it == 1 ) TH2 = zeros(size(R,1),m); end

    % save as time-history
    TH2(:,it) = R;

end % time-loop


%% Verification 2 - time history comparison

% % load TH1
% load(LPID+".mat");

% time history - trajectory
figure(1)
plot(TH1(2,:),TH1(3,:),'k'); hold on; % OG
plot(TH2(2,:),TH2(3,:),'r--');          % MILPE
title("trajectory"); xlabel("x"); ylabel("y");
legend("OG","MILPE");

figure(2); 
subplot(3,2,1); plot(TH1(1,:),TH1(2,:),'k'); hold on; plot(TH2(1,:),TH2(2,:),'r--'); title("x");    legend("OG","MILPE"); %xlim([0,3]);
subplot(3,2,2); plot(TH1(1,:),TH1(3,:),'k'); hold on; plot(TH2(1,:),TH2(3,:),'r--'); title("y");    legend("OG","MILPE"); %xlim([0,3]);
subplot(3,2,3); plot(TH1(1,:),TH1(4,:),'k'); hold on; plot(TH2(1,:),TH2(4,:),'r--'); title("u");    legend("OG","MILPE"); %xlim([0,3]);
subplot(3,2,4); plot(TH1(1,:),TH1(5,:),'k'); hold on; plot(TH2(1,:),TH2(5,:),'r--'); title("v");    legend("OG","MILPE"); %xlim([0,3]);
subplot(3,2,5); plot(TH1(1,:),TH1(6,:),'k'); hold on; plot(TH2(1,:),TH2(6,:),'r--'); title("udot"); legend("OG","MILPE"); %xlim([0,3]);
subplot(3,2,6); plot(TH1(1,:),TH1(7,:),'k'); hold on; plot(TH2(1,:),TH2(7,:),'r--'); title("vdot"); legend("OG","MILPE"); %xlim([0,3]);


end % FMLP

aaf(2,3,1)


