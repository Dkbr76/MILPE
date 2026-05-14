








% MILPE demo - 3 body problem (3BP)
% Mar.2026 ref: https://www.mdpi.com/3762868

clear all; close all;





%% Nomenclature
% displacement:     R = (x,y,z)
% velocity:         V = (u,v,w)
% acceleration:     A = (udot,vdot,wdot)
% subscript:        body number





%% OG 3BP simulation 1
% purpose: data aquisition for the eigenvector extraction
% IC reference: https://arxiv.org/pdf/math/0011268
disp("OG 3BP Sim #1");
te    =    1e+0;                      % [s] simulation duration
dt    =    1e-5;                      % [s] time-step
R1    =    [  0.97  -0.243   0.1];    % ic - disp of body 1 : gave a little tweak to OG IC from ref to generate 3D orbit
R2    =    [ -0.97   0.243  -0.1];    % ic - disp of body 2 : gave a little tweak to OG IC from ref to generate 3D orbit
R3    =    [  0      0       0.05];   % ic - disp of body 3 : gave a little tweak to OG IC from ref to generate 3D orbit
V1    =    [  0.466  0.432   0];      % ic - vel of body 1
V2    =    [  0.466  0.432   0];      % ic - vel of body 2
V3    =    [ -0.932 -0.864   0];      % ic - vel of body 3
ic    =    [ R1; R2; R3; V1; V2; V3]; % ic mat input
co    =    [ 1 1 1 1];                % 3BP param input - m1, m2, m3, G
[TH01,TH02,TH03] = TBP(te,dt,ic,co);  % run 3BP sim and obtain time-histories of body1-3


%% OG 3BP simulation 2
% purpose: Simulate 3BP longer time for the comparison against MILPE prediction
disp("OG 3BP Sim #2");
te    =    5e+1;                      % [s] simulation duration
dt    =    1e-5;                      % [s] time-step
R1    =    [  0.97  -0.243   0.1];    % ic - disp of body 1 : gave a little tweak to OG IC from ref to generate 3D orbit
R2    =    [ -0.97   0.243  -0.1];    % ic - disp of body 2 : gave a little tweak to OG IC from ref to generate 3D orbit
R3    =    [  0      0       0.05];   % ic - disp of body 3 : gave a little tweak to OG IC from ref to generate 3D orbit
V1    =    [  0.466  0.432   0];      % ic - vel of body 1
V2    =    [  0.466  0.432   0];      % ic - vel of body 2
V3    =    [ -0.932 -0.864   0];      % ic - vel of body 3
ic    =    [ R1; R2; R3; V1; V2; V3]; % ic mat input
co    =    [ 1 1 1 1];                % 3BP param input - m1, m2, m3, G
[TH11,TH12,TH13] = TBP(te,dt,ic,co);  % run 3BP sim and obtain time-histories of body1-3




%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation

% var alloc (committed overwriting var names)
x1    = TH01( 2,:); % body-1
y1    = TH01( 3,:);
z1    = TH01( 4,:);
udot1 = TH01( 8,:); 
vdot1 = TH01( 9,:);
wdot1 = TH01(10,:);

x2    = TH02( 2,:); % body-2
y2    = TH02( 3,:);
z2    = TH02( 4,:);
udot2 = TH02( 8,:); 
vdot2 = TH02( 9,:);
wdot2 = TH02(10,:);

x3    = TH03( 2,:); % body-3
y3    = TH03( 3,:);
z3    = TH03( 4,:);
udot3 = TH03( 8,:); 
vdot3 = TH03( 9,:);
wdot3 = TH03(10,:);

% input subspace(s)
R1  = [x1; y1; z1];
R2  = [x2; y2; z2];
R3  = [x3; y3; z3];
R12 = sqrt(dot(R1-R2,R1-R2));
R13 = sqrt(dot(R1-R3,R1-R3));
R23 = sqrt(dot(R2-R3,R2-R3));

Xx1 = [ x2./(R12.^3); x1./(R12.^3); x3./(R13.^3); x1./(R13.^3) ]; 
Xy1 = [ y2./(R12.^3); y1./(R12.^3); y3./(R13.^3); y1./(R13.^3) ];
Xz1 = [ z2./(R12.^3); z1./(R12.^3); z3./(R13.^3); z1./(R13.^3) ];

Xx2 = [ x1./(R12.^3); x2./(R12.^3); x3./(R23.^3); x2./(R23.^3) ];
Xy2 = [ y1./(R12.^3); y2./(R12.^3); y3./(R23.^3); y2./(R23.^3) ];
Xz2 = [ z1./(R12.^3); z2./(R12.^3); z3./(R23.^3); z2./(R23.^3) ];

Xx3 = [ x1./(R13.^3); x3./(R13.^3); x2./(R23.^3); x3./(R23.^3) ]; 
Xy3 = [ y1./(R13.^3); y3./(R13.^3); y2./(R23.^3); y3./(R23.^3) ];
Xz3 = [ z1./(R13.^3); z3./(R13.^3); z2./(R23.^3); z3./(R23.^3) ];

% output subspace(s)
Yx1 = [udot1];
Yy1 = [vdot1];
Yz1 = [wdot1];

Yx2 = [udot2];
Yy2 = [vdot2];
Yz2 = [wdot2];

Yx3 = [udot3];
Yy3 = [vdot3];
Yz3 = [wdot3];

% unified solution space(s) (Z)
UyUxPx1 = MILPE(Xx1,Yx1,0);
UyUxPy1 = MILPE(Xy1,Yy1,0);
UyUxPz1 = MILPE(Xz1,Yz1,0);

UyUxPx2 = MILPE(Xx2,Yx2,0);
UyUxPy2 = MILPE(Xy2,Yy2,0);
UyUxPz2 = MILPE(Xz2,Yz2,0);

UyUxPx3 = MILPE(Xx3,Yx3,0);
UyUxPy3 = MILPE(Xy3,Yy3,0);
UyUxPz3 = MILPE(Xz3,Yz3,0);







%% MILPE prediction
% purpose: Predict 3BP with approximated governing equation

% controller 
te    =    5e+1;                  % [s] simulation duration
dt    =    1e-5;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition - disp 
R1    =    [ 0.97 -0.243  0.1];   % gave a little tweak to OG IC from ref  
R2    =    [-0.97  0.243 -0.1];   % to generate 3D orbit
R3    =    [ 0     0      0.05];

x1  = R1(1);
y1  = R1(2);
z1  = R1(3);
x2  = R2(1);
y2  = R2(2);
z2  = R2(3);
x3  = R3(1);
y3  = R3(2);
z3  = R3(3);

% initial condition - vel 
V1    =    [ 0.466 0.432  0];    
V2    =    [ 0.466 0.432  0]; 
V3    =    [ -0.932 -0.864 0];

u1  = V1(1);
v1  = V1(2);
w1  = V1(3);
u2  = V2(1);
v2  = V2(2);
w2  = V2(3);
u3  = V3(1);
v3  = V3(2);
w3  = V3(3);

% time-loop
for it=1:m

    % echo
    if ( it == 1 ) cnt = 20; end % counter
    if ( mod(it,0.2*m) < 1 ) disp("MILPE......"+cnt+"%"); cnt = cnt+20; end

    % input subspace(s)
    R1  = [x1 y1 z1]; % just to utilize sqrt(dot) at below
    R2  = [x2 y2 z2];
    R3  = [x3 y3 z3];
    R12 = sqrt(dot(R1-R2,R1-R2));
    R13 = sqrt(dot(R1-R3,R1-R3));
    R23 = sqrt(dot(R2-R3,R2-R3));
    
    Xx1 = [ x2/(R12^3) x1/(R12^3) x3/(R13^3) x1/(R13^3) ]; 
    Xy1 = [ y2/(R12^3) y1/(R12^3) y3/(R13^3) y1/(R13^3) ];
    Xz1 = [ z2/(R12^3) z1/(R12^3) z3/(R13^3) z1/(R13^3) ];
    
    Xx2 = [ x1/(R12^3) x2/(R12^3) x3/(R23^3) x2/(R23^3) ]; 
    Xy2 = [ y1/(R12^3) y2/(R12^3) y3/(R23^3) y2/(R23^3) ];
    Xz2 = [ z1/(R12^3) z2/(R12^3) z3/(R23^3) z2/(R23^3) ];
    
    Xx3 = [ x1/(R13^3) x3/(R13^3) x2/(R23^3) x3/(R23^3) ]; 
    Xy3 = [ y1/(R13^3) y3/(R13^3) y2/(R23^3) y3/(R23^3) ];
    Xz3 = [ z1/(R13^3) z3/(R13^3) z2/(R23^3) z3/(R23^3) ];

    % output snapshot Y
    udot1 = UyUxPx1 * Xx1';
    vdot1 = UyUxPy1 * Xy1';
    wdot1 = UyUxPz1 * Xz1';
    
    udot2 = UyUxPx2 * Xx2';
    vdot2 = UyUxPy2 * Xy2';
    wdot2 = UyUxPz2 * Xz2';

    udot3 = UyUxPx3 * Xx3';
    vdot3 = UyUxPy3 * Xy3';
    wdot3 = UyUxPz3 * Xz3';

    % time-advance - vel
    u1 = u1 + udot1*dt;
    v1 = v1 + vdot1*dt;
    w1 = w1 + wdot1*dt;

    u2 = u2 + udot2*dt;
    v2 = v2 + vdot2*dt;
    w2 = w2 + wdot2*dt;

    u3 = u3 + udot3*dt;
    v3 = v3 + vdot3*dt;
    w3 = w3 + wdot3*dt;

    % time-advance - disp
    x1 = x1 + u1*dt;
    y1 = y1 + v1*dt;
    z1 = z1 + w1*dt;

    x2 = x2 + u2*dt;
    y2 = y2 + v2*dt;
    z2 = z2 + w2*dt;

    x3 = x3 + u3*dt;
    y3 = y3 + v3*dt;
    z3 = z3 + w3*dt;

    % store snapshot info
    SNP1 = [it*dt x1 y1 z1 u1 v1 w1 udot1 vdot1 wdot1]'; % body-1
    SNP2 = [it*dt x2 y2 z2 u2 v2 w2 udot2 vdot2 wdot2]'; % body-2
    SNP3 = [it*dt x3 y3 z3 u3 v3 w3 udot3 vdot3 wdot3]'; % body-3

    % init R0 at it==1
    if ( it == 1 ) 
        TH21 = zeros(size(SNP1,1),m); 
        TH22 = zeros(size(SNP2,1),m); 
        TH23 = zeros(size(SNP3,1),m); 
    end

    % save as time-history
    TH21(:,it) = SNP1; % body-1
    TH22(:,it) = SNP2; % body-2
    TH23(:,it) = SNP3; % body-3

end





%% Verification 1 - governing equation and NRMSE

m1 = co(1);
m2 = co(2);
m3 = co(3);
G  = co(end);

% OG governing equation
A1 = [ G*m2 -G*m2 G*m3 -G*m3;
       G*m2 -G*m2 G*m3 -G*m3;
       G*m2 -G*m2 G*m3 -G*m3];

A2 = [ G*m1 -G*m1 G*m3 -G*m3;
       G*m1 -G*m1 G*m3 -G*m3;
       G*m1 -G*m1 G*m3 -G*m3];

A3 = [ G*m1 -G*m1 G*m2 -G*m2;
       G*m1 -G*m1 G*m2 -G*m2;
       G*m1 -G*m1 G*m2 -G*m2];

% MILPE approximated governing equation
UyUxP1 = [UyUxPx1; UyUxPy1; UyUxPz1]
UyUxP2 = [UyUxPx2; UyUxPy2; UyUxPz2]
UyUxP3 = [UyUxPx3; UyUxPy3; UyUxPz3]

% governing equation difference
D1 = UyUxP1 - A1
D2 = UyUxP2 - A2
D3 = UyUxP3 - A3

% NRMSE
NRMSE1 = NRMSE(A1,UyUxP1)
NRMSE2 = NRMSE(A2,UyUxP2)
NRMSE3 = NRMSE(A3,UyUxP3)



%% Verification 2 - time history comparison

% time history - trajectory

figure(1) % body-1
plot3(TH11(2,:),TH11(3,:),TH11(4,:),'k'); hold on; % OG
plot3(TH21(2,:),TH21(3,:),TH21(4,:),'r');          % MILPE
xlabel('x');
ylabel('y');
zlabel('z');
legend("OG","MILPE");
title("Body1 trajectory");

figure(2) % body-2
plot3(TH12(2,:),TH12(3,:),TH12(4,:),'k'); hold on; % OG
plot3(TH22(2,:),TH22(3,:),TH22(4,:),'r');          % MILPE
xlabel('x');
ylabel('y');
zlabel('z');
legend("OG","MILPE");
title("Body2 trajectory");

figure(3) % body-3
plot3(TH13(2,:),TH13(3,:),TH13(4,:),'k'); hold on; % OG
plot3(TH23(2,:),TH23(3,:),TH23(4,:),'r');          % MILPE
xlabel('x');
ylabel('y');
zlabel('z');
legend("OG","MILPE");
title("Body3 trajectory");

figure(11); plot(TH11(1,:),TH11(2,:),'k'); hold on; plot(TH21(1,:),TH21(2,:),'r'); xlabel('t(s)'); ylabel('x'); legend("OG","MILPE"); title("Body1 x"); % body-1 x
figure(12); plot(TH11(1,:),TH11(3,:),'k'); hold on; plot(TH21(1,:),TH21(3,:),'r'); xlabel('t(s)'); ylabel('y'); legend("OG","MILPE"); title("Body1 y"); % body-1 y
figure(13); plot(TH11(1,:),TH11(4,:),'k'); hold on; plot(TH21(1,:),TH21(4,:),'r'); xlabel('t(s)'); ylabel('z'); legend("OG","MILPE"); title("Body1 z"); % body-1 z

figure(21); plot(TH12(1,:),TH12(2,:),'k'); hold on; plot(TH22(1,:),TH22(2,:),'r'); xlabel('t(s)'); ylabel('x'); legend("OG","MILPE"); title("Body2 x"); % body-2 x
figure(22); plot(TH12(1,:),TH12(3,:),'k'); hold on; plot(TH22(1,:),TH22(3,:),'r'); xlabel('t(s)'); ylabel('y'); legend("OG","MILPE"); title("Body2 y"); % body-2 y
figure(23); plot(TH12(1,:),TH12(4,:),'k'); hold on; plot(TH22(1,:),TH22(4,:),'r'); xlabel('t(s)'); ylabel('z'); legend("OG","MILPE"); title("Body2 z"); % body-2 z

figure(31); plot(TH13(1,:),TH13(2,:),'k'); hold on; plot(TH23(1,:),TH23(2,:),'r'); xlabel('t(s)'); ylabel('x'); legend("OG","MILPE"); title("Body3 x"); % body-3 x
figure(32); plot(TH13(1,:),TH13(3,:),'k'); hold on; plot(TH23(1,:),TH23(3,:),'r'); xlabel('t(s)'); ylabel('y'); legend("OG","MILPE"); title("Body3 y"); % body-3 y
figure(33); plot(TH13(1,:),TH13(4,:),'k'); hold on; plot(TH23(1,:),TH23(4,:),'r'); xlabel('t(s)'); ylabel('z'); legend("OG","MILPE"); title("Body3 z"); % body-3 z

aaf(4,3,1); 
