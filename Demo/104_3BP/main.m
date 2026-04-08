








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

% controller 
te    =    1;                     % [s] simulation duration
dt    =    1e-5;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition - disp 
R1    =    [ 0.97 -0.243 0.1];   % gave a little tweak to OG IC from ref  
R2    =    [-0.97  0.243 -0.1];  % to generate 3D orbit
R3    =    [0 0 0.05];

% initial condition - vel 
V1    =    [ 0.466 0.432  0];    
V2    =    [ 0.466 0.432  0]; 
V3    =    [ -0.932 -0.864 0];

% parameters 
m1    =    1;                   % [kg]
m2    =    1;                   % [kg]
m3    =    1;                   % [kg]
G     =    1;                   % [m3kg-1s-2] gravitational constant (normalized)

% time-loop
for it=1:m

    % echo
    if ( mod(it,1000) == 0 ) disp(it);  end

    % 3BP 
    R12 = norm(R1-R2);
    R13 = norm(R1-R3);
    R23 = norm(R2-R3);
    A1  = G*m2*(R2-R1)/(R12.^3) + G*m3*(R3-R1)/(R13.^3);
    A2  = G*m1*(R1-R2)/(R12.^3) + G*m3*(R3-R2)/(R23.^3);
    A3  = G*m1*(R1-R3)/(R13.^3) + G*m2*(R2-R3)/(R23.^3);

    % time-advance - acc 2 vel
    V1 = V1 + A1.*dt;
    V2 = V2 + A2.*dt;
    V3 = V3 + A3.*dt;

    % time-advance - vel 2 disp
    R1 = R1 + V1.*dt;
    R2 = R2 + V2.*dt;
    R3 = R3 + V3.*dt;

    % store snapshot info
    T1 = [it*dt R1 V1 A1]'; % body-1
    T2 = [it*dt R2 V2 A2]'; % body-2
    T3 = [it*dt R3 V3 A3]'; % body-3

    % init R0 at it==1
    if ( it == 1 ) 
        TH01 = zeros(size(T1,1),m); 
        TH02 = zeros(size(T2,1),m); 
        TH03 = zeros(size(T3,1),m); 
    end

    % save as time-history
    TH01(:,it) = T1; % body-1
    TH02(:,it) = T2; % body-2
    TH03(:,it) = T3; % body-3

end





%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation

% var alloc (committed overwriting var names)
x1    = TH01( 2,:); % body-1
y1    = TH01( 3,:);
z1    = TH01( 4,:);
u1    = TH01( 5,:); 
v1    = TH01( 6,:);
w1    = TH01( 7,:);
udot1 = TH01( 8,:); 
vdot1 = TH01( 9,:);
wdot1 = TH01(10,:);

x2    = TH02( 2,:); % body-2
y2    = TH02( 3,:);
z2    = TH02( 4,:);
u2    = TH02( 5,:); 
v2    = TH02( 6,:);
w2    = TH02( 7,:);
udot2 = TH02( 8,:); 
vdot2 = TH02( 9,:);
wdot2 = TH02(10,:);

x3    = TH03( 2,:); % body-3
y3    = TH03( 3,:);
z3    = TH03( 4,:);
u3    = TH03( 5,:); 
v3    = TH03( 6,:);
w3    = TH03( 7,:);
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
Zx1 = [Xx1;Yx1];
Zy1 = [Xy1;Yy1];
Zz1 = [Xz1;Yz1];

Zx2 = [Xx2;Yx2];
Zy2 = [Xy2;Yy2];
Zz2 = [Xz2;Yz2];

Zx3 = [Xx3;Yx3];
Zy3 = [Xy3;Yy3];
Zz3 = [Xz3;Yz3];

for i=1:9

    if(i==1) X = Xx1; Z = Zx1; end
    if(i==2) X = Xy1; Z = Zy1; end
    if(i==3) X = Xz1; Z = Zz1; end
    if(i==4) X = Xx2; Z = Zx2; end
    if(i==5) X = Xy2; Z = Zy2; end
    if(i==6) X = Xz2; Z = Zz2; end
    if(i==7) X = Xx3; Z = Zx3; end
    if(i==8) X = Xy3; Z = Zy3; end
    if(i==9) X = Xz3; Z = Zz3; end

    % SVD for eigenvector extraction
    [U, S, V] = svd(Z,'econ');
    
    % governing equation
    nX       =  size(X,1);              % number of variables in input subspace 
    Ux       =  U(1:nX     , 1:nX);     % projection matrix on input subspace   
    Uy       =  U(nX+1:end , 1:nX);     % projection matrix on output subspace  
    UyUxP    =  Uy*pinv(Ux);            % MILPE low-rank governing equation (Uy*Ux+)

    if(i==1) UyUxPx1 = UyUxP; end
    if(i==2) UyUxPy1 = UyUxP; end
    if(i==3) UyUxPz1 = UyUxP; end
    if(i==4) UyUxPx2 = UyUxP; end
    if(i==5) UyUxPy2 = UyUxP; end
    if(i==6) UyUxPz2 = UyUxP; end
    if(i==7) UyUxPx3 = UyUxP; end
    if(i==8) UyUxPy3 = UyUxP; end
    if(i==9) UyUxPz3 = UyUxP; end

end





%% OG 3BP simulation 2
% purpose: Simulate Lorenz longer time for the comparison against MILPE prediction

% controller 
te    =    50.0;                  % [s] simulation duration
dt    =    1e-5;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition - disp 
R1    =    [ 0.97 -0.243 0.1];   % gave a little tweak to OG IC from ref  
R2    =    [-0.97  0.243 -0.1];  % to generate 3D orbit
R3    =    [0 0 0.05];

% initial condition - vel 
V1    =    [ 0.466 0.432  0];    
V2    =    [ 0.466 0.432  0]; 
V3    =    [ -0.932 -0.864 0];

% parameters 
m1    =    1;                   % [kg]
m2    =    1;                   % [kg]
m3    =    1;                   % [kg]
G     =    1;                   % [m3kg-1s-2] gravitational constant (normalized)

% time-loop
for it=1:m

    % echo
    if ( mod(it,1000) == 0 ) disp(it);  end

    % 3BP 
    R12 = norm(R1-R2);
    R13 = norm(R1-R3);
    R23 = norm(R2-R3);
    A1  = G*m2*(R2-R1)/(R12.^3) + G*m3*(R3-R1)/(R13.^3);
    A2  = G*m1*(R1-R2)/(R12.^3) + G*m3*(R3-R2)/(R23.^3);
    A3  = G*m1*(R1-R3)/(R13.^3) + G*m2*(R2-R3)/(R23.^3);

    % time-advance - acc 2 vel
    V1 = V1 + A1.*dt;
    V2 = V2 + A2.*dt;
    V3 = V3 + A3.*dt;

    % time-advance - vel 2 disp
    R1 = R1 + V1.*dt;
    R2 = R2 + V2.*dt;
    R3 = R3 + V3.*dt;

    % store snapshot info
    T1 = [it*dt R1 V1 A1]'; % body-1
    T2 = [it*dt R2 V2 A2]'; % body-2
    T3 = [it*dt R3 V3 A3]'; % body-3

    % init R0 at it==1
    if ( it == 1 ) 
        TH11 = zeros(size(T1,1),m); 
        TH12 = zeros(size(T2,1),m); 
        TH13 = zeros(size(T3,1),m); 
    end

    % save as time-history
    TH11(:,it) = T1; % body-1
    TH12(:,it) = T2; % body-2
    TH13(:,it) = T3; % body-3

end





%% MILPE prediction
% purpose: Predict Lorenz with approximated governing equation

% controller 
te    =    50.0;                  % [s] simulation duration
dt    =    1e-5;                  % [s] time-step
m     =    floor((te)/dt+1);      % number of snapshots

% initial condition - disp 
R1    =    [ 0.97 -0.243 0.1];   % gave a little tweak to OG IC from ref  
R2    =    [-0.97  0.243 -0.1];  % to generate 3D orbit
R3    =    [0 0 0.05];

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
    if ( mod(it,10000) == 0 ) disp(it);  end

    % input subspace(s)
    R1  = [x1 y1 z1]; % just to use sqrt(dot) at below
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
    T1 = [it*dt x1 y1 z1 u1 v1 w1 udot1 vdot1 wdot1]'; % body-1
    T2 = [it*dt x2 y2 z2 u2 v2 w2 udot2 vdot2 wdot2]'; % body-2
    T3 = [it*dt x3 y3 z3 u3 v3 w3 udot3 vdot3 wdot3]'; % body-3

    % init R0 at it==1
    if ( it == 1 ) 
        TH21 = zeros(size(T1,1),m); 
        TH22 = zeros(size(T2,1),m); 
        TH23 = zeros(size(T3,1),m); 
    end

    % save as time-history
    TH21(:,it) = T1; % body-1
    TH22(:,it) = T2; % body-2
    TH23(:,it) = T3; % body-3

end





%% Verification 1 - governing equation and NRMSE

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
UyUxP1 = [UyUxPx1;UyUxPy1;UyUxPz1]
UyUxP2 = [UyUxPx2;UyUxPy2;UyUxPz2]
UyUxP3 = [UyUxPx3;UyUxPy3;UyUxPz3]

% governing equation difference
D1 = UyUxP1 - A1;
D2 = UyUxP2 - A2;
D3 = UyUxP3 - A3;

% NRMSE
NUMER = sqrt( sum(sum(D1.^2,2)) / ( size(D1,1)*size(D1,2) ) ); %numerator
DENOM = sqrt( sum(sum(A1.^2,2)) / ( size(A1,1)*size(A1,2) ) ); %denomenator
NRMSE1 = NUMER/DENOM

NUMER = sqrt( sum(sum(D2.^2,2)) / ( size(D2,1)*size(D2,2) ) ); %numerator
DENOM = sqrt( sum(sum(A2.^2,2)) / ( size(A2,1)*size(A2,2) ) ); %denomenator
NRMSE2 = NUMER/DENOM

NUMER = sqrt( sum(sum(D3.^2,2)) / ( size(D3,1)*size(D3,2) ) ); %numerator
DENOM = sqrt( sum(sum(A3.^2,2)) / ( size(A3,1)*size(A3,2) ) ); %denomenator
NRMSE3 = NUMER/DENOM




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

AAF(4,3,1); 
