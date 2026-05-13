








% MILPE principle - 1
% 2 0 2 6 . 0 5 
% Reconstruct cos(2*pi*0.5*t) with sin(2*pi*0.5*t): probably time-delay problem?

clear all; close all;

% ctrl
te = 10;                    % simulation duration
dt = 1e-2;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots

% time-vector
t  = [0:m]*dt;

% input time-history(s)
x1 = sin( t*2*pi*0.5 + pi*0.0);
x2 = sin( t*2*pi*0.5 + pi*0.1);
x3 = sin( t*2*pi*0.5 + pi*0.2);
x4 = sin( t*2*pi*0.5 + pi*0.3);
x5 = sin( t*2*pi*0.5 + pi*0.4);

% input-subspace
X  =  [x1;x2;x3;x4;x5];

% output-subspace
Y  = cos(t*2*pi*0.5);

% MILPE
Z         =  [X;Y];         % unified space
nX        =  size(X,1);     % dim of X
nMe       =  0;             % number of modes excluded
[U,S,V]   =  svd(Z,'econ'); U
Ux        =  U(   1:nX  ,  1:nX-nMe);
Uy        =  U(nX+1:end ,  1:nX-nMe);
UyUxP     =  Uy*pinv(Ux)

% sav time-history
TH0 = [t;X;Y];





% MILPE-Prediction

% time-loop
for it=1:m

    % echo
    if ( it == 1 ) cnt = 20; end % counter
    if ( mod(it,0.2*m) < 1 ) disp("MILPE......"+cnt+"%"); cnt = cnt+20; end

    % input snapshot X
    X = [ x1(it) x2(it) x3(it) x4(it) x5(it)];

    % output snapshot Y
    Y = UyUxP * X';

    % store snapshot info
    SNP = [it*dt X Y]';

    % init TH at it==1
    if ( it == 1 ) TH1 = zeros(size(SNP,1),m); end

    % save as time-history
    TH1(:,it) = SNP;

end
   
% time history - x-position
figure(11)
plot(TH0(1,:),TH0(  2,:),'c'); hold on; % OG - x1
plot(TH0(1,:),TH0(  3,:),'c'); hold on; % OG - x2
plot(TH0(1,:),TH0(  4,:),'c'); hold on; % OG - x3
plot(TH0(1,:),TH0(  5,:),'c'); hold on; % OG - x4
plot(TH0(1,:),TH0(  6,:),'c'); hold on; % OG - x5
plot(TH0(1,:),TH0(end,:),'k'); hold on; % OG - y
plot(TH1(1,:),TH1(end,:),'m--');        % MILPE - y
xlabel('t');
ylabel('vars');

aaf
