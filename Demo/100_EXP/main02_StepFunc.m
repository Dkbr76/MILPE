








% MILPE principle - 2
% 2 0 2 6 . 0 5 
% Reconstruct step-function with cos(t)

clear all; close all;






% ctrl
te = 2*pi;                  % simulation duration
dt = 1e-5;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots

% time-vector
t  = [0:m]*dt;

% input time-history(s)
x1  = cos(t* 1);              
x2  = cos(t* 3);       
x3  = cos(t* 5);              
x4  = cos(t* 7);              
x5  = cos(t* 9);              
x6  = cos(t*11);              

% input time-history(s)
y1  = sign(cos(t));

% input-subspace
X  =  [ x1; x2; x3; x4; x5; x6 ];

% output-subspace
Y  =  [ y1 ];

% MILPE
UyUxP = MILPE(X,Y,0)

% sav time-history
TH0 = [t;X;Y];






% MILPE-Prediction

% time-loop
for it=1:m

    % echo
    if ( it == 1 ) cnt = 20; end % counter
    if ( mod(it,0.2*m) < 1 ) disp("MILPE......"+cnt+"%"); cnt = cnt+20; end

    % input snapshot X
    X = [ x1(it)  x2(it)  x3(it)  x4(it)  x5(it)  x6(it) ];

    % output snapshot Y
    Y = UyUxP * X';

    % store snapshot info
    SNP = [it*dt X Y]';

    % init TH at it==1
    if ( it == 1 ) TH1 = zeros(size(SNP,1),m); end

    % save as time-history
    TH1(:,it) = SNP;

end
   



% verification - compare with Analytic DFT result

% pure analytic
Analytic    =  [ 4/pi -4/(3*pi) 4/(5*pi) -4/(7*pi) 4/(9*pi) -4/(11*pi) ]
D           =  Analytic - UyUxP

% DFT
[amp,phs,f] =  DFT(y1,dt);  % exec DFT for y1, frequency might dependent on circumstances
y1_RECON    =  DFT_RECON(amp,phs,f,dt,m,12);





% time history 
figure(11)
plot(TH0(1,:),TH0(  2,:),'c'); hold on; % OG - x1
plot(TH0(1,:),TH0(  3,:),'c'); hold on; % OG - x2
plot(TH0(1,:),TH0(  4,:),'c'); hold on; % OG - x3
plot(TH0(1,:),TH0(  5,:),'c'); hold on; % OG - x4
plot(TH0(1,:),TH0(  6,:),'c'); hold on; % OG - x5
plot(TH0(1,:),TH0(  7,:),'c'); hold on; % OG - x6
plot(TH0(1,:),TH0(end,:),'k'); hold on; % OG - y1
plot(TH1(1,:),TH1(end,:),'k--');        % MILPE - y1
plot(TH1(1,:),y1_RECON,'m--');            % DFT - y1
xlabel('t');
ylabel('vars');

aaf
