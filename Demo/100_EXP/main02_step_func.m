








% MILPE principle - 2
% 2 0 2 6 . 0 5 
% Reconstruct step-function with cos(t)

clear all; close all;






% ctrl
te = 1e+1;                  % simulation duration
dt = 1e-2;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots

% time-vector
t  = [0:m]*dt;

% input time-history(s)
x1  = cos(t);              
x2  = cos(t* 2);       
x3  = cos(t* 3);              
x4  = cos(t* 4);              
x5  = cos(t* 5);              
x6  = cos(t* 6);              
x7  = cos(t* 7);              
x8  = cos(t* 8);              
x9  = cos(t* 9);              
x10 = cos(t*10);              
x11 = cos(t*11);
x12 = cos(t*12);
x13 = cos(t*13);
x14 = cos(t*14);
x15 = cos(t*15);
x16 = cos(t*16);
x17 = cos(t*17);
x18 = cos(t*18);
x19 = cos(t*19);
x20 = cos(t*20);
x21 = cos(t*21);
x22 = cos(t*22);
x23 = cos(t*23);
x24 = cos(t*24);
x25 = cos(t*25);
x26 = cos(t*26);
x27 = cos(t*27);
x28 = cos(t*28);
x29 = cos(t*29);
x30 = cos(t*30);

% input-subspace
X  =  [ x1; x2; x3; x4; x5; x6; x7; x8; x9;x10; ...
       x11;x12;x13;x14;x15;x16;x17;x18;x19;x20; ...
       x21;x22;x23;x24;x25;x26;x27;x28;x29;x30];

% output-subspace
Y  =  2*double(cos(t) >= 0)-1;

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
    X = [  x1(it)  x2(it)  x3(it)  x4(it)  x5(it)  x6(it)  x7(it)  x8(it)  x9(it) x10(it) ...
          x11(it) x12(it) x13(it) x14(it) x15(it) x16(it) x17(it) x18(it) x19(it) x20(it) ...
          x21(it) x22(it) x23(it) x24(it) x25(it) x26(it) x27(it) x28(it) x29(it) x30(it) ];

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
plot(TH0(1,:),TH0(  7,:),'c'); hold on; % OG - x6
plot(TH0(1,:),TH0(  8,:),'c'); hold on; % OG - x7
plot(TH0(1,:),TH0(  9,:),'c'); hold on; % OG - x8
plot(TH0(1,:),TH0( 10,:),'c'); hold on; % OG - x9
plot(TH0(1,:),TH0(end,:),'k'); hold on; % OG - y
plot(TH1(1,:),TH1(end,:),'m--');        % MILPE - y
xlabel('t');
ylabel('vars');

aaf
