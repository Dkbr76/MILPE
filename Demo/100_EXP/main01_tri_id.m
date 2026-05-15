








% MILPE principle - 1
% 2 0 2 6 . 0 5 
% Reconstruct cos(t) with sin(t): 
%   - same as triagonal identity and probably related to time-delay problem?

clear all; close all;






% ctrl
te = 1e+1;                  % simulation duration
dt = 1e-2;                  % time-step
m  = floor((te)/dt+1);      % number of snapshots

% time-vector
t  = [0:m]*dt;

% input time-history(s)
x1 = sin(t);                % <- baseline
x2 = sin(t + 1e-13);        % variant

% input-subspace
X  = [x1;x2];

% output-subspace
Y  = cos(t);

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
    X = [ x1(it) x2(it) ];

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
plot(TH0(1,:),TH0(  2,:),'c'); hold on;   % OG - x1
plot(TH0(1,:),TH0(  3,:),'b--'); hold on; % OG - x2
plot(TH0(1,:),TH0(end,:),'k'); hold on;   % OG - y
plot(TH1(1,:),TH1(end,:),'m--');          % MILPE - y
xlabel('t');
ylabel('vars');

aaf




















