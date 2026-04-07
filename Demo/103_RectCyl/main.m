








% MILPE demo - CFD Pressure field regression and prediction 
% Mar.2026 https://www.mdpi.com/3762868
% YouTube: 
https://www.youtube.com/watch?v=gvJTZSLJ0lc
https://www.youtube.com/watch?v=Oc1LzhAe6wk

clear all; close all;





%% Data load
% purpose: load CFD time histories for eigenvector extraction and
% comparision with prediction
load("case1/v.mat"); % vel
load("case1/a.mat"); % acc
load("case1/Y.mat"); % pres





%% MILPE algorithm
% purpose: extract eigenvectors and construct approximated governing equation

% input subspace X
X = [v; a; a.*a; v.*v; v.*a; a.*a.*a; v.*v.*v];  % candidate #5 from Table A1

% unified solution space Z
Z = [X;Y];

% SVD for eigenvector extraction
[U, S, V] = svd(Z,'econ');

% governing equation
nXE       =  1;                        % number of mode excluded
nX        =  size(X,1);                % number of variables in input subspace 
Ux        =  U(1:nX     , 1:nX-nXE);   % projection matrix on input subspace   
Uy        =  U(nX+1:end , 1:nX-nXE);   % projection matrix on output subspace  
UyUxP     =  Uy*pinv(Ux);              % MILPE low-rank governing equation (Uy*Ux+)





%% MILPE prediction - Case #1 (forecasting)
% purpose: Predict pressure field with approximated governing equation
% forced-oscillation r=(0.1)*sin(0.4*pi*t)

% controller
te    =  10;                    % [s] simulation duration
dt    =  0.04;                  % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots

% time-loop
for it=1:m

    % echo
    if ( mod(it,100) == 0 ) disp(it);  end

    % input snapshot X
    X = [v(it) a(it) a(it)^2 v(it)^2 v(it)*a(it) a(it)^3 v(it)^3]; % candidate #5 from Table A1

    % output snapshot Y
    Ytmp = UyUxP * X';      

    % init Ym at it==1 (Ym: Y_MILPE)
    if ( it == 1 ) Ym = zeros(size(Ytmp,1),m); end

    % save as time-history
    Ym(:,it) = Ytmp;

end





%% Verification 1 - time history comparison for case #1

t = [0:m-1]*dt; % prepare time vector

% 3072 pressure signals (OG)
figure(11) 
for i=1:size( Y,1) 
    plot(t,Y(i,:)); hold on; 
end 
title("case #1 - 3072 pressure signals (OG)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% 3072 pressure signals (MILPE)
figure(12) 
for i=1:size(Ym,1) 
    plot(t,Ym(i,:)); hold on; 
end 
title("case #1 - 3072 pressure signals (MILPE)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% force (spatial sum of pressures)
figure(13) 
plot(t,sum( Y,1)); hold on; % OG
plot(t,sum(Ym,1)); hold on; % MILPE
title("case #1 - spatial sum of pressures");
xlabel("t[s]");
ylabel("\Sigma p[pa]");
legend("OG","MILPE");





%% MILPE prediction - Case #2 (prediction for untrained scenario)
% forced-oscillation r=(0.1)*sin(0.4*pi*t)+(0.05)*sin(0.2*pi*t)

% Data load
load("case2/v.mat"); % vel
load("case2/a.mat"); % acc
load("case2/Y.mat"); % pres

% controller
te    =  20;                    % [s] simulation duration
dt    =  0.04;                  % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots

% time-loop
for it=1:m

    % echo
    if ( mod(it,100) == 0 ) disp(it);  end

    % input snapshot X
    X = [v(it) a(it) a(it)^2 v(it)^2 v(it)*a(it) a(it)^3 v(it)^3]; % candidate #5 from Table A1

    % output snapshot Y
    Ytmp = UyUxP * X';      % use UyUxP constructed from case #1

    % init Ym at it==1 (Ym: Y_MILPE)
    if ( it == 1 ) Ym = zeros(size(Ytmp,1),m); end

    % save as time-history
    Ym(:,it) = Ytmp;

end





%% Verification 2 - time history comparison for case #2

t = [0:m-1]*dt; % prepare time vector

% 3072 pressure signals (OG)
figure(21) 
for i=1:size( Y,1) 
    plot(t,Y(i,:)); hold on; 
end 
title("case #2 - 3072 pressure signals (OG)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% 3072 pressure signals (MILPE)
figure(22) 
for i=1:size(Ym,1) 
    plot(t,Ym(i,:)); hold on; 
end 
title("case #2 - 3072 pressure signals (MILPE)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% force (spatial sum of pressures)
figure(23) 
plot(t,sum( Y,1)); hold on; % OG
plot(t,sum(Ym,1)); hold on; % MILPE
title("case #2 - spatial sum of pressures");
xlabel("t[s]");
ylabel("\Sigma p[pa]");
legend("OG","MILPE");
