








% MILPE demo - CFD Pressure field regression and prediction 
% Mar.2026 https://www.mdpi.com/3762868

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

% MILPE
UyUxP  =  MILPE(X,Y,1);              % MILPE low-rank governing equation (Uy*Ux+)






%% MILPE prediction - Case #1 (forecasting)
% purpose: Predict pressure field with approximated governing equation
% forced-oscillation r=(0.1)*sin(0.4*pi*t)

Y_MILPE = UyUxP * X;






%% LSQ
a_LSQ = L2(X,Y);
Y_LSQ = a_LSQ * X;






%% Verification 1 - time history comparison for case #1

% controller
te    =  10;                    % [s] simulation duration
dt    =  0.04;                  % [s] time-step
m     =  floor((te)/dt+1);      % number of snapshots
t     =  [0:m-1]*dt; % prepare time vector

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
for i=1:size(Y_MILPE,1) 
    plot(t,Y_MILPE(i,:)); hold on; 
end 
title("case #1 - 3072 pressure signals (MILPE)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% 3072 pressure signals (LSQ)
figure(13) 
for i=1:size(Y_LSQ,1) 
    plot(t,Y_LSQ(i,:)); hold on; 
end 
title("case #1 - 3072 pressure signals (LSQ)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% force (spatial sum of pressures)
figure(14) 
plot(t,sum( Y,1)); hold on; % OG
plot(t,sum(Y_MILPE,1)); hold on; % MILPE
plot(t,sum(Y_LSQ,1)); hold on; % MILPE
title("case #1 - spatial sum of pressures");
xlabel("t[s]");
ylabel("\Sigma p[pa]");
legend("OG","MILPE","LSQ");





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
t     =  [0:m-1]*dt; % prepare time vector

% input subspace X
X = [v; a; a.*a; v.*v; v.*a; a.*a.*a; v.*v.*v];  % candidate #5 from Table A1

Y_MILPE = UyUxP * X;





% LSQ
Y_LSQ   = a_LSQ * X;






%% Verification 2 - time history comparison for case #2

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
for i=1:size(Y_MILPE,1) 
    plot(t,Y_MILPE(i,:)); hold on; 
end 
title("case #2 - 3072 pressure signals (MILPE)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% 3072 pressure signals (LSQ)
figure(23) 
for i=1:size(Y_LSQ,1) 
    plot(t,Y_LSQ(i,:)); hold on; 
end 
title("case #2 - 3072 pressure signals (LSQ)");
xlabel("t[s]");
ylabel("p[pa]");
ylim([-40,40]);

% force (spatial sum of pressures)
figure(24) 
plot(t,sum( Y,1)); hold on; % OG
plot(t,sum(Y_MILPE,1)); hold on; % MILPE
plot(t,sum(Y_LSQ,1)); hold on; % LSQ
title("case #2 - spatial sum of pressures");
xlabel("t[s]");
ylabel("\Sigma p[pa]");
legend("OG","MILPE","LSQ");

aaf(2,4,1)
