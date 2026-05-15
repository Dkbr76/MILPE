%
%
%
%
%
%
%
%
%
function [amp,phs,f] = DFT(X,dt)

    % DFT
    m        =   size(X,2);
    y        =   fft(X,m,2);                            % dft
    y        =   y(:, 1:floor(m/2));                    % use half - round down
    amp      =   abs(y) / floor(m/2);                   % [mm] amp
    amp(:,1) =   amp(:,1) / 2;                          % halve a0 <- will replace a0 with mean
    phs      =   angle(y);                              % [rad] phs    
    f        =   (1/dt) * (0:floor(m/2)-1) / m;         % [1/s] f-vec

end
%
%
%
%
%
%
%
%
%
%
