%
%
%
%
%
%
%
%
%
function y = DFT_RECON(amp,phs,f,dt,m,order)

    w = 2*pi*f;         % wave frequency
    n = size(amp,1);    % num of var
    y = zeros(n,m);     % init

    % Recon
    for j=1:n % var loop

        y(j,:) = y(j,:) + amp(j,1);

        for i=1:m % time loop

            for k=2:order
                y(j,i) = y(j,i) + amp(j,k)*cos(w(j,k)*i*dt+phs(j,k));
            end
        
        end
    end

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
