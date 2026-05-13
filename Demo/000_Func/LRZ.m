%
%
%
%
%
%
% Lorenz function
%
% inputs: 
%  - simulation duration    (te)
%  - time-step              (dt) 
%  - initial condition      (ic, 1x3 array) 
%  - LRZ coeff              (co, 1x3 array)
%
% output:
%  - time-history           (TH, 7xm array)
%
%
%
%
%
function TH = LRZ(te,dt,ic,co) 

    m = floor((te)/dt+1);      % number of snapshots
    
    x  = ic(1);
    y  = ic(2);
    z  = ic(3);

    L1 = co(1);
    L2 = co(2);
    L3 = co(3);

    % time-loop
    for it=1:m
    
        % echo
        if ( it == 1 ) cnt = 20; end % counter
        if ( mod(it,0.2*m) < 1 ) disp("LRZ......"+cnt+"%"); cnt = cnt+20; end
    
        % Lorenz 
        u = L1*(y-x);
        v = x*(L3-z)-y;
        w = x*y-L2*z;
    
        % time-advance
        x = x + u*dt;
        y = y + v*dt;
        z = z + w*dt;
    
        % store snapshot (7x1 array) info 
        SNP = [it*dt x y z u v w]';
    
        % init R0 at it==1
        if ( it == 1 ) TH = zeros(size(SNP,1),m); end
    
        % save as time-history
        TH(:,it) = SNP;

    end 

end % end of func


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
%
%


