%
%
%
%
%
%
% 3BP function
%
% inputs: 
%  - simulation duration    (te)
%  - time-step              (dt) 
%  - initial condition      (ic, 6x3 array) 
%  - 3BP principals         (co, 1x4 array)
%
% output:
%  - time-histories         (THx, 7xm array)
%
%
%
%
%
function [TH1,TH2,TH3] = TBP(te,dt,ic,co) 

    m    =  floor((te)/dt+1);       % number of snapshots
    
    R1   =  ic(1,:);
    R2   =  ic(2,:);
    R3   =  ic(3,:);
    V1   =  ic(4,:);
    V2   =  ic(5,:);
    V3   =  ic(6,:);

    m1   =  co(1);                   % [kg]
    m2   =  co(2);                   % [kg]
    m3   =  co(3);                   % [kg]
    G    =  co(4);                   % [m3kg-1s-2] gravitational constant (normalized)

    % time-loop
    for it=1:m
    
        % echo
        if ( it == 1 ) cnt = 20; end % counter
        if ( mod(it,0.2*m) < 1 ) disp("3BP......"+cnt+"%"); cnt = cnt+20; end
        
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
        SNP1 = [it*dt R1 V1 A1]'; % body-1
        SNP2 = [it*dt R2 V2 A2]'; % body-2
        SNP3 = [it*dt R3 V3 A3]'; % body-3
    
        % init R0 at it==1
        if ( it == 1 ) 
            TH1 = zeros(size(SNP1,1),m); 
            TH2 = zeros(size(SNP2,1),m); 
            TH3 = zeros(size(SNP3,1),m); 
        end
    
        % save as time-history
        TH1(:,it) = SNP1; % body-1
        TH2(:,it) = SNP2; % body-2
        TH3(:,it) = SNP3; % body-3
    
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


