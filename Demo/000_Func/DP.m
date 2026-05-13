%
%
%
%
%
%
% Double-pendulum function
%
% inputs: 
%  - simulation duration    (te)
%  - time-step              (dt) 
%  - initial condition      (ic, 1x4 array) 
%  - DP principals          (co, 1x5 array)
%
% output:
%  - time-history           (TH, 7xm array)
%
%
%
%
%
function TH = DP(te,dt,ic,co) 

    m    =  floor((te)/dt+1);       % number of snapshots
    d2r  =  pi/180.0;               % multiplication factor for deg to rad conversion
    
    t1   =  ic(1)*d2r;              % [rad]   theta1, angle of the 1st pendulum
    t2   =  ic(2)*d2r;              % [rad]   theta2, angle of the 2nd pendulum     
    t1d  =  ic(3)*d2r;              % [rad/s] theta1_dot, angular velocity of the 1st pendulum 
    t2d  =  ic(3)*d2r;              % [rad/s] theta2_dot, angular velocity of the 2nd pendulum

    L1   =  co(1);                  % [m]     length of 1st arm
    L2   =  co(2);                  % [m]     length of 2nd arm
    m1   =  co(3);                  % [kg]    mass of 1st pendulum
    m2   =  co(4);                  % [kg]    mass of 2nd pendulum
    g    =  co(5);                  % [m/s2]  gravitational acceleration

    % time-loop
    for it=1:m
    
        % echo
        if ( it == 1 ) cnt = 20; end % counter
        if ( mod(it,0.2*m) < 1 ) disp("LRZ......"+cnt+"%"); cnt = cnt+20; end
        
        % Double-Pendulum 
        a1  =  -g*(2*m1+m2)*sin(t1);
        b1  =  -m2*g*sin(t1-2*t2);
        c1  =  -2*sin(t1-t2)*m2*(t2d^2*L2);
        d1  =  -2*sin(t1-t2)*m2*(t1d^2*L1*cos(t1-t2));
        e1  =  L1*(2*m1+m2-m2*cos(2*t1-2*t2));
    
        a2  =  2*sin(t1-t2)*t1d^2*L1*(m1+m1);
        b2  =  2*sin(t1-t2)*g*(m1+m2)*cos(t1);
        c2  =  2*sin(t1-t2)*t2d^2*L2*m2*cos(t1-t2);
        e2  =  L2*(2*m1+m2-m2*cos(2*t1-2*t2));
        
        t1dd  = (a1+b1+c1+d1)/e1;
        t2dd  = (a2+b2+c2)/e2;
    
        % Time marching (1st Euler)
        t1d = t1d + t1dd*dt;
        t2d = t2d + t2dd*dt;
        t1  = t1  + t1d *dt;
        t2  = t2  + t2d *dt;
        
        % store snapshot info
        SNP = [it*dt t1 t2 t1d t2d t1dd t2dd]';
    
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


