%
%
%
%
%
%
%
%
%
function R = NRMSE(a,b)

    % NRMSE
    D     = b - a;
    NUMER = sqrt( sum(sum(D.^2,2)) / ( size(D,1)*size(D,2) ) ); %numerator
    DENOM = sqrt( sum(sum(a.^2,2)) / ( size(a,1)*size(a,2) ) ); %denomenator
    R     = NUMER/DENOM;

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
