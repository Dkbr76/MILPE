%
%
%
%
%
%
%
%
%
function a = L2(X,Y)

    a = Y*X'*pinv(X*X');

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
