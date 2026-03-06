
%%
%
%
%
%
%
%                     
%           _____ ______   ___  ___       ________  _______        
%          |\   _ \  _   \|\  \|\  \     |\   __  \|\  ___ \     
%          \ \  \\\__\ \  \ \  \ \  \    \ \  \|\  \ \   __/|    
%           \ \  \\|__| \  \ \  \ \  \    \ \   ____\ \  \_|/__  
%            \ \  \    \ \  \ \  \ \  \____\ \  \___|\ \  \_|\ \ 
%             \ \__\    \ \__\ \__\ \_______\ \__\    \ \_______\
%              \|__|     \|__|\|__|\|_______|\|__|     \|_______|
%                                                      
%            MILPE GitHub-Dkbr76
%
%
%
%
%
%
%
%
%
%            Apr.2024
%
%            #MILPE #SVD #Koopman #ROM
%            #SystemIdentification 
%            #LinearRepresentation 
%
%
%
%
%

function UyUxP = MILPE(X,Y,nMe) 

    % Pre
    Z         = [X;Y];
    nX        = size(X,1);

    % MILPE
    [U, S, V] = svd(Z,'econ');
    Ux        = U(   1:nX  ,  1:nX-nMe);
    Uy        = U(nX+1:end ,  1:nX-nMe);
    UyUxP     = Uy*pinv(Ux);

end


% X     =  Input subspace   (Row: vars. Col: snapshots)
% Y     =  Output subspace  (Row: vars. Col: snapshots)
% Z     =  Unified space    (Row: vars. Col: snapshots)
% nX    =  Number of vars in input subspace        
% nMe   =  Number of modes excluded               
% Ux    =  Projection matrix on input subspace    
% Uy    =  Projection matrix on output subspace    
% UyUxP =  Approx. Governing Equation




