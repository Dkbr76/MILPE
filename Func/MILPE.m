
%%
%
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
%            MILPE 
%            Apr.2024
%
%
%            Author info:    
%            Dkbr(Dkbr76@gmail.com)
%
%            MILPE Preprint:
%            https://www.researchgate.net/profile/Dk-Br   
%           
%            Description:
%            MILPE is an algorithm that effectively reconstructs the
%            original governing equation of data. The purpose is to seek 
%            a multivariate regression model (or multi-input--multi-output 
%            function) which has extrapolation capability, ultimately. 
%            In order to have a sufficient regression capability, 
%            the Hilbert space where the eigenvectors are extracted from 
%            should be multidimensional linear.
%
%            Keywords:
%            #MILPE #SVD #LinearRepresentation 
%            #SystemIdentification #ROM #Koopman
%
%
%
%
%
%

%ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø%
%                                                                     %
%                                                                     %
%                     (ex) MILPE as a Function                        %
%                                                                     %
%                                                                     %
%ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø,¸,ø¤°º¤ø%


function UyUxP = MILPE(Z,nX,nMe) 

    [U, S, V] = svd(Z,'econ');
    Ux        = U(   1:nX  ,  1:nX-nMe);
    Uy        = U(nX+1:end ,  1:nX-nMe);
    UyUxP     = Uy*pinv(Ux);

end

% Z      =  Unified Solution space (Row: vars. Col: snapshots)
% nX     =  Number of vars in input subspace        
% nMe    =  Number of modes excluded               
% Ux     =  Projection matrix on input  subspace    
% Uy     =  Projection matrix on output subspace    
% UyUxP  =  Abb. of Governing Equation Approximated by MILPE (=Uy*Ux+)




