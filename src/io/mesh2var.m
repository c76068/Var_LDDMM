function [X,Xi] = mesh2var(pts,tri,eps)
% [X,tf,Xi] = FCATOMS(pts,f,tri,signal_type) compute the dirac representation of meshes (1D or 2D)
%
% Input :
%  pts :  matrix with vertices (one column per dimension, matrix nxd)
%  tri : connectivity matrix. if tri is of size (M x 2) this is a 1d current 
%  	(so each line of tri contains the two indexes of the points defining 
%  	an oriented segment), if tri is (M x 3) this to be considered 
%  	as a 2d current (so each line of tri contains the three indexes of 
%  	the points defining an oriented triangle).
%  eps : meshes with area <= eps will be removed, default is 1e-5
%
% Output
% X : the matrix of the centers of the faces   (M x d)
% tf: is the column vector of functional values interpolated on the X (M x 1)
% Xi: is the matric of p-vectors (tangent (1d current) or normal 2d current) (M x 2 or 3)
%
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, G. Nardi, A. Trouve (2012-2015)
    
    if nargin<3
        eps = 1e-5;
    end

    [T,M]=size(tri);
    d = size(pts,2);

    %tangent vectors
    Xi=cell(M-1,1);
    for m=1:(M-1)
        Xi{m}=pts(tri(:,m+1),:)-pts(tri(:,1),:); 
    end

    %centers
    X=sum(reshape(pts(tri(:),:)',d,T,M ),3)'/M; % center of the faces
    
    %remove faces with 0 areas
    if M-1==1
        Area = vecnorm(Xi{1},2,2);
    elseif M-1==2
        normal = cross(Xi{1},Xi{2});
        Area = vecnorm(normal,2,2)/2;
    end
    
    ind = Area>eps;
    X = X(ind,:);
    for m=1:(M-1)
        Xi{m} = Xi{m}(ind,:);
    end
    
end
