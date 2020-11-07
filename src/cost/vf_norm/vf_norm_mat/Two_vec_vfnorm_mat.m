function [g,PXY] = Two_vec_vfnorm_mat(X,Y,objfun)
center_faceX=X.center; center_faceY=Y.center;
[Nx,d]=size(center_faceX);
Ny=size(center_faceY,1);



%compute squared distances 
distance_center_faceXX = zeros(Nx);
distance_center_faceYY=zeros(Ny);
distance_center_faceXY=zeros(Nx,Ny);


for l=1:d
    distance_center_faceXX = distance_center_faceXX+(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1)).^2;
    distance_center_faceYY = distance_center_faceYY+(repmat(center_faceY(:,l),1,Ny)-repmat(center_faceY(:,l)',Ny,1)).^2;
    distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)).^2;
end

%compute inner products between 2-vectors


scp_normalsXX = Two_vec_Inn(X,X);
scp_normalsYY = Two_vec_Inn(Y,Y);
scp_normalsXY = Two_vec_Inn(X,Y);


%norms of 2-vectors
norm_normalsX = sqrt(diag(scp_normalsXX));
norm_normalsY = sqrt(diag(scp_normalsYY));

% Area
AreaXX = (norm_normalsX * norm_normalsX');
AreaYY = (norm_normalsY * norm_normalsY');
AreaXY = (norm_normalsX * norm_normalsY');

%normalized the inner products

scp_unit_normalsXX = scp_normalsXX./AreaXX;
scp_unit_normalsYY = scp_normalsYY./AreaYY;
scp_unit_normalsXY = scp_normalsXY./AreaXY;

% Geometric kernel
Kernel_geomXX = radial_function_geom(distance_center_faceXX,0,objfun);
Kernel_geomYY = radial_function_geom(distance_center_faceYY,0,objfun);
Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);


% tangent space kernel
Kernel_tanXX = radial_function_sphere(scp_unit_normalsXX,0,objfun);
Kernel_tanYY = radial_function_sphere(scp_unit_normalsYY,0,objfun);
Kernel_tanXY = radial_function_sphere(scp_unit_normalsXY,0,objfun);



% norm(x)=
PXX = sum(sum(AreaXX .* Kernel_geomXX .* Kernel_tanXX ));

% morm(y) =
PYY =sum(sum(AreaYY .* Kernel_geomYY .* Kernel_tanYY));

%prs(x,y) =
PXY = sum(sum(AreaXY .* Kernel_geomXY .* Kernel_tanXY));

%end

g= PXX + PYY - 2* PXY;
end
