function [g,PXY] = One_vec_vfnorm_mat(X,Y,objfun)

center_faceX=X.center; center_faceY=Y.center;
normalsX = X.vector{1}; normalsY = Y.vector{1};

[Nx,d]=size(center_faceX);
Ny=size(center_faceY,1);
% Compute unit normals
norm_normalsX = sqrt(sum(normalsX .^2,2));
norm_normalsY = sqrt(sum(normalsY .^2,2));
% norm_normalsX = weightX;
% norm_normalsY = weightY;

unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));

% unit_normalsX = normalsX;
% unit_normalsY = normalsY;
%compute squared distances and angles
distance_center_faceXX = zeros(Nx);
distance_center_faceYY=zeros(Ny);
distance_center_faceXY=zeros(Nx,Ny);

scp_unit_normalsXX = zeros(Nx);
scp_unit_normalsYY = zeros(Ny);
scp_unit_normalsXY = zeros(Nx,Ny);

for l=1:d
    distance_center_faceXX = distance_center_faceXX+(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1)).^2;
    distance_center_faceYY = distance_center_faceYY+(repmat(center_faceY(:,l),1,Ny)-repmat(center_faceY(:,l)',Ny,1)).^2;
    distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)).^2;
    
    scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Nx).*repmat(unit_normalsX(:,l)',Nx,1));
    scp_unit_normalsYY = scp_unit_normalsYY + (repmat(unit_normalsY(:,l),1,Ny).*repmat(unit_normalsY(:,l)',Ny,1));
    scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ny).*repmat(unit_normalsY(:,l)',Nx,1));
end

% Geometric kernel
Kernel_geomXX = radial_function_geom(distance_center_faceXX,0,objfun);
Kernel_geomYY = radial_function_geom(distance_center_faceYY,0,objfun);
Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);


% tangent space kernel
Kernel_tanXX = radial_function_sphere(scp_unit_normalsXX,0,objfun);
Kernel_tanYY =  radial_function_sphere(scp_unit_normalsYY,0,objfun);
Kernel_tanXY = radial_function_sphere(scp_unit_normalsXY,0,objfun);

% Area
AreaXX = (norm_normalsX * norm_normalsX');
AreaYY = (norm_normalsY * norm_normalsY');
AreaXY = (norm_normalsX * norm_normalsY');

% norm(x)=
PXX = sum(sum(AreaXX .* Kernel_geomXX .* Kernel_tanXX ));

% morm(y) =
PYY =sum(sum(AreaYY .* Kernel_geomYY .* Kernel_tanYY));

%prs(x,y) =
PXY = sum(sum(AreaXY .* Kernel_geomXY .* Kernel_tanXY));

%end

g= PXX + PYY - 2* PXY;
end
