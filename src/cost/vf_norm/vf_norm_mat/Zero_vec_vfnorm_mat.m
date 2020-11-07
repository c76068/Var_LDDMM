function [g,PXY] = Zero_vec_vfnorm_mat(X,Y,objfun)

   center_faceX=X.center; center_faceY=Y.center;
   [Nx,d]=size(center_faceX);
   Ny=size(center_faceY,1);
   
   if isfield(X,'weight') == 0
      X.weight = ones(Nx,1);
   end
   
   if isfield(Y,'weight') == 0
      Y.weight = ones(Ny,1);
   end 
   
   norm_normalsX = X.weight;
   norm_normalsY = Y.weight;
   
   %compute squared distances 
   distance_center_faceXX = zeros(Nx);
   distance_center_faceYY=zeros(Ny);
   distance_center_faceXY=zeros(Nx,Ny);

 
   for l=1:d
      distance_center_faceXX = distance_center_faceXX+(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1)).^2;
      distance_center_faceYY = distance_center_faceYY+(repmat(center_faceY(:,l),1,Ny)-repmat(center_faceY(:,l)',Ny,1)).^2;
      distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)).^2;   
   end
   
   % Geometric kernel
   Kernel_geomXX = radial_function_geom(distance_center_faceXX,0,objfun);
   Kernel_geomYY = radial_function_geom(distance_center_faceYY,0,objfun);
   Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);
   
   % Area
   AreaXX = (norm_normalsX * norm_normalsX');
   AreaYY = (norm_normalsY * norm_normalsY');
   AreaXY = (norm_normalsX * norm_normalsY');

   % norm(x)=
   PXX = sum(sum(AreaXX .* Kernel_geomXX));

   % morm(y) =
   PYY =sum(sum(AreaYY .* Kernel_geomYY));

   %prs(x,y) =
   PXY = sum(sum(AreaXY .* Kernel_geomXY));

   %end

   g= PXX + PYY - 2* PXY;
end