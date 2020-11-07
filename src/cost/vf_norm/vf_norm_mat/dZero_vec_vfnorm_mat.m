function g = dZero_vec_vfnorm_mat(X,Y,objfun)
   center_faceX=X.center; center_faceY=Y.center;
   [Tx,d]=size(center_faceX);
   Ty=size(center_faceY,1);
   
   if isfield(X,'weight') == 0
      X.weight = ones(Tx,1);
   end
   
   if isfield(Y,'weight') == 0
      Y.weight = ones(Ty,1);
   end 
   
   norm_normalsX = X.weight;
   norm_normalsY = Y.weight;
   
   distance2_center_faceXX = zeros(Tx);
   distance2_center_faceXY=zeros(Tx,Ty);
   
   for l=1:d
       distance2_center_faceXX = distance2_center_faceXX+(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1)).^2;
       distance2_center_faceXY = distance2_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;

   end
   
%    Kernel_geomXX  = radial_function_geom(distance2_center_faceXX,0,objfun);
%    Kernel_geomXY  = radial_function_geom(distance2_center_faceXY,0,objfun);
   dKernel_geomXX = radial_function_geom(distance2_center_faceXX,1,objfun);
   dKernel_geomXY = radial_function_geom(distance2_center_faceXY,1,objfun);         
   
   AreaXX = (norm_normalsX * norm_normalsX');
   AreaXY = (norm_normalsX * norm_normalsY');

   DXX  = AreaXX .* dKernel_geomXX;
   DXY  = AreaXY .* dKernel_geomXY;
        
   Dcenter_faceXg=zeros(Tx,d);   
        
   for l=1:d
      Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1))),2)...
                -sum(DXY .* (repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)),2)); % scalar kernel case
   end
   
   g.center = Dcenter_faceXg;
end