
% function [Dcenter_faceXg,varargout]= dvarifoldnorm(center_faceX,normalsX, weightX, Target,objfun,defo)
function g = dOne_vec_vfnorm_mat(X,Y,objfun)


center_faceX = X.center; center_faceY = Y.center;
normalsX = X.vector{1}; normalsY = Y.vector{1};


% Compute unit normals
        [Tx,d]=size(center_faceX);
        Ty=size(center_faceY,1);
        % Compute unit normals
        norm_normalsX = sqrt(sum(normalsX .^2,2));
        norm_normalsY = sqrt(sum(normalsY .^2,2));
        
        unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
        unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));
        
% compute squared distances and angles
 
        distance2_center_faceXX = zeros(Tx);
        distance2_center_faceXY=zeros(Tx,Ty);
               
        scp_unit_normalsXX = zeros(Tx);
        scp_unit_normalsXY = zeros(Tx,Ty);
    
        for l=1:d
            distance2_center_faceXX = distance2_center_faceXX+(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1)).^2;
            distance2_center_faceXY = distance2_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;

            scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Tx).*repmat(unit_normalsX(:,l)',Tx,1));
            scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ty).*repmat(unit_normalsY(:,l)',Tx,1));
        end

% compute  Geometric kernel      
        Kernel_geomXX  = radial_function_geom(distance2_center_faceXX,0,objfun);
        Kernel_geomXY  = radial_function_geom(distance2_center_faceXY,0,objfun);
        dKernel_geomXX = radial_function_geom(distance2_center_faceXX,1,objfun);
        dKernel_geomXY = radial_function_geom(distance2_center_faceXY,1,objfun);                
        
% compute tangent space kernel
        Kernel_tanXX  = radial_function_sphere(scp_unit_normalsXX,0,objfun);
        Kernel_tanXY  = radial_function_sphere(scp_unit_normalsXY,0,objfun);
        dKernel_tanXX = radial_function_sphere(scp_unit_normalsXX,1,objfun);
        dKernel_tanXY = radial_function_sphere(scp_unit_normalsXY,1,objfun);

% compute Area 
        AreaXX = (norm_normalsX * norm_normalsX');
        AreaXY = (norm_normalsX * norm_normalsY');

%------------------------------------
% Compute derivative wrt center_face 
%------------------------------------
        DXX  = AreaXX .* dKernel_geomXX .* Kernel_tanXX;
        DXY  = AreaXY .* dKernel_geomXY .* Kernel_tanXY;
        
        Dcenter_faceXg=zeros(Tx,d);   
        for l=1:d
            Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1))),2)...
                -sum(DXY .* (repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)),2)); % scalar kernel case
        end
        
%--------------------------------
% Compute derivative wrt normals 
%--------------------------------
        MXX = Kernel_geomXX .* dKernel_tanXX.*AreaXX;
        MXY = Kernel_geomXY .* dKernel_tanXY.*AreaXY;                    
        mXX = Kernel_geomXX .* Kernel_tanXX;
        mXY = Kernel_geomXY .* Kernel_tanXY;                
         
%         DnormalsXg = 2*(  repmat(mXX * norm_normalsX  - (MXX .* scp_unit_normalsXX  ) * norm_normalsX,1,d) .* unit_normalsX ...
%                 	+ MXX * normalsX ...
%                 	- repmat(mXY * norm_normalsY  - (MXY .* scp_unit_normalsXY  ) * norm_normalsY,1,d) .* unit_normalsX...
%                 	- MXY * normalsY);
          DnormalsXg = 2*(MXX*unit_normalsX-MXY*unit_normalsY);
          
          
          DweightX = 2*(mXX*norm_normalsX-mXY*norm_normalsY);
          
          g.center = Dcenter_faceXg;
          g.vector{1} = PJ(DnormalsXg,unit_normalsX)./repmat(norm_normalsX,1,d)+ unit_normalsX.*repmat(DweightX,1,d); 
          
          
        

%end
        
         
end

function IN = Inn(S,T)
[~,d] = size(S);
IN = repmat(sum(S.*T,2),1,d);
end

function PJ = PJ(S,T)
PJ=S-Inn(S,T).*T;
end