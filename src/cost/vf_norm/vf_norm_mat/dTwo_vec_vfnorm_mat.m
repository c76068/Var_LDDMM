function g = dTwo_vec_vfnorm_mat(X,Y,objfun)
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
AreaXY = (norm_normalsX * norm_normalsY');

%normalized the inner products

scp_unit_normalsXX = scp_normalsXX./AreaXX;
scp_unit_normalsXY = scp_normalsXY./AreaXY;

% Geometric kernel
Kernel_geomXX = radial_function_geom(distance_center_faceXX,0,objfun);
Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);
dKernel_geomXX = radial_function_geom(distance_center_faceXX,1,objfun);
dKernel_geomXY = radial_function_geom(distance_center_faceXY,1,objfun);

% tangent space kernel
Kernel_tanXX = radial_function_sphere(scp_unit_normalsXX,0,objfun);
Kernel_tanXY = radial_function_sphere(scp_unit_normalsXY,0,objfun);
dKernel_tanXX = radial_function_sphere(scp_unit_normalsXX,1,objfun);
dKernel_tanXY = radial_function_sphere(scp_unit_normalsXY,1,objfun);


%Inner products between X.vecotr and {X.vector Y.vector}
Inn = [X.vector{1};X.vector{2}]*[X.vector{1}' X.vector{2}' Y.vector{1}' Y.vector{2}'];

%quotient of norms
Q_tanXX = (1./norm_normalsX)*norm_normalsX';
Q_tanYX = (1./norm_normalsX)*norm_normalsY';


%project X.vector{i} to X.vector{j}^{\perp} times norm(X.vector{j})^2
P_1 = repmat(diag(Inn(Nx+1:2*Nx,Nx+1:2*Nx)),1,d).*X.vector{1}...
    - repmat(diag(Inn(1:Nx,Nx+1:2*Nx)),1,d).*X.vector{2};
P_2 = repmat(diag(Inn(1:Nx,1:Nx)),1,d).*X.vector{2}...
    - repmat(diag(Inn(1:Nx,Nx+1:2*Nx)),1,d).*X.vector{1};


%----------------------------------
% Compute derivative wrt 2-vectors 
%----------------------------------
A = Kernel_geomXX.*dKernel_tanXX; 
B = Kernel_geomXY.*dKernel_tanXY;
C = Kernel_geomXX.*(Kernel_tanXX - dKernel_tanXX.*scp_unit_normalsXX).*Q_tanXX;
D = Kernel_geomXY.*(Kernel_tanXY - dKernel_tanXY.*scp_unit_normalsXY).*Q_tanYX;
E = repmat(sum(C,2)-sum(D,2),1,d);

g.vector{1} = 2*((A.*Inn(Nx+1:2*Nx,Nx+1:2*Nx))*X.vector{1}-(A.*Inn(Nx+1:2*Nx,1:Nx))*X.vector{2}...
               -( (B.*Inn(Nx+1:2*Nx,2*Nx+Ny+1:2*Nx+2*Ny))*Y.vector{1} - (B.*Inn(Nx+1:2*Nx,2*Nx+1:2*Nx+Ny))*Y.vector{2} )...
               + E.*P_1);
g.vector{2} = 2*((A.*Inn(1:Nx,1:Nx))*X.vector{2}-(A.*Inn(1:Nx,Nx+1:2*Nx))*X.vector{1}...
               -( (B.*Inn(1:Nx,2*Nx+1:2*Nx+Ny))*Y.vector{2} - (B.*Inn(1:Nx,2*Nx+Ny+1:2*Nx+2*Ny))*Y.vector{1} )...
               + E.*P_2);
    

%------------------------------------
% Compute derivative wrt center_face 
%------------------------------------

DXX  = AreaXX .* dKernel_geomXX .* Kernel_tanXX;
        DXY  = AreaXY .* dKernel_geomXY .* Kernel_tanXY;
        
        Dcenter_faceXg=zeros(Nx,d);   
        for l=1:d
            Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1))),2)...
                -sum(DXY .* (repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)),2)); % scalar kernel case
        end
        
g.center =  Dcenter_faceXg;       
end