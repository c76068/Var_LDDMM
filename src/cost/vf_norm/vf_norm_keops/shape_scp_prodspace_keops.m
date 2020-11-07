function res=shape_scp_prodspace_keops(X,Y,objfun)
% function res=shape_scp_prodspace_keops(cX,cY,XiX,XiY,kernel_geom,kernel_grass,kernel_size_geom,kernel_size_grass)
%compute the inner product beween two discrete varifolds X and Y
cX = X.center'; 
cY = Y.center'; 
[d,Nx]=size(cX);
[~,Ny]=size(cY);
kernel_geom = objfun.kernel_geom;
kernel_size_geom = objfun.kernel_size_geom;

if isfield(X,'vector') == 0
    m = 0;     
else 
    
    XiX = X.vector;
    XiY = Y.vector;
    m=length(XiX);
    kernel_grass = objfun.kernel_grass;
    kernel_size_grass = objfun.kernel_size_grass;
end


 switch kernel_geom
           case 'gaussian'
               expr_geom=['Exp(-SqDist(x,y)*a)'];
            
           case 'cauchy'
               expr_geom=['IntCst(1)/(IntCst(1)+SqDist(x,y)*a)'];
            
 end

 
switch m

    case 0
      
        if isfield(X,'weight') == 0
           X.weight = ones(Nx,1);
        end
   
        if isfield(Y,'weight') == 0
           Y.weight = ones(Ny,1);
        end
        
        norm_XiX = X.weight;
        norm_XiY = Y.weight;
        
        K = Kernel([expr_geom '*r'],'a=Pm(0,1)',['x=Vx(1,' num2str(eval('d')) ')'],['y=Vy(2,' num2str(eval('d')) ')'],'r=Vy(3,1)');
    
        res = sum(K(1/kernel_size_geom^2,cX,cY,norm_XiY').*norm_XiX');

    case 1
       %tic
       norm_XiX=sqrt(sum(XiX{1}.^2,2));
       norm_XiY=sqrt(sum(XiY{1}.^2,2));
    
       XitildeX=XiX{1}./repmat(norm_XiX,1,d);
       XitildeY=XiY{1}./repmat(norm_XiY,1,d);
       %toc
%     
%        switch kernel_geom
%            case 'gaussian'
%                expr_geom=['Exp(-SqDist(x,y)*a)'];
%             
%            case 'cauchy'
%                expr_geom=['IntCst(1)/(IntCst(1)+SqDist(x,y)*a)'];
%             
%        end
    
       switch kernel_grass
           case 'linear'
               expr_grass=['(u|v)'];
            
           case 'gaussian_oriented'
               expr_grass=['Exp(IntCst(2)*b*((u|v)-IntCst(1)))'];
            
           case 'binet'
               expr_grass=['Square((u|v))'];
            
           case 'gaussian_unoriented'
               expr_grass=['Exp(IntCst(2)*b*(Square((u|v))-IntCst(1)))'];
            
       end
    
       K = Kernel([expr_geom '*' expr_grass '*r'],'a=Pm(0,1)','b=Pm(1,1)',['x=Vx(2,' num2str(eval('d')) ')'],['y=Vy(3,' num2str(eval('d')) ')'],['u=Vx(4,' num2str(eval('d')) ')'],['v=Vy(5,' num2str(eval('d')) ')'],'r=Vy(6,1)');
    
       res = sum(K(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XitildeX',XitildeY',norm_XiY').*norm_XiX');

    
    case 2
    
       norm_XiX=pvec_Norm(XiX);
       norm_XiY=pvec_Norm(XiY);
    
       XitildeX=XiX;
       XitildeY=XiY;
    
       XitildeX{1}=XiX{1}./repmat(norm_XiX,1,d);
       XitildeY{1}=XiY{1}./repmat(norm_XiY,1,d);
       %toc
    
%        switch kernel_geom
%            case 'gaussian'
%                expr_geom=['Exp(-SqDist(x,y)*a)'];
%             
%            case 'cauchy'
%                expr_geom=['IntCst(1)/(IntCst(1)+SqDist(x,y)*a)'];
%             
%        end
    
       switch kernel_grass
           case 'linear'
               expr_grass=['((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))'];
            
           case 'gaussian_oriented'
               expr_grass=['Exp(IntCst(2)*c*((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1)-IntCst(1)))'];
            
           case 'binet'
               expr_grass=['Square((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))'];
            
           case 'gaussian_unoriented'
               expr_grass=['Exp(IntCst(2)*c*(Square((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))-IntCst(1)))'];
            
       end
    
       K = Kernel([expr_geom '*' expr_grass '*r'],'a=Pm(0,1)','c=Pm(1,1)',['x=Vx(2,' num2str(eval('d')) ')'],['y=Vy(3,' num2str(eval('d')) ')'],['u1=Vx(4,' num2str(eval('d')) ')'],['u2=Vx(5,' num2str(eval('d')) ')'],['v1=Vy(6,' num2str(eval('d')) ')'],['v2=Vy(7,' num2str(eval('d')) ')'],'r=Vy(8,1)');

       res = sum(K(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XitildeX{1}',XitildeX{2}',XitildeY{1}',XitildeY{2}',norm_XiY').*norm_XiX');
    
    
end

end