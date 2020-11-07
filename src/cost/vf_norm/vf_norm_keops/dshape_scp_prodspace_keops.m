function [g] = dshape_scp_prodspace_keops(X,Y,objfun)
% function [Dcenter_faceXg,DtanvectXg,DsignalXg]=dshape_scp_prodspace_keops(cX,cY,XiX,XiY,kernel_geom,kernel_signal,kernel_grass,kernel_size_geom,kernel_size_signal,kernel_size_grass)

% cX = X.center'; 
% cY = Y.center'; 
% 
% kernel_geom = objfun.kernel_geom;
% kernel_size_geom = objfun.kernel_size_geom;
% kernel_grass = objfun.kernel_grass;
% kernel_size_grass = objfun.kernel_size_grass;
% 
% XiX = X.vector;
% XiY = Y.vector;
% [d,Nx]=size(cX);
% m=length(XiX);

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
       GxK = GradKernel(K,'x','p1=Vx(1)');
       Dcenter_faceXg = GxK(1/kernel_size_geom^2,cX,cY,norm_XiY',ones(1,Nx)).*repmat(norm_XiX',d,1);

       
    case 1  
       %tic
       norm_XiX=sqrt(sum(XiX{1}.^2,2));
       norm_XiY=sqrt(sum(XiY{1}.^2,2));
    
       XitildeX=XiX{1}./repmat(norm_XiX,1,d);
       XitildeY=XiY{1}./repmat(norm_XiY,1,d);
    %toc
    
       switch kernel_geom
           case 'gaussian'
               expr_geom=['Exp(-SqDist(x,y)*a)'];
            
           case 'cauchy'
               expr_geom=['IntCst(1)/(IntCst(1)+SqDist(x,y)*a)'];
            
       end
    
%     switch kernel_signal
%         case 'gaussian'
%             expr_signal=['Exp(-SqDist(f,g)*b)'];
%             
%         case 'cauchy'
%             expr_signal=['IntCst(1)/(IntCst(1)+SqDist(f,g)*b)'];
%             
%     end
    
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
    
       GxK = GradKernel(K,'x','p1=Vx(1)');
       GuK = GradKernel(K,'u','p3=Vx(1)');
    
       Dcenter_faceXg = GxK(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XitildeX',XitildeY',norm_XiY',ones(1,Nx)).*repmat(norm_XiX',d,1);

       Gu=GuK(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XitildeX',XitildeY',norm_XiY',ones(1,Nx));
         

       DtanvectXg{1} = Gu - repmat(sum(Gu.*XitildeX',1),d,1).*XitildeX' +...
           repmat(K(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XitildeX',XitildeY',norm_XiY'),d,1).*XitildeX';
    
       DtanvectXg{1} = DtanvectXg{1}';
    
    case 2
    
       norm_XiX=pvec_Norm(XiX);
       norm_XiY=pvec_Norm(XiY);
    
       XitildeY=XiY;
    
       XitildeY{1}=XiY{1}./repmat(norm_XiY,1,d);
    
       switch kernel_geom
           case 'gaussian'
               expr_geom=['Exp(-SqDist(x,y)*a)'];
                   
           case 'cauchy'
               expr_geom=['IntCst(1)/(IntCst(1)+SqDist(x,y)*a)'];
       end
    
%     switch kernel_signal
%         case 'gaussian'
%             expr_signal=['Exp(-SqDist(f,g)*b)'];
%             
%         case 'cauchy'
%             expr_signal=['IntCst(1)/(IntCst(1)+SqDist(f,g)*b)'];
%             
%     end
    
       switch kernel_grass
           case 'linear'
               expr_grass=['(((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))/r1)'];
            
           case 'gaussian_oriented'
               expr_grass=['Exp(IntCst(2)*c*(((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))/r1-IntCst(1)))'];
            
           case 'binet'
               expr_grass=['(Square((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))/Square(r1))'];
            
           case 'gaussian_unoriented'
               expr_grass=['Exp(IntCst(2)*c*(Square((u1|v1)*(u2|v2)-(u1|v2)*(u2|v1))/Square(r1)-IntCst(1)))'];
            
       end
    
%     K = Kernel([expr_geom '*' expr_signal '*' expr_grass '*r1*r2'],'a=Pm(0,1)','b=Pm(1,1)','c=Pm(2,1)',['x=Vx(3,' num2str(eval('d')) ')'],['y=Vy(4,' num2str(eval('d')) ')'],['u1=Vx(5,' num2str(eval('d+1')) ')'],['u2=Vx(6,' num2str(eval('d+1')) ')'],['v1=Vy(7,' num2str(eval('d+1')) ')'],['v2=Vy(8,' num2str(eval('d+1')) ')'],'r1=Vx(9,1)','r2=Vy(10,1)');
       K = Kernel([expr_geom '*' expr_grass '*r1*r2'],'a=Pm(0,1)','c=Pm(1,1)',['x=Vx(2,' num2str(eval('d')) ')'],['y=Vy(3,' num2str(eval('d')) ')'],['u1=Vx(4,' num2str(eval('d')) ')'],['u2=Vx(5,' num2str(eval('d')) ')'],['v1=Vy(6,' num2str(eval('d')) ')'],['v2=Vy(7,' num2str(eval('d')) ')'],'r1=Vx(8,1)','r2=Vy(9,1)');

       GxK = GradKernel(K,'x','p1=Vx(1)');
%     GfK = GradKernel(K,'f','p2=Vx(1)');
       Gu1K = GradKernel(K,'u1','p3=Vx(1)');
       Gu2K = GradKernel(K,'u2','p3=Vx(1)');
       GrK = GradKernel(K,'r1','p4=Vx(1)');

       Dcenter_faceXg = GxK(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XiX{1}',XiX{2}',XitildeY{1}',XitildeY{2}',norm_XiX',norm_XiY',ones(1,Nx));
%     DsignalXg = GfK(1/kernel_size_geom^2,1/kernel_size_signal^2,1/kernel_size_grass^2,cX,cY,XiX{1}',XiX{2}',XitildeY{1}',XitildeY{2}',norm_XiX',norm_XiY',ones(1,Nx));
    
       Gr=GrK(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XiX{1}',XiX{2}',XitildeY{1}',XitildeY{2}',norm_XiX',norm_XiY',ones(1,Nx));
       NX1=sum(XiX{1}.^2,2);
       NX2=sum(XiX{2}.^2,2);
       InnX=sum(XiX{1}.*XiX{2},2);
    
       Gu=Gu1K(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XiX{1}',XiX{2}',XitildeY{1}',XitildeY{2}',norm_XiX',norm_XiY',ones(1,Nx));
       DtanvectXg{1} = Gu'+repmat((Gr')./norm_XiX,1,d).*(repmat(NX2,1,d).*XiX{1}-repmat(InnX,1,d).*XiX{2});
    
       Gu=Gu2K(1/kernel_size_geom^2,1/kernel_size_grass^2,cX,cY,XiX{1}',XiX{2}',XitildeY{1}',XitildeY{2}',norm_XiX',norm_XiY',ones(1,Nx));
       DtanvectXg{2} = Gu'+repmat((Gr')./norm_XiX,1,d).*(repmat(NX1,1,d).*XiX{2}-repmat(InnX,1,d).*XiX{1});
    
    
end

g.center = Dcenter_faceXg';
if m>0
   g.vector = DtanvectXg;
end
% DsignalXg = DsignalXg';


end